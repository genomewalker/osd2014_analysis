library(phyloseq)
library(vegan)
library(tidyverse)
library(RPostgreSQL)

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

tara_emapper_rand_summary <- tbl(my_db, "tara_emapper_rand_summary") %>%
  collect(n = Inf) %>%
  select(label, group_nam, mean) %>%
  rename(abundance = mean) %>%
  mutate(group_nam = ifelse(grepl("^COG", group_nam), group_nam, paste0("ENOG41", group_nam)))

tara_emapper_assm <- tbl(my_db, "tara_emapper_assm") %>%
  collect(n = Inf)

myround <- function(x) { trunc(x + 0.5) }

tara_cdata <- tbl(my_db, "tara_metadata_prok") %>%
  collect(n = Inf) %>%
  filter(grepl("SRF", label))
tara_cdata <- as.data.frame(tara_cdata)
rownames(tara_cdata) <- tara_cdata$label


tara_emapper_rand_summary_wide <- tara_emapper_rand_summary %>%
  spread(label, abundance, fill = 0) %>%
  as.data.frame()
rownames(tara_emapper_rand_summary_wide) <- tara_emapper_rand_summary_wide$group_nam
tara_emapper_rand_summary_wide$group_nam <- NULL


tara_emapper_assm_wide <- tara_emapper_assm %>%
  filter(label %in% rownames(tara_cdata)) %>%
  spread(label, abundance, fill = 0) %>%
  as.data.frame()
rownames(tara_emapper_assm_wide) <- tara_emapper_assm_wide$group_nam
tara_emapper_assm_wide$group_nam <- NULL

eggnog_annotation <- tbl(my_db, "eggnog4_annotation") %>%
  collect(n = Inf)

eggnog_funcat <- tbl(my_db, "eggnog4_funcat") %>%
  collect(n = Inf)

eggnog_annotation_comb <- eggnog_annotation %>%
  left_join(eggnog_funcat) %>%
  group_by(taxonomic_level, group_nam, protein_count, species_count, consensus_functional_description) %>%
  summarise(cog_functional_category = paste(cog_functional_category, collapse = "|"),
            cog_description = paste(cog_description, collapse = "|"),
            cog_category = paste(cog_category, collapse = "|")) %>%
  as.data.frame()
rownames(eggnog_annotation_comb) <- eggnog_annotation_comb$group_nam


tara_emapper_rand_phyloseq <- phyloseq(otu_table(as.matrix(myround(tara_emapper_rand_summary_wide)), taxa_are_rows = TRUE), sample_data(tara_cdata), tax_table(as.matrix(eggnog_annotation_comb)))
tara_emapper_assm_phyloseq <- phyloseq(otu_table(as.matrix(myround(tara_emapper_assm_wide)), taxa_are_rows = TRUE), sample_data(tara_cdata), tax_table(as.matrix(eggnog_annotation_comb)))

preprocess <- function(X) {

  n_samples <- X$label %>% unique() %>% length()

  summary_tcounts <- X %>%
    rename(counts = abundance) %>%
    mutate(counts = myround(counts)) %>%
    group_by(group_nam) %>%
    summarise(total_counts = sum(counts), occ = sum(counts > 0)) %>%
    ungroup()

  # Filtering for alpha diversity analyses ----------------------------------

  # We will remove absolute singletons

  summary_alpha <- summary_tcounts %>%
    select(group_nam, total_counts, occ) %>%
    filter(total_counts > 1)

  summary_alpha_discarded <- summary_tcounts %>%
    select(group_nam, total_counts, occ) %>%
    filter(total_counts == 1, occ == 1)

  # Summary without the singletons ------------------------------------------

  summary <- X %>%
    rename(counts = abundance) %>%
    mutate(counts = myround(counts)) %>%
    filter(group_nam %in% summary_alpha$group_nam) %>%
    group_by(label) %>%
    mutate(prop = counts/sum(counts)) %>%
    ungroup() %>%
    group_by(group_nam) %>%
    mutate(mean_prop = sum(prop)/n_samples) %>%
    ungroup() %>%
    filter(counts > 0) %>%
    arrange(counts)

  sample_sum_plot <- summary %>%
    group_by(label) %>%
    summarise(counts = sum(counts), group_nam = length(unique(group_nam))) %>%
    mutate(label = fct_reorder(label, -counts)) %>%
    ggplot(aes(counts, group_nam)) +
    geom_point() +
    theme_bw() +
    scale_y_continuous(labels = scales::comma) +
    scale_x_continuous(labels = scales::comma) +
    xlab("Counts") +
    ylab("Number of eggNOGs")

  results <- list(summary = summary, plot = sample_sum_plot)

}

tara_emapper_rand_phyloseq_preprop <- preprocess(tara_emapper_rand_summary)
tara_emapper_assm_phyloseq_preprop <- preprocess(tara_emapper_assm)


# We will filter low abundant eggNOGs for beta diversity analyses ------------

meanprop_filtering <- function(X, Y) {
  df <- as_tibble(Y) %>%
    ungroup() %>%
    dplyr::select(label, group_nam, counts, mean_prop) %>%
    filter(X <= mean_prop) %>%
    summarise(mean_prop = X,
              n_clstrs = length(unique(group_nam)),
              n_samples = length(unique(label)),
              no_zero = sum(counts > 0),
              sparsity = 1 - (no_zero/(n_clstrs * n_samples)),
              total_counts = sum(counts))
  return(df)
}

meanprops <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)

# apply function
tara_emapper_rand_stats_filtered_prev <- bind_rows(lapply(X = meanprops, FUN = meanprop_filtering, Y = tara_emapper_rand_phyloseq_preprop$summary))
tara_emapper_rand_stats_filtered_prev <- tara_emapper_rand_stats_filtered_prev %>%
  dplyr::select(mean_prop, n_samples, n_clstrs, total_counts, sparsity) %>%
  rename(Samples = n_samples, eggNOGs = n_clstrs, Counts = total_counts, Sparsity = sparsity) %>%
  gather(Cluster_stats, value, -mean_prop) %>%
  mutate(value = ifelse(is.na(value), 0, value))

tara_emapper_assm_stats_filtered_prev <- bind_rows(lapply(X = meanprops, FUN = meanprop_filtering, Y = tara_emapper_assm_phyloseq_preprop$summary))
tara_emapper_assm_stats_filtered_prev <- tara_emapper_assm_stats_filtered_prev %>%
  dplyr::select(mean_prop, n_samples, n_clstrs, total_counts, sparsity) %>%
  rename(Samples = n_samples, eggNOGs = n_clstrs, Counts = total_counts, Sparsity = sparsity) %>%
  gather(Cluster_stats, value, -mean_prop) %>%
  mutate(value = ifelse(is.na(value), 0, value))


# Filtering trersholds for the RAND emapper results ------------------------

spar <- tara_emapper_rand_stats_filtered_prev %>% filter(Cluster_stats == "Sparsity")
f_spar <- approxfun(spar$mean_prop, spar$value)

samp <- tara_emapper_rand_stats_filtered_prev %>% filter(Cluster_stats == "Samples")
f_samp <- approxfun(samp$mean_prop, samp$value)

counts <- tara_emapper_rand_stats_filtered_prev %>% filter(Cluster_stats == "Counts")
f_counts <- approxfun(counts$mean_prop, counts$value)

clusters <- tara_emapper_rand_stats_filtered_prev %>% filter(Cluster_stats == "eggNOGs")
f_clusters <- approxfun(clusters$mean_prop, clusters$value)

p <- 1e-5
p1 <- 1e-6
p_spar_1 <- ggplot(spar, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p1, y = 0, xend = p1, yend = f_spar(p1)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_spar(p1), xend = p1, yend = f_spar(p1)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_continuous(labels = scales::percent) +
  geom_point(aes(p1, f_spar(p1)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("Sparsity") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("Sparsity:", scales::percent(round(f_spar(p1),3))))

p_counts_1 <- ggplot(counts, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p1, y = 0, xend = p1, yend = f_counts(p1)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_counts(p1), xend = p1, yend = f_counts(p1)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_continuous(labels=function(x){paste0(x/1e6,"M")}) +
  geom_point(aes(p1, f_counts(p1)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("Counts") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("Counts:", scales::comma(round(f_counts(p1),3))))

p_clusters_1 <- ggplot(clusters, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p1, y = 0, xend = p1, yend = f_clusters(p1)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_clusters(p1), xend = p1, yend = f_clusters(p1)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_log10(labels=scales::comma) +
  geom_point(aes(p1, f_clusters(p1)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("eggNOGs") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("eggNOGs:", scales::comma(round(f_clusters(p1),3))))

p_spar <- ggplot(spar, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p, y = 0, xend = p, yend = f_spar(p)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_spar(p), xend = p, yend = f_spar(p)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_continuous(labels = scales::percent) +
  geom_point(aes(p, f_spar(p)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("Sparsity") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("Sparsity:", scales::percent(round(f_spar(p),3))))

p_counts <- ggplot(counts, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p, y = 0, xend = p, yend = f_counts(p)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_counts(p), xend = p, yend = f_counts(p)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_continuous(labels=function(x){paste0(x/1e6,"M")}) +
  geom_point(aes(p, f_counts(p)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("Counts") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("Counts:", scales::comma(round(f_counts(p),3))))

p_clusters <- ggplot(clusters, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p, y = 0, xend = p, yend = f_clusters(p)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_clusters(p), xend = p, yend = f_clusters(p)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_log10(labels=scales::comma) +
  geom_point(aes(p, f_clusters(p)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("eggNOG") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("eggNOGss:", scales::comma(round(f_clusters(p),3))))

tara_emapper_rand_stats_filtered_prev_plot <- ggarrange(p_spar_1, p_counts_1, p_clusters_1, p_spar, p_counts, p_clusters, labels = "AUTO", nrow = 2, ncol = 3)
ggsave(tara_emapper_rand_stats_filtered_prev_plot, filename = "osd2014_shotgun/figures/tara_emapper_rand_filt_mean_prop.pdf", width = 11.69, height = 8.27)


# Filtering trersholds for the ASSM emapper results ------------------------

spar <- tara_emapper_assm_stats_filtered_prev %>% filter(Cluster_stats == "Sparsity")
f_spar <- approxfun(spar$mean_prop, spar$value)

samp <- tara_emapper_assm_stats_filtered_prev %>% filter(Cluster_stats == "Samples")
f_samp <- approxfun(samp$mean_prop, samp$value)

counts <- tara_emapper_assm_stats_filtered_prev %>% filter(Cluster_stats == "Counts")
f_counts <- approxfun(counts$mean_prop, counts$value)

clusters <- tara_emapper_assm_stats_filtered_prev %>% filter(Cluster_stats == "eggNOGs")
f_clusters <- approxfun(clusters$mean_prop, clusters$value)

p <- 1e-5
p1 <- 1e-6
p_spar_1 <- ggplot(spar, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p1, y = 0, xend = p1, yend = f_spar(p1)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_spar(p1), xend = p1, yend = f_spar(p1)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_continuous(labels = scales::percent) +
  geom_point(aes(p1, f_spar(p1)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("Sparsity") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("Sparsity:", scales::percent(round(f_spar(p1),3))))

p_counts_1 <- ggplot(counts, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p1, y = 0, xend = p1, yend = f_counts(p1)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_counts(p1), xend = p1, yend = f_counts(p1)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_continuous(labels=function(x){paste0(x/1e6,"M")}) +
  geom_point(aes(p1, f_counts(p1)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("Counts") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("Counts:", scales::comma(round(f_counts(p1),3))))

p_clusters_1 <- ggplot(clusters, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p1, y = 0, xend = p1, yend = f_clusters(p1)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_clusters(p1), xend = p1, yend = f_clusters(p1)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_log10(labels=scales::comma) +
  geom_point(aes(p1, f_clusters(p1)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("eggNOGs") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("eggNOGs:", scales::comma(round(f_clusters(p1),3))))

p_spar <- ggplot(spar, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p, y = 0, xend = p, yend = f_spar(p)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_spar(p), xend = p, yend = f_spar(p)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_continuous(labels = scales::percent) +
  geom_point(aes(p, f_spar(p)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("Sparsity") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("Sparsity:", scales::percent(round(f_spar(p),3))))

p_counts <- ggplot(counts, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p, y = 0, xend = p, yend = f_counts(p)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_counts(p), xend = p, yend = f_counts(p)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_continuous(labels=function(x){paste0(x/1e6,"M")}) +
  geom_point(aes(p, f_counts(p)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("Counts") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("Counts:", scales::comma(round(f_counts(p),3))))

p_clusters <- ggplot(clusters, aes(mean_prop, value)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = p, y = 0, xend = p, yend = f_clusters(p)), lty = "dotted") +
  geom_segment(aes(x = 0, y = f_clusters(p), xend = p, yend = f_clusters(p)), lty = "dotted") +
  scale_x_log10(breaks = meanprops) +
  scale_y_log10(labels=scales::comma) +
  geom_point(aes(p, f_clusters(p)), size = 2, shape = 21, fill = "white") +
  xlab("Mean proportion threshold") +
  ylab("eggNOG") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("eggNOGss:", scales::comma(round(f_clusters(p),3))))

tara_emapper_assm_stats_filtered_prev_plot <- ggarrange(p_spar_1, p_counts_1, p_clusters_1, p_spar, p_counts, p_clusters, labels = "AUTO", nrow = 2, ncol = 3)
ggsave(tara_emapper_assm_stats_filtered_prev_plot, filename = "osd2014_shotgun/figures/tara_emapper_assm_filt_mean_prop.pdf", width = 11.69, height = 8.27)




# Distribution of reads per sample ----------------------------------------
tara_emapper_rand_sample_counts_plot <- bind_rows(sample_sums(tara_emapper_rand_phyloseq) %>%
                                          map_df(~ data_frame(counts = .x), .id = "label") %>%
                                          mutate(class = "original"),
                                        tara_emapper_rand_phyloseq_preprop$summary %>%
                                          filter(mean_prop >= 1e-5) %>%
                                          group_by(label) %>%
                                          summarise(counts = sum(counts)) %>%
                                          mutate(class = "beta"),
                                        tara_emapper_rand_summary %>%
                                          filter(group_nam %in% tara_emapper_rand_phyloseq_preprop$summary$group_nam) %>%
                                          group_by(label) %>%
                                          summarise(counts = sum(abundance)) %>%
                                          mutate(class = "alpha")) %>%
  mutate(class = fct_relevel(class, c("original", "alpha", "beta")), label = fct_reorder(label, -counts)) %>%
  ggplot(aes(label, counts, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_light() +
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")

ggsave(tara_emapper_rand_sample_counts_plot, filename = "osd2014_shotgun/figures/tara_emapper_rand_sample_counts_plot.pdf", width = 11.69, height = 8.27)

tara_emapper_assm_sample_counts_plot <- bind_rows(sample_sums(tara_emapper_assm_phyloseq) %>%
                                                    map_df(~ data_frame(counts = .x), .id = "label") %>%
                                                    mutate(class = "original"),
                                                  tara_emapper_assm_phyloseq_preprop$summary %>%
                                                    filter(mean_prop >= 1e-5) %>%
                                                    group_by(label) %>%
                                                    summarise(counts = sum(counts)) %>%
                                                    mutate(class = "beta"),
                                                  tara_emapper_assm %>%
                                                    filter(group_nam %in% tara_emapper_assm_phyloseq_preprop$summary$group_nam) %>%
                                                    group_by(label) %>%
                                                    summarise(counts = sum(abundance)) %>%
                                                    mutate(class = "alpha")) %>%
  mutate(class = fct_relevel(class, c("original", "alpha", "beta")), label = fct_reorder(label, -counts)) %>%
  ggplot(aes(label, counts, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_light() +
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")

ggsave(tara_emapper_assm_sample_counts_plot, filename = "osd2014_shotgun/figures/tara_emapper_assm_sample_counts_plot.pdf", width = 11.69, height = 8.27)



# Transform data with CSS -------------------------------------------------

# Rows should correspond to features and columns to samples.
cssTrans <- function(f.physeq.p = f.physeq.p, norm = norm, log = log){

  if (taxa_are_rows(f.physeq.p)) {
    f.physeq.p <- (f.physeq.p)
  }else{
    f.physeq.p <- t(f.physeq.p)
  }

  OTU <- as((otu_table(f.physeq.p, taxa_are_rows = TRUE)), "matrix")
  MGS <- newMRexperiment(
    counts = (OTU)
  )
  MGS <- cumNorm(MGS, p = cumNormStat(MGS))
  f.norm.p <- f.physeq.p
  otu_table(f.norm.p) <- otu_table((as.matrix(MRcounts(
    MGS,
    norm = norm,
    log = log,
    sl = median(unlist(normFactors(MGS)))
  ))), taxa_are_rows = T)
  return(f.norm.p)
}

# From https://github.com/DenefLab/MicrobeMiseq/blob/master/R/miseqR.R
# Better rounding function than R's base round
myround <- function(x) { trunc(x + 0.5) }

# Scales reads by
# 1) taking proportions
# 2) multiplying by a given library size of n
# 3) rounding
# Default for n is the minimum sample size in your library
# Default for round is floor
scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "floor") {

  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq,
                                          function(x) {(n * x/sum(x))}
  )

  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- myround(otu_table(physeq.scale))
  }

  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}



# Create phyloseq objects for alpha and beta
tara_emapper_assm_phyloseq_alpha <- prune_taxa(tara_emapper_assm_phyloseq_preprop$summary$group_nam, tara_emapper_assm_phyloseq)
tara_emapper_assm_phyloseq_alpha_css <- cssTrans(tara_emapper_assm_phyloseq_alpha, norm = T, log = T)
tara_emapper_assm_phyloseq_beta_css <- prune_taxa(tara_emapper_assm_phyloseq_preprop$summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, tara_emapper_assm_phyloseq_alpha_css)

nlib <- plyr::round_any(min(sample_sums(tara_emapper_assm_phyloseq_alpha)), 1000, f = floor)
tara_emapper_assm_phyloseq_alpha_scaled <- scale_reads(tara_emapper_assm_phyloseq_alpha, n = nlib)
tara_emapper_assm_phyloseq_beta_scaled <- prune_taxa(tara_emapper_assm_phyloseq_preprop$summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, tara_emapper_assm_phyloseq_alpha_scaled)

tara_emapper_assm_phyloseq_beta <- prune_taxa(tara_emapper_assm_phyloseq_preprop$summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, tara_emapper_assm_phyloseq_alpha)

# DESEQ2 VST and normalised counts
tara_emapper_assm_phyloseq_alpha_deseq <- phyloseq_to_deseq2(tara_emapper_assm_phyloseq_alpha, ~1)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(tara_emapper_assm_phyloseq_alpha_deseq), 1, gm_mean)
diagdds <- estimateSizeFactors(tara_emapper_assm_phyloseq_alpha_deseq, type = "poscounts")
diagdds <- estimateDispersions(diagdds)
diagvst <- varianceStabilizingTransformation(diagdds, blind = FALSE)
diagvst <- assay(diagvst)
diagvst[diagvst < 0] <- 0
tara_emapper_assm_phyloseq_alpha_vst <- tara_emapper_assm_phyloseq_alpha
tara_emapper_assm_phyloseq_alpha_norm <- tara_emapper_assm_phyloseq_alpha
otu_table(tara_emapper_assm_phyloseq_alpha_vst) <- otu_table(diagvst, taxa_are_rows = TRUE)
otu_table(tara_emapper_assm_phyloseq_alpha_norm) <- otu_table(counts(diagdds, normalized = TRUE), taxa_are_rows = TRUE)

tara_emapper_assm_phyloseq_beta_vst <- prune_taxa(tara_emapper_assm_phyloseq_preprop$summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, tara_emapper_assm_phyloseq_alpha_vst)
tara_emapper_assm_phyloseq_beta_norm <- prune_taxa(tara_emapper_assm_phyloseq_preprop$summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, tara_emapper_assm_phyloseq_alpha_norm)



# Create phyloseq objects for alpha and beta
tara_emapper_rand_phyloseq_alpha <- prune_taxa(tara_emapper_rand_phyloseq_preprop$summary$group_nam, tara_emapper_rand_phyloseq)
tara_emapper_rand_phyloseq_alpha_css <- cssTrans(tara_emapper_rand_phyloseq_alpha, norm = T, log = T)
tara_emapper_rand_phyloseq_beta_css <- prune_taxa(tara_emapper_rand_phyloseq_preprop$summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, tara_emapper_rand_phyloseq_alpha_css)

nlib <- plyr::round_any(min(sample_sums(tara_emapper_rand_phyloseq_alpha)), 1000, f = floor)
tara_emapper_rand_phyloseq_alpha_scaled <- scale_reads(tara_emapper_rand_phyloseq_alpha, n = nlib)
tara_emapper_rand_phyloseq_beta_scaled <- prune_taxa(tara_emapper_rand_phyloseq_preprop$summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, tara_emapper_rand_phyloseq_alpha_scaled)

tara_emapper_rand_phyloseq_beta <- prune_taxa(tara_emapper_rand_phyloseq_preprop$summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, tara_emapper_rand_phyloseq_alpha)

# DESEQ2 VST and normalised counts
tara_emapper_rand_phyloseq_alpha_deseq <- phyloseq_to_deseq2(tara_emapper_rand_phyloseq_alpha, ~1)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(tara_emapper_rand_phyloseq_alpha_deseq), 1, gm_mean)
diagdds <- estimateSizeFactors(tara_emapper_rand_phyloseq_alpha_deseq, type = "poscounts")
diagdds <- estimateDispersions(diagdds)
diagvst <- varianceStabilizingTransformation(diagdds, blind = FALSE)
diagvst <- assay(diagvst)
diagvst[diagvst < 0] <- 0
tara_emapper_rand_phyloseq_alpha_vst <- tara_emapper_rand_phyloseq_alpha
tara_emapper_rand_phyloseq_alpha_norm <- tara_emapper_rand_phyloseq_alpha
otu_table(tara_emapper_rand_phyloseq_alpha_vst) <- otu_table(diagvst, taxa_are_rows = TRUE)
otu_table(tara_emapper_rand_phyloseq_alpha_norm) <- otu_table(counts(diagdds, normalized = TRUE), taxa_are_rows = TRUE)

tara_emapper_rand_phyloseq_beta_vst <- prune_taxa(tara_emapper_rand_phyloseq_preprop$summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, tara_emapper_rand_phyloseq_alpha_vst)
tara_emapper_rand_phyloseq_beta_norm <- prune_taxa(tara_emapper_rand_phyloseq_preprop$summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, tara_emapper_rand_phyloseq_alpha_norm)



as(otu_table(tara_emapper_rand_phyloseq_beta_vst), "matrix") %>%
  as_tibble(rownames = "group_nam") %>%
  gather(label, abundance_rand, -group_nam) %>%
  inner_join(as(otu_table(tara_emapper_assm_phyloseq_beta_vst), "matrix") %>%
               as_tibble(rownames = "group_nam") %>%
               gather(label, abundance_assm, -group_nam)) %>%
  ggplot(aes(abundance_rand, abundance_assm)) +
  geom_hex() +
  xlim(c(0,15)) + ylim(c(0,20)) +
  viridis::scale_fill_viridis() +
  geom_abline(intercept = 0, slope = 1) +
  theme_light() +
  xlab("eggNOG in subsampled reads (VST)") +
  ylab("eggNOG in assemblies (VST)")


tara_assm_rand_upset <- as(otu_table(tara_emapper_rand_phyloseq_beta_vst), "matrix") %>%
  as_tibble(rownames = "group_nam") %>%
  gather(label, abundance_rand, -group_nam) %>%
  full_join(as(otu_table(tara_emapper_assm_phyloseq_beta_vst), "matrix") %>%
               as_tibble(rownames = "group_nam") %>%
               gather(label, abundance_assm, -group_nam)) %>%
  mutate(assm = ifelse(is.na(abundance_assm), 0, 1),
         rand = ifelse(is.na(abundance_rand), 0, 1)) %>%
  select(group_nam, assm, rand) %>%
  group_by(group_nam) %>%
  summarise(assm = sum(assm), rand = sum(rand)) %>%
  mutate(assm = ifelse(assm > 0, 1, 0),
         rand = ifelse(rand > 0, 1, 0)) %>%
  dplyr::rename(Assembly = assm, Subsampled = rand) %>%
  gather(class, value, -group_nam) %>%
  spread(class, value, fill = 0) %>%
  as.data.frame() %>% column_to_rownames("group_nam") #%>% t() %>% as.data.frame()
  UpSetR::upset(test, order.by = "freq")


save(tara_emapper_rand_phyloseq_alpha, tara_emapper_rand_phyloseq_alpha_css, tara_emapper_rand_phyloseq_beta_css,
     tara_emapper_rand_phyloseq_alpha_scaled, tara_emapper_rand_phyloseq_beta_scaled, tara_emapper_rand_phyloseq_beta,
     tara_emapper_rand_phyloseq_alpha_norm, tara_emapper_rand_phyloseq_alpha_vst,
     tara_emapper_rand_phyloseq_beta_vst, tara_emapper_rand_phyloseq_beta_norm,
     tara_emapper_assm_phyloseq_alpha, tara_emapper_assm_phyloseq_alpha_css, tara_emapper_assm_phyloseq_beta_css,
     tara_emapper_assm_phyloseq_alpha_scaled, tara_emapper_assm_phyloseq_beta_scaled, tara_emapper_assm_phyloseq_beta,
     tara_emapper_assm_phyloseq_alpha_norm, tara_emapper_assm_phyloseq_alpha_vst,
     tara_emapper_assm_phyloseq_beta_vst, tara_emapper_assm_phyloseq_beta_norm,
     file = "osd2014_shotgun/data/tara_emapper_rand_physeq_filt_objects.Rdata")

# library(RPostgreSQL)  # loads the PostgreSQL driver
# library(tidyverse)
# load("~/Downloads/tara_osd2014_rand_eggnog.Rda")
#
# drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
# con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
# dbWriteTable(con, c("osd_analysis", "tara_emapper_assm"), value=anot_nog_dt_filt,overwrite=TRUE,row.names=FALSE)
