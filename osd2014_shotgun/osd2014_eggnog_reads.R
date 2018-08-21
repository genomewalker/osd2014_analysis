library(DESeq2)
library(RPostgreSQL)
library(phyloseq)
library(metagenomeSeq)
library(ggpubr)
library(tidyverse)

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")


# Create phyloseq objects -------------------------------------------------
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)

osd2014_eggnog_reads <- tbl(my_db, "osd2014_reads_emapper_results") %>%
  collect(n = Inf)  %>%
  filter(label %in% osd2014_amp_mg_intersect$label)
osd2014_eggnog_reads <- osd2014_eggnog_reads %>%
  spread(label, value = abun, fill = 0) %>%
  as.data.frame()
rownames(osd2014_eggnog_reads) <- osd2014_eggnog_reads$group_nam
osd2014_eggnog_reads$group_nam <- NULL

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

osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)
osd2014_cdata <- as.data.frame(osd2014_cdata)
rownames(osd2014_cdata) <- osd2014_cdata$label

osd2014_eggnog_reads_phyloseq <- phyloseq(otu_table(as.matrix(osd2014_eggnog_reads), taxa_are_rows = TRUE), sample_data(osd2014_cdata), tax_table(as.matrix(eggnog_annotation_comb)))

# eggNOG table pre-processing ------------------------------------------------

osd2014_eggnog_reads_phyloseq_df <- psmelt(osd2014_eggnog_reads_phyloseq) %>%
  select(Sample, OTU, Abundance) %>%
  rename(label = Sample, group_nam = OTU, counts = Abundance) %>%
  tbl_df

osd2014_eggnog_summary_tcounts <-  osd2014_eggnog_reads_phyloseq_df %>%
  group_by(group_nam) %>%
  summarise(total_counts = sum(counts), occ = sum(counts > 0)) %>%
  ungroup()

# Filtering for alpha diversity analyses ----------------------------------

# We will remove absolute singletons

osd2014_eggnog_summary_alpha <- osd2014_eggnog_summary_tcounts %>%
  select(group_nam, total_counts, occ) %>%
  filter(total_counts > 1)

osd2014_eggnog_summary_alpha_discarded <- osd2014_eggnog_summary_tcounts %>%
  select(group_nam, total_counts, occ) %>%
  filter(total_counts == 1, occ == 1)

# Summary without the singletons ------------------------------------------

osd2014_eggnog_summary <- osd2014_eggnog_reads_phyloseq_df %>%
  filter(group_nam %in% osd2014_eggnog_summary_alpha$group_nam) %>%
  group_by(label) %>%
  mutate(prop = counts/sum(counts)) %>%
  ungroup() %>%
  group_by(group_nam) %>%
  mutate(mean_prop = sum(prop)/nsamples(osd2014_eggnog_reads_phyloseq)) %>%
  ungroup() %>%
  filter(counts > 0) %>%
  arrange(counts)

osd2014_eggnog_sample_sum_plot <- osd2014_eggnog_summary %>%
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
osd2014_eggnog_stats_filtered_prev <- bind_rows(lapply(X = meanprops, FUN = meanprop_filtering, Y = osd2014_eggnog_summary))
osd2014_eggnog_stats_filtered_prev <- osd2014_eggnog_stats_filtered_prev %>%
  dplyr::select(mean_prop, n_samples, n_clstrs, total_counts, sparsity) %>%
  rename(Samples = n_samples, eggNOGs = n_clstrs, Counts = total_counts, Sparsity = sparsity) %>%
  gather(Cluster_stats, value, -mean_prop) %>%
  mutate(value = ifelse(is.na(value), 0, value))


spar <- osd2014_eggnog_stats_filtered_prev %>% filter(Cluster_stats == "Sparsity")
f_spar <- approxfun(spar$mean_prop, spar$value)

samp <- osd2014_eggnog_stats_filtered_prev %>% filter(Cluster_stats == "Samples")
f_samp <- approxfun(samp$mean_prop, samp$value)

counts <- osd2014_eggnog_stats_filtered_prev %>% filter(Cluster_stats == "Counts")
f_counts <- approxfun(counts$mean_prop, counts$value)

clusters <- osd2014_eggnog_stats_filtered_prev %>% filter(Cluster_stats == "eggNOGs")
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
  ylab("ASVs") +
  theme_light() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        plot.title = element_text(size = 8, face = "bold", vjust = -1)) +
  ggtitle(paste("ASVs:", scales::comma(round(f_clusters(p),3))))

osd2014_eggnog_stats_filtered_prev_plot <- ggarrange(p_spar_1, p_counts_1, p_clusters_1, p_spar, p_counts, p_clusters, labels = "AUTO", nrow = 2, ncol = 3)
ggsave(osd2014_eggnog_stats_filtered_prev_plot, filename = "osd2014_shotgun/figures/osd2014_filt_mean_prop.pdf", width = 11.69, height = 8.27)


# Distribution of reads per sample ----------------------------------------
osd2014_sample_counts_plot <- bind_rows(sample_sums(osd2014_eggnog_reads_phyloseq) %>%
                                          map_df(~ data_frame(counts = .x), .id = "label") %>%
                                          mutate(class = "original"),
                                        osd2014_eggnog_summary %>%
                                          filter(mean_prop >= 1e-5) %>%
                                          group_by(label) %>%
                                          summarise(counts = sum(counts)) %>%
                                          mutate(class = "beta"),
                                        osd2014_eggnog_summary %>%
                                          filter(asv %in% osd2014_eggnog_summary_alpha$asv) %>%
                                          group_by(label) %>%
                                          summarise(counts = sum(counts)) %>%
                                          mutate(class = "alpha")) %>%
  mutate(class = fct_relevel(class, c("original", "alpha", "beta")), label = fct_reorder(label, -counts)) %>%
  ggplot(aes(label, counts, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_light() +
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")

ggsave(osd2014_sample_counts_plot, filename = "osd2014_eggnog_asv/figures/osd2014_sample_counts_plot.pdf", width = 11.69, height = 8.27)



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
osd2014_eggnog_phyloseq_alpha <- prune_taxa(osd2014_eggnog_summary_alpha$group_nam, osd2014_eggnog_reads_phyloseq)

osd2014_eggnog_phyloseq_alpha_css <- cssTrans(osd2014_eggnog_phyloseq_alpha, norm = T, log = T)
osd2014_eggnog_phyloseq_beta_css <- prune_taxa(osd2014_eggnog_summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, osd2014_eggnog_phyloseq_alpha_css)

nlib <- plyr::round_any(min(sample_sums(osd2014_eggnog_phyloseq_alpha)), 1000, f = floor)
osd2014_eggnog_phyloseq_alpha_scaled <- scale_reads(osd2014_eggnog_phyloseq_alpha, n = nlib)
osd2014_eggnog_phyloseq_beta_scaled <- prune_taxa(osd2014_eggnog_summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, osd2014_eggnog_phyloseq_alpha_scaled)

osd2014_eggnog_phyloseq_beta <- prune_taxa(osd2014_eggnog_summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, osd2014_eggnog_phyloseq_alpha)

# DESEQ2 VST and normalised counts
osd2014_eggnog_phyloseq_alpha_deseq <- phyloseq_to_deseq2(osd2014_eggnog_phyloseq_alpha, ~1)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(osd2014_eggnog_phyloseq_alpha_deseq), 1, gm_mean)
diagdds <- estimateSizeFactors(osd2014_eggnog_phyloseq_alpha_deseq, type = "poscounts")
diagdds <- estimateDispersions(diagdds)
diagvst <- varianceStabilizingTransformation(diagdds, blind = FALSE)
diagvst <- assay(diagvst)
diagvst[diagvst < 0] <- 0
osd2014_eggnog_phyloseq_alpha_vst <- osd2014_eggnog_phyloseq_alpha
osd2014_eggnog_phyloseq_alpha_norm <- osd2014_eggnog_phyloseq_alpha
otu_table(osd2014_eggnog_phyloseq_alpha_vst) <- otu_table(diagvst, taxa_are_rows = TRUE)
otu_table(osd2014_eggnog_phyloseq_alpha_norm) <- otu_table(counts(diagdds, normalized = TRUE), taxa_are_rows = TRUE)

osd2014_eggnog_phyloseq_beta_vst <- prune_taxa(osd2014_eggnog_summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, osd2014_eggnog_phyloseq_alpha_vst)
osd2014_eggnog_phyloseq_beta_norm <- prune_taxa(osd2014_eggnog_summary %>% filter(mean_prop >= 1e-5) %>% .$group_nam, osd2014_eggnog_phyloseq_alpha_norm)



save.image(file = "osd2014_shotgun/data/osd2014_eggnog_preprocess.Rdata", compress = TRUE)
save(osd2014_eggnog_phyloseq_alpha, osd2014_eggnog_phyloseq_alpha_css, osd2014_eggnog_phyloseq_beta_css,
     osd2014_eggnog_phyloseq_alpha_scaled, osd2014_eggnog_phyloseq_beta_scaled, osd2014_eggnog_phyloseq_beta,
     osd2014_eggnog_phyloseq_alpha_norm, osd2014_eggnog_phyloseq_alpha_vst,
     osd2014_eggnog_phyloseq_beta_vst, osd2014_eggnog_phyloseq_beta_norm, file = "osd2014_shotgun/data/osd2014_eggnog_physeq_filt_objects.Rdata")
