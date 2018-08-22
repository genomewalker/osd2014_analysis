library(phyloseq)
library(vegan)
library(tidyverse)
library(breakaway)
library(RPostgreSQL)
library(pbmcapply)
library(magrittr)
library(ggpubr)
library(pbmcapply)

# BEGIN: WARNING!!!! -------------------------------------------------------------
# You can access to the data used in this analysis in several ways:
# 1. You have a copy of the PostgreSQL DB
# 2. You downloaded the .Rdata files from http://osd2014.metagenomics.eu/ and placed them
#    in the data folder
# 3. You can load the files remotely, it might take a while when the file is very large
# END: WARNING!!!! -------------------------------------------------------------


# BEGIN: WARNING!!: This will load all the data and results for the analysis --------
# Uncomment if you want to use it. Some of the analysis step might require long
# computational times and you might want to use a computer with many cores/CPUs

# load("osd2014_16S_asv/data/osd2014_16S_alpha_diversity.Rdata", verbose = TRUE)
# load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_alpha_diversity.Rdata"), verbose = TRUE)

# END: WARNING!! ---------------------------------------------------------------



# BEGIN: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------

# Load necessary data -----------------------------------------------------
# Use if you have the postgres DB in place

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
st_100_order_terrestrial <- tbl(my_db, "osd2014_st_order_coastal") %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label)
osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)
osd2014_cdata_filt <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label, meow_province %in% osd2014_meow_regions$meow_province)
osd2014_sample_cohesion <- tbl(my_db, "osd2014_sample_cohesion") %>%
  collect(n = Inf)
st_100 <- tbl(my_db, "osd2014_st_100") %>%
  collect(n = Inf)
st_100_long <- st_100 %>%
  gather(variable, value, -label) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  mutate(class = ifelse(variable == "saline_water_gt_25km", "marine", "other")) %>%
  dplyr::select(label, class, value) %>%
  group_by(label, class) %>%
  summarise(prop=sum(value)) %>%
  filter(class == "marine") %>%
  arrange(label)


# If downloaded file at osd2014_16S_asv/data/ use:
load(file = "osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata", verbose = TRUE)
load("osd2014_16S_asv/data/osd2014_16S_asv_divnet_results.Rdata", verbose = TRUE)

# Basic contextual data
load("osd2014_16S_asv/data/osd2014_basic_cdata.Rdata", verbose = TRUE)


# If remote use
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_asv_divnet_results.Rdata"), verbose = TRUE)

# Basic contextual data
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_basic_cdata.Rdata"), verbose = TRUE)
# Load necessary data -----------------------------------------------------

# END: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------



get_div <- function(i, physeq){
  # Subsample
  rarefied_physeq <- rarefy_even_depth(physeq, sample.size = min_lib, verbose = FALSE, replace = TRUE)

  # Calculate observed richness for that group and store value in a df column
  richness <- estimate_richness(rarefied_physeq, measures = "Observed")[ ,1]
  # Calculate Simpson's E for that group and store value in a df column
  simpson <- (estimate_richness(rarefied_physeq, measures = "InvSimpson")[ ,1])
  simpson_e <- estimate_richness(rarefied_physeq, measures = "SimpsonE")[ ,1]
  shannon <- estimate_richness(rarefied_physeq, measures = "Shannon")[ ,1]
  pielous <- estimate_richness(rarefied_physeq, measures = "Pielou")[ ,1]
  df <- data.frame(label = sample_names(rarefied_physeq), iter = rep(i, nsamples(rarefied_physeq)), richness = richness, invsimpson = simpson, shannon = shannon, pielou = pielous, simpson_e = simpson_e)
  return(df)
}


trials <- 1000


# traditional approach ----------------------------------------------------
physeq_filt_prev <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_alpha) > 0, osd2014_dada2_phyloseq_alpha)
min_lib <- min(sample_sums(physeq_filt_prev))
osd2014_asv_alpha <- pbmclapply(1:trials, get_div, physeq = physeq_filt_prev, mc.cores = 4)


lcb <- function(x){quantile(x, probs = 0.025)}
ucb <- function(x){quantile(x, probs = 0.975)}

osd2014_asv_alpha_summary <- bind_rows(osd2014_asv_alpha) %>%
  select(label, richness, invsimpson, shannon, pielou, simpson_e) %>%
  group_by(label) %>%
  summarise_all(funs(mean, sd, lcb, ucb))

get_alpha <- function(X, ns = ns) {
  ns <- sum(X[,1]*X[,2])
  samples <- replicate(1000, resample_estimate(X, my_function = shannon, my_sample_size = ns))
  estimates <- data_frame(est = mean(samples), seest = sd(samples),
                          lcb = quantile(samples, 0.025), ucb = quantile(samples, 0.975), method = "shannon")
  return(estimates)
}
osd2014_frequencytablelist_asv <- build_frequency_count_tables(t(otu_table(osd2014_dada2_phyloseq_alpha)))
alpha_ests <- pbmclapply(X = osd2014_frequencytablelist_asv, get_alpha)
osd2014_asv_alpha_summary_bw <- bind_rows(alpha_ests, .id = "label")


osd2014_asv_obsrich <- osd2014_asv_alpha_summary %>%
  select(label, contains("richness")) %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label)) %>%
  ggplot(aes(label, richness_mean,  ymin = richness_mean - richness_sd, ymax = richness_mean + richness_sd, group = 1)) +
  geom_linerange() +
  geom_line(size = 0.1) +
  geom_point(size = 0.3, color = "red") +
  theme_light() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top") +
  xlab("") +
  ylab("Observed ASVs")



osd2014_asv_shannon <- osd2014_asv_alpha_summary %>%
  select(label, contains("invsimpson")) %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label)) %>%
  ggplot(aes(label, invsimpson_mean,  ymin = invsimpson_mean - invsimpson_sd, ymax = invsimpson_mean + invsimpson_sd, group = 1)) +
  geom_linerange() +
  geom_line(size = 0.1) +
  geom_point(size = 0.3, color = "red") +
  theme_light() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top") +
  xlab("") +
  ylab("Inv Simpson")


osd2014_asv_invsimpson <- osd2014_asv_alpha_summary %>%
  select(label, contains("invsimpson")) %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label)) %>%
  ggplot(aes(label, invsimpson_mean,  ymin = invsimpson_mean - invsimpson_sd, ymax = invsimpson_mean + invsimpson_sd, group = 1)) +
  geom_linerange() +
  geom_line(size = 0.1) +
  geom_point(size = 0.3, color = "red") +
  theme_light() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top") +
  xlab("") +
  ylab("Inv Simpson")


osd2014_asv_simpson_e <- osd2014_asv_alpha_summary %>%
  select(label, contains("simpson_e")) %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label)) %>%
  ggplot(aes(label, simpson_e_mean,  ymin = simpson_e_mean - simpson_e_sd, ymax = simpson_e_mean + simpson_e_sd, group = 1)) +
  geom_linerange() +
  geom_line(size = 0.1) +
  geom_point(size = 0.3, color = "red") +
  theme_light() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top") +
  xlab("") +
  ylab("Inv Simpson")

osd2014_asv_pielou <- osd2014_asv_alpha_summary %>%
  select(label, contains("pielou")) %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label)) %>%
  ggplot(aes(label, pielou_mean,  ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd, group = 1)) +
  geom_linerange() +
  geom_line(size = 0.1) +
  geom_point(size = 0.3, color = "red") +
  theme_light() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top") +
  xlab("") +
  ylab("Pielou's J")



# Rarefaction curves ------------------------------------------------------
calculate_rarefaction_curves <- function(physeq, measures, depths, parallel=FALSE, ncpus=1) {
  estimate_rarified_richness <- function(physeq, measures, depth) {
    if(max(sample_sums(physeq)) < depth) return()
    physeq <- prune_samples(sample_sums(physeq) >= depth, physeq)

    rarified_physeq <- rarefy_even_depth(physeq, depth, verbose = FALSE)

    alpha_diversity <- estimate_richness(rarified_physeq, measures = measures)

    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- alpha_diversity %>%
      as_tibble() %>%
      mutate(Sample = sample_names(physeq)) %>%
      gather(Measure, Alpha_diversity, -Sample)

    molten_alpha_diversity
  }

  if (parallel){
    #if parallel setup the cluster
    library(doParallel)
    print("Running Calculation in Parallel...")
    cl <- makeCluster(ncpus)
    registerDoParallel(cl)
  }

  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- plyr::ldply(depths, estimate_rarified_richness,
                                        physeq = physeq, measures = measures,
                                        .id = 'Depth',
                                        .paropts = list(.packages = c('phyloseq', 'vegan','reshape2', 'tidyverse')),
                                        .progress = ifelse(interactive(), 'text', 'none'),
                                        .parallel = parallel)

  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]

  #add a standard deviation column
  rarefaction_curve_data_summary <- plyr::ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'),
                                                plyr::summarise,
                                                Alpha_diversity_mean = mean(Alpha_diversity),
                                                Alpha_diversity_sd = sd(Alpha_diversity))

  #add the sample data
  rarefaction_curve_data_summary_verbose <- rarefaction_curve_data_summary
  return(rarefaction_curve_data_summary_verbose)
}

osd2014_dada2_phyloseq_alpha_rarecurv <- calculate_rarefaction_curves(physeq = osd2014_dada2_phyloseq_alpha, measures = c('Observed', 'Shannon', 'InvSimpson'), depths =  rep(c(1, 10, 100, 1:100 * 1000), each = 100), parallel = TRUE, ncpus = 4)

osd2014_asv_rarecurve <- ggplot(data = osd2014_dada2_phyloseq_alpha_rarecurv %>% inner_join(osd2014_cdata %>% dplyr::rename(Sample = label)), mapping = aes(x = Depth,
                                                                                                                                                            y = Alpha_diversity_mean,
                                                                                                                                                            ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                                                                                                                                                            ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                                                                                                                                                            group = Sample)) +
  geom_line() +
  geom_ribbon(alpha = 0.6) +
  facet_grid(facets = Measure ~ meow_province, scales = 'free_y') + theme_bw()

ggsave(osd2014_asv_rarecurve, filename = "osd2014_16S_asv/figures/osd2014_asv_rarecurve.pdf", width = 11.69, height = 8.27)

# ASV profiles ------------------------------------------------------------
osd2014_dada2_phyloseq_alpha_prop <- transform_sample_counts(osd2014_dada2_phyloseq_alpha, function(x) x/sum(x))

osd2014_dada2_phyloseq_phylum_prop_df <- psmelt(osd2014_dada2_phyloseq_alpha_prop) %>%
  as_tibble() %>%
  select(Sample, Abundance, Phylum, asv_name) %>%
  rename(label = Sample, prop = Abundance, phylum = Phylum) %>%
  filter(prop > 0) %>%
  group_by(label, phylum) %>%
  mutate(phyl_agg = sum(prop)) %>%
  ungroup() %>%
  group_by(phylum) %>%
  mutate(phyl_agg_glob = sum(prop)/nsamples(osd2014_dada2_phyloseq_alpha)) %>%
  ungroup()

osd2014_dada2_phyloseq_phylum_prop_df <- osd2014_dada2_phyloseq_phylum_prop_df %>%
  mutate(phylum = ifelse(phyl_agg <= 0.01, "Other", as.character(phylum)),
         label = fct_relevel(label, st_100_order_terrestrial$label),
         phylum = fct_reorder(phylum, -phyl_agg_glob),
         phylum = fct_relevel(phylum, 'Other', after = Inf))

colors <- c(
  "#cf4149", #Proteobacteria
  "#66b14a", #Bacteroidetes
  "#8f62ca", #Firmicutes
  "#6b8bcd", #Cyanobacteria
  "#c2ad4b", #Chlamydiae
  "#ce7531", #Actinobacteria
  "#4baf90", #Euryarchaeota
  "#c85d9d", #Verrucomicrobia
  "#542437", #No rel
  "#c26f65", #Planctomycetes
  "#A5C990", #Chlorobi
  "#767834", #Marinimicr
  "#559279", #Tenericutes
  "#D3BBC3", #Lentisp
  "#7f8c8d") #Other

osd2014_dada2_phyloseq_phylum_prop_plot <- ggplot(osd2014_dada2_phyloseq_phylum_prop_df %>% group_by(label, phylum) %>% summarise(N = sum(prop)), aes(label, N, fill = phylum, group = phylum)) +
  #geom_area(position = "stack") +
  geom_col(color = "grey20", size = 0.2, width = 1) +
  scale_x_discrete(expand = c(0, 0), limits = st_100_order_terrestrial$label) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_fill_manual(values = (colors)) +
  ylab("Proportion") +
  xlab("") +
  labs(fill="") +
  #theme_bw() +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top")
osd2014_asv_barplot <- ggpubr::ggarrange(osd2014_asv_invsimpson, osd2014_asv_obsrich, osd2014_dada2_phyloseq_phylum_prop_plot, ncol = 1, nrow = 3, heights = c(0.15,0.15, 0.6))
ggsave(osd2014_asv_barplot, filename = "osd2014_16S_asv/figures/osd2014_asv_barplot.pdf", width = 11.69, height = 8.27)


# Collector curves --------------------------------------------------------

specacc <- function(comm, permutations){
  x <- comm
  x <- as.matrix(x)
  x <- x[, colSums(x) > 0, drop = FALSE]
  n <- nrow(x)
  p <- ncol(x)

  accumulator <- function(X, x, ind, permat) {
    rowSums(apply(x[permat[X,], ], 2, cumsum) > 0)
  }
  specaccum <- sdaccum <- sites <- perm <- NULL
  permat <- vegan:::getPermuteMatrix(permutations, n)
  perm <- t(do.call( "rbind", (pbmclapply(X = 1:nrow(permat), accumulator, x = x, perma = permat, mc.cores = 4))))
  specaccum <- apply(perm, 1, mean)
  sdaccum <- apply(perm, 1, sd)

  data_frame(richness_mean = specaccum, richness_sd = sdaccum) %>%
    mutate(nsample = row_number())
}

osd2014_dada2_phyloseq_alpha_rarecurve <- specacc(comm = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_alpha), permutations = 999)

rarecurves <- vector(mode = "list")
osd2014_dada2_phyloseq_alpha_r <- rarefy_even_depth(prune_samples(osd2014_cdata %>% dplyr::filter(meow_province %in% osd2014_meow_regions$meow_province) %>% .$label, osd2014_dada2_phyloseq_alpha))
for (i in osd2014_meow_regions$meow_province){
  osd2014_dada2_phyloseq_alpha_f <- prune_samples(osd2014_cdata %>% dplyr::filter(meow_province == i) %>% .$label, osd2014_dada2_phyloseq_alpha_r)
  rarecurves[[i]] <- specacc(comm = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_alpha_f), permutations = 9999) %>% mutate(meow_province = i)
}


osd2014_dada2_phyloseq_alpha_rarecurve_plot_all <- ggplot(bind_rows(rarecurves), aes(nsample, richness_mean)) +
  geom_ribbon(aes(ymin = richness_mean - richness_sd, ymax = richness_mean + richness_sd), alpha = 0.4) +
  geom_line() +
  facet_wrap(~meow_province) +
  theme_light() +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  ylab("# of ASVs") +
  xlab("# of samples")

osd2014_dada2_phyloseq_alpha_rarecurve_plot <- ggplot(osd2014_dada2_phyloseq_alpha_rarecurve, aes(nsample, richness_mean)) +
  geom_ribbon(aes(ymin = richness_mean - richness_sd, ymax = richness_mean + richness_sd), alpha = 0.4) +
  geom_line() +
  theme_light() +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  ylab("# of ASVs") +
  xlab("# of samples")

ggsave(osd2014_dada2_phyloseq_alpha_rarecurve_plot, filename = "osd2014_16S_asv/figures/osd2014_dada2_phyloseq_alpha_rarecurve_plot.pdf", width = 11.69, height = 8.27)

# DivNet ------------------------------------------------------------------

library(tidyverse)
library(DivNet)
library(phyloseq)

osd2014_dada2_phyloseq_alpha_diversity <- prune_samples(osd2014_cdata$label, osd2014_dada2_phyloseq_alpha)
osd2014_dada2_phyloseq_alpha_diversity <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_alpha_diversity) > 0, osd2014_dada2_phyloseq_alpha_diversity)
taxa_sums <- taxa_sums(osd2014_dada2_phyloseq_alpha_diversity)
taxa_unobserved <- microbiome::prevalence(osd2014_dada2_phyloseq_alpha_diversity)
taxa_unobserved <- ifelse(taxa_unobserved < 0.9, 0, 1)
base_asv <- (which.max(taxa_sums*taxa_unobserved))

# osd2014_divnet_asv_provinces <-  osd2014_dada2_phyloseq_alpha_diversity %>%
#   divnet(X = "meow_province", ncores = 4, base = base_asv)


osd2014_dada2_phyloseq_alpha_genus <- tax_glom(osd2014_dada2_phyloseq_alpha_diversity, taxrank = "Genus")
osd2014_dada2_phyloseq_alpha_genus <- prune_samples((osd2014_cdata$label), osd2014_dada2_phyloseq_alpha_genus)
osd2014_dada2_phyloseq_alpha_genus <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_alpha_genus) > 0, osd2014_dada2_phyloseq_alpha_genus)

osd2014_asv_alpha_genus_provinces <-  osd2014_dada2_phyloseq_alpha_genus %>%
  divnet(X = "meow_province", ncores = 4)

osd2014_asv_alpha_genus_provinces1 <-  osd2014_dada2_phyloseq_alpha_genus %>%
  divnet(ncores = 4)

osd2014_asv_alpha_genus_provinces_plugin <- pbmclapply(1:trials, get_div, physeq = osd2014_dada2_phyloseq_alpha_genus, mc.cores = 4)

osd2014_frequencytablelist_genus_bw <- build_frequency_count_tables(t(otu_table(osd2014_dada2_phyloseq_alpha_genus)))
alpha_ests <- pbmclapply(X = osd2014_frequencytablelist_genus_bw, get_alpha)
osd2014_asv_alpha_summary_genus_bw <- bind_rows(alpha_ests, .id = "label")


osd2014_asv_alpha_genus_provinces_plugin_summary <- bind_rows(osd2014_asv_alpha_genus_provinces_plugin) %>%
  select(label, richness, invsimpson, shannon, pielou, simpson_e) %>%
  group_by(label) %>%
  summarise_all(funs("mean", "sd"))

library(breakaway)
osd2014_frequencytablelist_asv <- build_frequency_count_tables(t(otu_table(osd2014_dada2_phyloseq_alpha)))
osd2014_breakaway_estimates_asv <- lapply(osd2014_frequencytablelist_asv, breakaway)
names(osd2014_breakaway_estimates_asv) <- names(osd2014_frequencytablelist_asv)
osd2014_objective_bayes_negbin_estimates_asv <- lapply(osd2014_frequencytablelist_asv, objective_bayes_negbin, answers = TRUE, output = FALSE)

osd2014_frequencytablelist_genus <- build_frequency_count_tables(t(otu_table(prune_samples((osd2014_cdata$label), tax_glom(osd2014_dada2_phyloseq_alpha_diversity, taxrank = "Genus")))))
osd2014_breakaway_estimates_genus <- lapply(osd2014_frequencytablelist_genus, breakaway)
names(osd2014_breakaway_estimates_genus) <- names(osd2014_frequencytablelist_genus)
osd2014_objective_bayes_negbin_estimates_genus <- lapply(osd2014_frequencytablelist_genus, objective_bayes_negbin, answers = TRUE, output = FALSE)

osd2014_objective_bayes_negbin_estimates_asv_summary <- osd2014_breakaway_estimates_asv %>%
  map_df("est") %>% gather(label, est) %>%
  inner_join(osd2014_objective_bayes_negbin_estimates_asv %>%
               map_df("mean") %>% gather(label, mean)) %>%
  inner_join(osd2014_objective_bayes_negbin_estimates_asv %>%
               map_df("semeanest") %>% gather(label, semeanest)) %>%
  inner_join(osd2014_objective_bayes_negbin_estimates_asv %>%
               map_df("ci", gather) %>% t %>% as_tibble(rownames = "label") %>%
               rename(l95 = V1, u95 = V2))


osd2014_objective_bayes_negbin_estimates_genus_summary <- osd2014_objective_bayes_negbin_estimates_genus %>%
  map_df("est") %>% gather(label, est) %>%
  inner_join(osd2014_objective_bayes_negbin_estimates_genus %>%
               map_df("mean") %>% gather(label, mean)) %>%
  inner_join(osd2014_objective_bayes_negbin_estimates_genus %>%
               map_df("semeanest") %>% gather(label, semeanest)) %>%
  inner_join(osd2014_objective_bayes_negbin_estimates_genus %>%
               map_df("ci", gather) %>% t %>% as_tibble(rownames = "label") %>%
               rename(l95 = V1, u95 = V2))


osd2014_shannon_estimates_provinces <- data.frame("label" = names(osd2014_asv_alpha_genus_provinces$shannon),
                                                  "DivNet.est" = osd2014_asv_alpha_genus_provinces$shannon,
                                                  "DivNet.var" = osd2014_asv_alpha_genus_provinces$`shannon-variance`) %>%
  left_join(osd2014_sample_cohesion) %>%
  inner_join(osd2014_asv_alpha_summary_genus_bw %>% select(label, est, seest) %>% rename(shannon = est, shannon_se = seest)) %>%
  inner_join(osd2014_objective_bayes_negbin_estimates_genus_summary) %>%
  inner_join(osd2014_cdata) %>%
  as_tibble()


osd2014_shannon_estimates_provinces_asv <-  data.frame("label" = names(osd2014_divnet_asv_provinces$shannon),
                                                       "DivNet.est" = osd2014_divnet_asv_provinces$shannon,
                                                       "DivNet.var" = osd2014_divnet_asv_provinces$`shannon-variance`) %>%
  left_join(osd2014_sample_cohesion) %>%
  inner_join(osd2014_asv_alpha_summary %>% select(label, shannon_mean, shannon_sd) %>% rename(shannon = shannon_mean, shannon_se = shannon_sd)) %>%
  inner_join(osd2014_objective_bayes_negbin_estimates_asv_summary) %>%
  inner_join(osd2014_cdata) %>%
  inner_join(osd2014_asv_alpha_summary %>% select(label, shannon_mean, shannon_sd) %>% rename(shannon = shannon_mean, shannon_se = shannon_sd)) %>%
  as_tibble()



# osd2014_invsimpson_estimates_provinces <- data.frame("meow_province" = osd2014_dada2_phyloseq_alpha_genus %>% sample_data %>%
#                                                     get_variable("meow_province"),
#                                                   "sample" = osd2014_dada2_phyloseq_alpha_genus %>%
#                                                     estimate_richness(measures = "Simpson") %$% Simpson,
#                                                   "DivNet.est" = osd2014_divnet_genus_provinces$simpson,
#                                                   "DivNet.var" = osd2014_divnet_genus_provinces$`simpson-variance`)

colors1 <- rev(c("#C059D1", "#D2D0B7", "#B99ED0", "#ABE172", "#7AD4CB",  "#DB8763"))

osd2014_shannon_estimates_provinces <- osd2014_shannon_estimates_provinces %>%
  mutate(meow_province = fct_relevel(meow_province, rev(c("Tropical Northwestern Atlantic", "Warm Temperate Northwest Atlantic", "Cold Temperate Northwest Atlantic",
                                                          "Lusitanian", "Mediterranean Sea", "Northern European Seas"))))
osd2014_shannon_estimates_provinces_asv <- osd2014_shannon_estimates_provinces_asv %>%
  mutate(meow_province = fct_relevel(meow_province, rev(c("Tropical Northwestern Atlantic", "Warm Temperate Northwest Atlantic", "Cold Temperate Northwest Atlantic",
                                                          "Lusitanian", "Mediterranean Sea", "Northern European Seas"))))

beta_plot <- ggplot(osd2014_shannon_estimates_provinces_asv, aes(y = shannon, x = (meow_province), fill = meow_province)) +
  #geom_point(aes(y = sample, color = meow_provinces, fill = meow_provinces), shape = 21, alpha = 0.5) +
  # geom_density_ridges(data = osd2014_shannon_estimates_provinces, aes(shannon, meow_province, color = meow_province, fill = meow_province),
  #                     scale = 0.7, alpha = .7, rel_min_height = 0.0 ,panel_scaling = T, size = 0.2, color = "black",
  #                     jittered_points = TRUE, point_size = 0.5
  # ) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 1,
                        color = "black", alpha = 1, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4) +
  geom_errorbar(aes(x = (meow_province), ymin = DivNet.est - 2*sqrt(DivNet.var),
                    ymax = DivNet.est + 2*sqrt(DivNet.var)), col = "red", width = .1) +
  geom_point(aes(y = DivNet.est, x = meow_province), col = "red", shape = 21, fill = "grey", size = 1, stroke = 0.5) +
  ylab("") +
  xlab("Shannon diversity estimate\n(genus level)") +
  scale_fill_manual(values = colors1) +
  scale_color_manual(values = colors1) +
  theme_bw() +
  theme(legend.position = "none")

osd2014_shannon_estimates_provinces_asv <- osd2014_shannon_estimates_provinces_asv %>%
  mutate(meow_province = fct_relevel(meow_province, (c("Tropical Northwestern Atlantic", "Warm Temperate Northwest Atlantic", "Cold Temperate Northwest Atlantic",
                                                       "Lusitanian", "Mediterranean Sea", "Northern European Seas"))))

colors1 <- (c("#C059D1", "#D2D0B7", "#B99ED0", "#ABE172", "#7AD4CB",  "#DB8763"))


lm_y <- lm(test$start_lat ~ test$mean)
gam_y <- (mgcv::gam(test$start_lat ~ s(test$mean, bs = "cs")))
AIC(lm_y, gam_y)
anova(lm_y, gam_y, test="F")

lm_y <- lm(test$water_temperature ~ test$mean)
gam_y <- (mgcv::gam(test$water_temperature ~ s(test$mean, bs = "cs")))
AIC(lm_y, gam_y)
anova(lm_y, gam_y, test="F")


rich_lat <- osd2014_shannon_estimates_provinces_asv %>%
  filter(label %in% osd2014_cdata_filt$label) %>%
  ggplot(aes((start_lat), mean)) +
  geom_smooth( method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = mean - semeanest, ymax = mean + semeanest), width = 0.5) +
  geom_point(aes(fill = meow_province), shape = 21, color = "black", size = 1.5, alpha = 1, stroke = 0.5) +
  theme_bw() +
  scale_fill_manual(values = colors1) +
  ylab("Richness") +
  xlab("Latitude") +
  theme(legend.position = "none")

rich_temp <- osd2014_shannon_estimates_provinces_asv %>%
  filter(label %in% osd2014_cdata_filt$label) %>%
  ggplot(aes((water_temperature), mean)) +
  geom_smooth( method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = mean - semeanest, ymax = mean + semeanest), width = 0.5) +
  geom_point(aes(fill = meow_province), shape = 21, color = "black", size = 1.5, alpha = 1, stroke = 0.5) +
  theme_bw() +
  scale_fill_manual(values = colors1) +
  ylab("Richness") +
  xlab("Water temperature (C)") +
  theme(legend.position = "none")


shan_lat <- osd2014_shannon_estimates_provinces_asv %>%
  filter(label %in% osd2014_cdata_filt$label) %>%
  ggplot(aes((start_lat), shannon)) +
  geom_smooth( method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = shannon - shannon_se, ymax = shannon + shannon_se), width = 0) +
  geom_point(aes(fill = meow_province), shape = 21, color = "black", size = 1.5, alpha = 1, stroke = 0.5) +
  theme_bw() +
  scale_fill_manual(values = colors1) +
  ylab("Shannon") +
  xlab("Latitude") +
  theme(legend.position = "none")

shan_temp <- osd2014_shannon_estimates_provinces_asv %>%
  filter(label %in% osd2014_cdata_filt$label) %>%
  ggplot(aes((water_temperature), shannon)) +
  geom_smooth( method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = shannon - shannon_se, ymax = shannon + shannon_se), width = 0) +
  geom_point(aes(fill = meow_province), shape = 21, color = "black", size = 1.5, alpha = 1, stroke = 0.5) +
  theme_bw() +
  scale_fill_manual(values = colors1) +
  ylab("Shannon") +
  xlab("Water temperature (C)") +
  theme(legend.position = "none")

testDiversity(osd2014_asv_alpha_genus_provinces, "shannon")
testDiversity(osd2014_asv_alpha_genus_provinces, "invsimpson")


interpolate_alpha(osd2014_shannon_estimates_provinces_asv %>% select(label, est) %>% rename(mean = est) %>% filter(label %in% osd2014_cdata_filt$label))
interpolate_alpha(osd2014_shannon_estimates_provinces_asv %>% select(label, shannon) %>% rename(mean = shannon) %>% filter(label %in% osd2014_cdata_filt$label))


p1 <- osd2014_dada2_phyloseq_phylum_prop_plot
ggsave(plot = p1, filename = "osd2014_16S_asv/figures/osd2014_dada2_phyloseq_alpha_barplot.pdf", width = 11.69, height = 8.27)

p2 <- ggarrange(rich_lat , rich_temp , shan_lat , shan_temp)
ggsave(p2, filename = "osd2014_16S_asv/figures/osd2014_dada2_phyloseq_alpha_scatter.pdf", width = 5, height = 4)

ggsave(beta_plot, filename = "osd2014_16S_asv/figures/osd2014_dada2_phyloseq_alpha_betadiv_genus.pdf", width = 4, height = 3)

# BEGIN: Save objects ------------------------------------------------------------
# WARNING!!! You might not want to run this code --------------------------
save.image("osd2014_16S_asv/data/osd2014_16S_alpha_diversity.Rdata")
# END: Save objects ------------------------------------------------------------
