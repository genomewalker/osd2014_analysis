library(phyloseq)
library(vegan)
library(tidyverse)
library(breakaway)
library(RPostgreSQL)

load(file = "osd2014_shotgun/data/osd2014_eggnog_physeq_filt_objects.Rdata")

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)

st_100_order_terrestrial <- tbl(my_db, "osd2014_st_order_terrestrial") %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label)

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

source("https://raw.githubusercontent.com/colinbrislawn/phyloseq/f33259c9291434b965afa99b477dad0c751dad95/R/extend_vegan.R")
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

library(phyloseq)
library(pbmcapply)
trials <- 1000
physeq_filt_prev <- prune_taxa(taxa_sums(osd2014_eggnog_phyloseq_alpha) > 0, osd2014_eggnog_phyloseq_alpha)
min_lib <-   min(sample_sums(physeq_filt_prev))
osd2014_eggnog_alpha <- pbmclapply(1:trials, get_div, physeq = physeq_filt_prev, mc.cores = 4)

osd2014_eggnog_alpha_summary <- bind_rows(osd2014_asv_alpha) %>%
  select(label, richness, invsimpson, shannon, pielou, simpson_e) %>%
  group_by(label) %>%
  summarise_all(funs("mean", "sd"))

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


osd2014_asv_invsimpson <- d16S_alpha_summary %>%
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


osd2014_asv_simpson_e <- d16S_alpha_summary %>%
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



osd2014_asv_pielou<- d16S_alpha_summary %>%
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
  rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary,
                                                  data.frame(sample_data(physeq)),
                                                  by.x = 'Sample',
                                                  by.y = 'row.names')
  return(rarefaction_curve_data_summary_verbose)
}

osd2014_dada2_phyloseq_alpha_rarecurv <- calculate_rarefaction_curves(physeq = osd2014_dada2_phyloseq_alpha, measures = c('Observed', 'Shannon', 'InvSimpson'), depths =  rep(c(1, 10, 100, 1:100 * 1000), each = 100), parallel = TRUE, ncpus = 4)

osd2014_asv_rarecurve <- ggplot(data = osd2014_dada2_phyloseq_alpha_rarecurv, mapping = aes(x = Depth,
                                                                                            y = Alpha_diversity_mean,
                                                                                            ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                                                                                            ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                                                                                            group = Sample)) +
  geom_line() +
  geom_ribbon(alpha = 0.6) +
  facet_grid(facets = Measure ~ osd_meow_ecoregion, scales = 'free_y') + theme_bw()

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

osd2014_dada2_phyloseq_phylum_prop_plot <- ggplot(osd2014_dada2_phyloseq_phylum_prop_df, aes(label, prop, fill = phylum, group = phylum)) +
  #geom_area(position = "stack") +
  geom_bar(stat = "identity", width = 1, size = 0.000001) +
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

osd2014_asv_barplot <- ggarrange(osd2014_asv_invsimpson, osd2014_asv_obsrich, osd2014_dada2_phyloseq_phylum_prop_plot, ncol = 1, nrow = 3, heights = c(0.15,0.15, 0.6))
ggsave(osd2014_asv_barplot, filename = "osd2014_16S_asv/figures/osd2014_asv_barplot.pdf", width = 11.69, height = 8.27)
save.image("osd2014_16S_asv/data/osd2014_16S_alpha_diversity.Rdata")
