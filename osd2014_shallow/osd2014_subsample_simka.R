library(tidyverse)

simka_1M <- read.table(file = "~/Downloads/k21/allVSall_1M/mat_abundance_braycurtis.csv.gz", header = T, row.names = 1, sep = ";", check.names = F) %>%
  as.dist() %>%
  broom::tidy() %>%
  as_tibble() %>%
  separate(item1, into=c("sample1", "depth", "rep"), sep = '\\.', remove = TRUE, extra = "drop") %>%
  separate(item2, into=c("sample2"), sep = '\\.', remove = TRUE, extra = "drop")

simka_5M <- read.table(file = "~/Downloads/k21/allVSall_5M/mat_abundance_braycurtis.csv.gz", header = T, row.names = 1, sep = ";", check.names = F) %>%
  as.dist() %>%
  broom::tidy() %>%
  as_tibble() %>%
  separate(item1, into=c("sample1", "depth", "rep"), sep = '\\.', remove = TRUE, extra = "drop") %>%
  separate(item2, into=c("sample2"), sep = '\\.', remove = TRUE, extra = "drop")

simka_10M <- read.table(file = "~/Downloads/k21/allVSall_10M/mat_abundance_braycurtis.csv.gz", header = T, row.names = 1, sep = ";", check.names = F) %>%
  as.dist() %>%
  broom::tidy() %>%
  as_tibble() %>%
  separate(item1, into=c("sample1", "depth", "rep"), sep = '\\.', remove = TRUE, extra = "drop") %>%
  separate(item2, into=c("sample2"), sep = '\\.', remove = TRUE, extra = "drop")

simka_15M <- read.table(file = "~/Downloads/k21/allVSall_15M/mat_abundance_braycurtis.csv.gz", header = T, row.names = 1, sep = ";", check.names = F) %>%
  as.dist() %>%
  broom::tidy() %>%
  as_tibble() %>%
  separate(item1, into=c("sample1", "depth", "rep"), sep = '\\.', remove = TRUE, extra = "drop") %>%
  separate(item2, into=c("sample2"), sep = '\\.', remove = TRUE, extra = "drop")

simka_20M <- read.table(file = "~/Downloads/k21/allVSall_20M/mat_abundance_braycurtis.csv.gz", header = T, row.names = 1, sep = ";", check.names = F) %>%
  as.dist() %>%
  broom::tidy() %>%
  as_tibble() %>%
  separate(item1, into=c("sample1", "depth", "rep"), sep = '\\.', remove = TRUE, extra = "drop") %>%
  separate(item2, into=c("sample2"), sep = '\\.', remove = TRUE, extra = "drop")

bind_rows(simka_1M, simka_5M, simka_10M, simka_15M, simka_20M) %>%
  mutate(depth = fct_relevel(depth, c("1M", "5M", "10M", "15M", "20M")),
         sample1 = gsub("-NPL022_R1_clipped", "", sample1),
         sample2 = gsub("-NPL022_R1_clipped", "", sample2)) %>%
  filter(sample1 == sample2) %>%
  ungroup() %>%
  group_by(sample1, depth) %>%
  dplyr::summarise(mean = mean(distance), sd = sd(distance)) %>%
  ggplot(aes(depth, mean)) +
  geom_errorbar(aes(ymin=mean - sd, ymax=mean + sd),  width=0.1) +
  geom_line(aes(group = sample1)) +
  geom_point(shape = 21, fill = "grey") +
  facet_wrap(~sample1) +
  theme_linedraw() +
  ylab("Simka average BC dissimilarity")


bind_rows(simka_1M, simka_5M, simka_10M, simka_15M, simka_20M) %>%
  mutate(depth = fct_relevel(depth, c("1M", "5M", "10M", "15M", "20M")),
         sample1 = gsub("-NPL022_R1_clipped", "", sample1),
         sample2 = gsub("-NPL022_R1_clipped", "", sample2)) %>%
  filter(sample1 != sample2) %>%
  rowwise() %>%
  mutate(comp = paste(sample1,sample2)) %>%
  #ungroup() %>%
  group_by(comp, depth) %>%
  summarise(mean = mean(distance), sd = sd(distance)) %>%
  ggplot(aes(depth, mean, color = comp)) +
  geom_errorbar(aes(ymin=mean - sd, ymax=mean + sd),  width=0.1) +
  geom_line(aes(group = comp)) +
  geom_point(shape = 21, fill = "grey") +
  scale_y_continuous(limits = c(0.5,1)) +
  facet_wrap(~comp) +
  theme_linedraw() +
  theme(legend.position = "bottom")


mitags <- read_tsv("~/Downloads/all_feature-table.tsv", col_names = TRUE) %>%
  as.data.frame() %>%
  column_to_rownames("OTU_ID")

mitags_1M <- mitags[, grep('.1M.', names(mitags))]
mitags_5M <- mitags[, grep('.5M.', names(mitags))]
mitags_10M <- mitags[, grep('.10M.', names(mitags))]
mitags_15M <- mitags[, grep('.15M.', names(mitags))]
mitags_20M <- mitags[, grep('.20M.', names(mitags))]


vst_trans <- function(X){
  mitags_s <- data.frame(label = rownames(X))
  row.names(mitags_s) <- mitags_s$label
  mitags_phy <- phyloseq(otu_table(X, taxa_are_rows = F), sample_data(mitags_s))
  mitags_phy <- prune_taxa(taxa_sums(mitags_phy) > 0, mitags_phy)
  # DESEQ2 VST and normalised counts
  mitags_phy_deseq <- phyloseq_to_deseq2(mitags_phy, ~1)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans <- apply(counts(mitags_phy_deseq), 1, gm_mean)
  diagdds <- estimateSizeFactors(mitags_phy_deseq, type = "poscounts")
  diagdds <- estimateDispersions(diagdds, fitType='local')
  diagvst <- varianceStabilizingTransformation(diagdds, blind = FALSE)
  diagvst <- assay(diagvst)
  diagvst[diagvst < 0] <- 0
  mitags_phy_deseq_vst <- mitags_phy

  otu_table(mitags_phy_deseq_vst) <- otu_table(diagvst, taxa_are_rows = TRUE)
  return(mitags_phy_deseq_vst)
}

mitags_1M_dist <- phyloseq::distance(vst_trans(t(mitags_1M)), method = "bray") %>% broom::tidy() %>% as_tibble() %>%
  separate(item1, into=c("sample1", "depth", "rep"), sep = '\\.', remove = TRUE, extra = "drop") %>%
  separate(item2, into=c("sample2"), sep = '\\.', remove = TRUE, extra = "drop")

mitags_5M_dist <- phyloseq::distance(vst_trans(t(mitags_5M)), method = "bray") %>% broom::tidy() %>% as_tibble() %>%
  separate(item1, into=c("sample1", "depth", "rep"), sep = '\\.', remove = TRUE, extra = "drop") %>%
  separate(item2, into=c("sample2"), sep = '\\.', remove = TRUE, extra = "drop")

mitags_10M_dist <- phyloseq::distance(vst_trans(t(mitags_10M)), method = "bray") %>% broom::tidy() %>% as_tibble() %>%
  separate(item1, into=c("sample1", "depth", "rep"), sep = '\\.', remove = TRUE, extra = "drop") %>%
  separate(item2, into=c("sample2"), sep = '\\.', remove = TRUE, extra = "drop")

mitags_15M_dist <- phyloseq::distance(vst_trans(t(mitags_15M)), method = "bray") %>% broom::tidy() %>% as_tibble() %>%
  separate(item1, into=c("sample1", "depth", "rep"), sep = '\\.', remove = TRUE, extra = "drop") %>%
  separate(item2, into=c("sample2"), sep = '\\.', remove = TRUE, extra = "drop")

mitags_20M_dist <- phyloseq::distance(vst_trans(t(mitags_20M)), method = "bray") %>% broom::tidy() %>% as_tibble() %>%
  separate(item1, into=c("sample1", "depth", "rep"), sep = '\\.', remove = TRUE, extra = "drop") %>%
  separate(item2, into=c("sample2"), sep = '\\.', remove = TRUE, extra = "drop")


mitags_dist <- bind_rows(mitags_1M_dist, mitags_5M_dist,
          mitags_10M_dist, mitags_15M_dist,
          mitags_20M_dist) %>%
  mutate(depth = fct_relevel(depth, c("1M", "5M", "10M", "15M", "20M")),
         sample1 = gsub("-NPL022_R1_clipped", "", sample1),
         sample2 = gsub("-NPL022_R1_clipped", "", sample2))
mitags_dist %>%
  filter(sample1 != sample2) %>%
  rowwise() %>%
  mutate(comp = paste(sample1,sample2)) %>%
  #ungroup() %>%
  group_by(depth) %>%
  summarise(mean = mean(distance), sd = sd(distance)) %>%
  ggplot(aes(depth, mean, color = comp)) +
  geom_errorbar(aes(ymin=mean - sd, ymax=mean + sd),  width=0.1) +
  #geom_boxplot() +
  geom_line(aes(group = comp)) +
  geom_point(shape = 21, fill = "grey") +
  #scale_y_continuous(limits = c(0.5,1)) +
  facet_wrap(~comp) +
  theme_linedraw() +
  theme(legend.position = "bottom")




mitags_dist %>%
  filter(sample1 == sample2) %>%
  ungroup() %>%
  group_by(sample1, depth) %>%
  summarise(mean = mean(distance), sd = sd(distance)) %>%
  ggplot(aes(depth, mean)) +
  geom_errorbar(aes(ymin=mean - sd, ymax=mean + sd),  width=0.1) +
  geom_line(aes(group = sample1)) +
  geom_point(shape = 21, fill = "grey") +
  facet_wrap(~sample1) +
  theme_linedraw()


mitags_dist_simka <- bind_rows(simka_1M, simka_5M, simka_10M, simka_15M, simka_20M) %>%
  mutate(depth = fct_relevel(depth, c("1M", "5M", "10M", "15M", "20M")),
         sample1 = gsub("-NPL022_R1_clipped", "", sample1),
         sample2 = gsub("-NPL022_R1_clipped", "", sample2))

mitags_dist_simka %>% inner_join(mitags_dist %>% dplyr::rename(distance_mt = distance)) %>%
  filter(sample1 == sample2) %>%
  ggplot(aes(distance_mt, distance)) +
  geom_point(shape = 21, fill ="red", alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0,1)) + ylim(c(0,1)) +
  facet_wrap(~depth, nrow = 1)
