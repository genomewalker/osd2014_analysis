library(vegan)
library(tidyverse)
library(phyloseq)
library(spaa)
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata", verbose = TRUE)

osd2014_16s_otuXsample_physeq_filt_prev_beta_df <- as(otu_table(osd2014_dada2_phyloseq_beta), "matrix") %>% as.data.frame()


osd2014_dada2_phyloseq_beta_stats <- psmelt(osd2014_dada2_phyloseq_beta) %>%
  tbl_df() %>%
  select(Sample, OTU, Abundance) %>%
  group_by(Sample) %>%
  mutate(proportion_n = (Abundance/sum(Abundance))) %>%
  ungroup() %>%
  group_by(OTU) %>%
  mutate(proportion = (Abundance/sum(Abundance)), n = length(unique(Sample))) %>%
  mutate(mean_proportion = mean(proportion_n)) %>%
  dplyr::select(OTU, proportion, proportion_n, mean_proportion, n) %>%
  mutate(B = 1/(sum(proportion^2)), B_a = (B-1)/(n-1), B_a = ifelse(is.infinite(B_a), 0, B_a)) %>%
  #select(-prop ortion) %>%
  unique()

comm.tab <- phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta)
comm.tab <- comm.tab[,which(colSums(comm.tab)>0)]
comm.tab <- comm.tab[which(rowSums(comm.tab)>0),]


get_random <- function(x, y){
  null <- nullmodel(y, method = "quasiswap_count")
  res <- simulate(null, nsim=1)
  l <- as_data_frame(niche.width(res,  method = "levins"))
  return(l)
}

# levin.index.simul<-plyr::ldply(1:100, get_random, y = comm.tab, .parallel = T)
levin.index.simul <- pbmcapply::pbmclapply(1:1000, get_random, y = comm.tab, mc.cores = 4)
levin.index.simul <- lapply(1, get_random, y = comm.tab)

levin.index.simul <- data.table::rbindlist(levin.index.simul)

levin.index.real <- as.numeric(niche.width(comm.tab, method = "levins"))

colnames(levin.index.simul) <- colnames(comm.tab)
levin.index.simul <- as.data.frame(levin.index.simul)
media <- apply(levin.index.simul, 2, mean)
ci <- apply(levin.index.simul, 2, quantile, probs = c(0.025, 0.975))
results <- data.frame(OTU = colnames(comm.tab), observed = levin.index.real, mean.simulated = media, lowCI = ci[1, ], uppCI = ci[2, ], sign = NA)

results <- results %>% mutate(sign = case_when( observed > uppCI ~ 'Broad',
                                                observed < lowCI ~ 'Narrow',
                                                observed >= lowCI & observed <= uppCI ~ 'Non significant'))%>%
  tbl_df()

results <- results %>%
  left_join(osd2014_dada2_phyloseq_beta_stats) %>%
  mutate(sign = fct_relevel(sign, c("Broad", "Non significant", "Narrow")))

all_plot <- ggplot(results %>% select(mean_proportion, observed, sign) %>% unique, aes(mean_proportion, observed, fill = sign)) +
  geom_point(size=1.7, alpha = 0.9, shape=21, color = "#3A3A3A") + theme_light() +
  scale_x_log10() +
  xlab("Mean proportion (log10)") +
  ylab("Levin's niche breadth (B)") +
  #ggtitle("Environmental unknowns PC distribution") +
  guides(fill = guide_legend(override.aes = list(size=2), nrow = 1)) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        legend.background = element_rect(fill=alpha('white', 0.4)),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.text.align=0) +
  scale_fill_manual(values = c("#588157", "#9A9A9A", "#3891A6"))

save.image("osd2014_16S_asv/data/osd2014_16S_niche_breadth.Rdata")
