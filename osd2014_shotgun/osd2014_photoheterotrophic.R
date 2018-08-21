
library(tidyverse)
library(magrittr)
library(metagenomeSeq)
library(vegan)
library(pvclust)
library(pheatmap)
library(ggdendro)
library(grid)
library(gridExtra)
library(zoo)
library(cowplot)
library(ggtree)
library(geoR)
library(sp)
library(gstat)

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

osd2014_pufm_counts <- tbl(my_db, "osd2014_pufm_counts") %>%
  collect(n = Inf)
osd2014_pr_counts <- tbl(my_db, "osd2014_pr_counts") %>%
  collect(n = Inf)
osd2014_pr_groups <- tbl(my_db, "osd2014_pr_groups") %>%
  collect(n = Inf)
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)
osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)
st_100_order_terrestrial <- tbl(my_db, "osd2014_st_order_coastal") %>%
  collect(n = Inf)


# PufM genes

## Explore data:
missing_samples <- tibble(label = c("OSD152_2014-06-20_5m_NPL022", "OSD15_2014-06-21_0m_NPL022", "OSD80_2014-06-21_2m_NPL022",
                                    "OSD105_2014-06-21_4.12m_NPL022", "OSD146_2014-06-21_5m_NPL022"), group = "None", counts = 1, prop = 1)

osd2014_pufm_counts_prop <- osd2014_pufm_counts %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  droplevels() %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(prop = counts/sum(counts)) %>%
  ungroup() %>%
  arrange(label) %>%
  bind_rows(missing_samples) %>%
  complete(label, group) %>%
  filter(!is.na(group)) %>%
  dplyr::mutate(counts = ifelse(is.na(counts), 0 , counts), prop = ifelse(is.na(prop), 0 , prop), label = fct_relevel(label, st_100_order_terrestrial$label), group = fct_relevel(group, c("Phylogroup C", "Phylogroup D", "Phylogroup E", "Phylogroup F", "Phylogroup G", "Sphingomonadales", "Phylogroup K", "Rubrivivax-like", "Phylogroup I", "Phylogroup J", "Other", "Unaffiliated", "None")))


# Most AAnPB are members of the Alpha-, Beta- and Gammaproteobacteria. Yutin et al. (2007)
# proposed a phylogenetic framework for assigning AAnPB into 12 phylogroups
# (labelled phylogroups A through L), defined according to the structure of the puf operon
# and the phy- logeny of the pufM gene. Different phylogroups have been shown to be abundant
# in different marine environments. For example,

# Alphaproteobacterial phylogroups E and F,  which contain Rhodobacter-like bacteria, and G, which contains Roseobacter-like bacteria,
# are common in nutrient rich coastal waters (Salka et al., 2008; Lehours et al., 2010),

# phylogroups C and D have been mostly observed in offshore waters (Ritchie and Johnson, 2012).

# Phylogroup K, which contains members of the Gammaproteobacteria, also has a widespread distribution, #
# and has been demonstrated to dominate AAnPB composition in samples from the Atlantic and Pacific Oceans
# as well as the Arctic and Mediterranean Seas during summer (Lehours et al., 2010; Ferrera et al., 2014).

pufm_plot <- ggplot(osd2014_pufm_counts_prop, aes(label, prop, fill = group)) +
  geom_col(color = "grey20", size = 0.2, width = 1) +
  scale_x_discrete(expand = c(0, 0), limits = st_100_order_terrestrial$label) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_fill_manual(values = c("Phylogroup C" = "#E9573F", "Phylogroup D" = "#8CC152", "Phylogroup E" =  "#4A89DC", "Phylogroup F" = "#3BAFDA", "Phylogroup G" = "#DA4453", "Sphingomonadales" = "#CFCC9B", "Phylogroup K" = "#F6BB42", "Rubrivivax-like" = "#37BC9B", "Phylogroup I" = "#967ADC", "Phylogroup J" = "#D770AD", "Other" = "#434A54", "Unaffiliated" = "#AAB2BD", "None" = "#FFFFFF")) +
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


# PR analysis -------------------------------------------------------------

missing_samples <- tibble(label = c("OSD152_2014-06-20_5m_NPL022", "OSD15_2014-06-21_0m_NPL022", "OSD80_2014-06-21_2m_NPL022",
                                    "OSD105_2014-06-21_4.12m_NPL022", "OSD146_2014-06-21_5m_NPL022"), group = "None", counts = 1, prop = 1)

osd2014_pr_counts_prop <- osd2014_pr_counts %>%
  inner_join(osd2014_pr_groups %>% rename(supergroup = group, group = OTU)) %>%
  select(label, supergroup, counts) %>%
  dplyr::group_by(label, supergroup) %>%
  summarise(counts = sum(counts)) %>%
  ungroup() %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  droplevels() %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(prop = counts/sum(counts)) %>%
  ungroup() %>%
  arrange(label) %>%
  bind_rows(missing_samples) %>%
  complete(label, supergroup) %>%
  filter(!is.na(supergroup)) %>%
  dplyr::mutate(counts = ifelse(is.na(counts), 0 , counts), prop = ifelse(is.na(prop), 0 , prop), label = fct_relevel(label, st_100_order_terrestrial$label),
                supergroup = fct_relevel(supergroup, c("Marine SAR11-like clade (Alphaproteobacteria IV-8)",
                                                  "Freshwater SAR11-like clade (Alphaproteobacteria IV-8)",
                                                  "Other Alphaproteobacteria type IV",
                                                  "Rhodobacterales-like clade (Alphaproteobacteria III-3)",
                                                  "Other Alphaproteobacteria type III",
                                                  "Other proteorhodopsins",
                                                  "Betaproteobacteria type III (III-2; III3)",
                                                  "Gammaproteobacteria type I (I-1; I-2; I-3)",
                                                  "Gammaproteobacteria type III (III-5)",
                                                  "SAR86-like clade (Gammaproteobacteria type IV-5)")))


# Most AAnPB are members of the Alpha-, Beta- and Gammaproteobacteria. Yutin et al. (2007)
# proposed a phylogenetic framework for assigning AAnPB into 12 phylogroups
# (labelled phylogroups A through L), defined according to the structure of the puf operon
# and the phy- logeny of the pufM gene. Different phylogroups have been shown to be abundant
# in different marine environments. For example,

# Alphaproteobacterial phylogroups E and F,  which contain Rhodobacter-like bacteria, and G, which contains Roseobacter-like bacteria,
# are common in nutrient rich coastal waters (Salka et al., 2008; Lehours et al., 2010),

# phylogroups C and D have been mostly observed in offshore waters (Ritchie and Johnson, 2012).

# Phylogroup K, which contains members of the Gammaproteobacteria, also has a widespread distribution, #
# and has been demonstrated to dominate AAnPB composition in samples from the Atlantic and Pacific Oceans
# as well as the Arctic and Mediterranean Seas during summer (Lehours et al., 2010; Ferrera et al., 2014).

pr_plot <- ggplot(osd2014_pr_counts_prop, aes(label, prop, fill = supergroup)) +
  geom_col(color = "grey20", size = 0.2, width = 1) +
  scale_x_discrete(expand = c(0, 0), limits = st_100_order_terrestrial$label) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_fill_manual(values = c("Marine SAR11-like clade (Alphaproteobacteria IV-8)" = "#E9573F",
                               "Freshwater SAR11-like clade (Alphaproteobacteria IV-8)" = "#8CC152",
                               "Other Alphaproteobacteria type IV" =  "#4A89DC",
                               "Rhodobacterales-like clade (Alphaproteobacteria III-3)" = "#3BAFDA",
                               "Other proteorhodopsins" = "#434A54",
                               "Other Alphaproteobacteria type III" = "#F6BB42",
                               "Betaproteobacteria type III (III-2; III3)" = "#DA4453",
                               "Gammaproteobacteria type I (I-1; I-2; I-3)" = "#967ADC",
                               "Gammaproteobacteria type III (III-5)" = "#D770AD",
                               "SAR86-like clade (Gammaproteobacteria type IV-5)" = "#37BC9B")) +
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

ggpubr::ggarrange(pufm_plot, pr_plot, nrow = 2)
ggsave(plot = last_plot(), filename = "osd2014_shotgun/figures/osd2014_photoheterotrophic.pdf", width = 11.69, height = 8.27)
save.image(file = "osd2014_shotgun/data/osd2014_photoheterotrophic.Rda")
