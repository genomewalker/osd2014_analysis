library(phyloseq)
library(tidyverse)
library(broom)
library(Hmisc)
library(ggparl)
library(ggpol)
library(DESeq2)

# BEGIN: WARNING!!!! -------------------------------------------------------------
# You can access to the data used in this analysis in several ways:
# 1. You have a copy of the PostgreSQL DB
# 2. You downloaded the .Rdata files from http://osd2014.metagenomics.eu/ and placed them
#    in the data folder
# 3. You can load the files remotely, it might take a while when the file is very large
# END: WARNING!!!! -------------------------------------------------------------


# BEGIN: Load necessary data -----------------------------------------------------
# Use if you have the postgres DB in place
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)
osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)

# If downloaded file at osd2014_16S_asv/data/ use:
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_orig.Rdata", verbose = TRUE)
osd2014_dada2_phyloseq_emp <- osd2014_dada2_phyloseq_alpha
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata")
osd2014_dada2_phyloseq_alma <- osd2014_dada2_phyloseq_alpha
load("osd2014_16S_asv/data/osd2014_mitag_qiime97_physeq_filt_objects.Rdata", verbose = TRUE)
osd2014_phyloseq_mitag_qiime97 <- osd2014_mitag_qiime97_phyloseq_alpha

# Basic contextual data
load("osd2014_16S_asv/data/osd2014_basic_cdata.Rdata", verbose = TRUE)

# Data for comparing miTAGs, ASVs and functional unknown fraction
load("osd2014_shotgun/data/osd2014_eggnog_physeq_filt_objects.Rdata", verbose = TRUE)
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata", verbose = TRUE)
load("osd2014_16S_asv/data/osd2014_mitag_asv_physeq_filt_objects.Rdata", verbose = TRUE)
load("osd2014_shotgun/data/osd2014_unks_reads_physeq_filt_objects.Rdata", verbose = TRUE)
load("osd2014_shotgun/data/osd2014_unks_assm_physeq_filt_objects.Rdata", verbose = TRUE)
osd2014_simka_bc <- read_delim("osd2014_shotgun/data/allVSall_paired_k21/mat_abundance_braycurtis.csv.gz", col_names = TRUE, delim = ";") %>%
  column_to_rownames("X1") %>%
  as.matrix()

# If remote use
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_orig.Rdata"), verbose = TRUE)
osd2014_dada2_phyloseq_emp <- osd2014_dada2_phyloseq_alpha
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata"), verbose = TRUE)
osd2014_dada2_phyloseq_alma <- osd2014_dada2_phyloseq_alpha
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_mitag_qiime97_physeq_filt_objects.Rdata"), verbose = TRUE)
osd2014_phyloseq_mitag_qiime97 <- osd2014_mitag_qiime97_phyloseq_alpha


# Basic contextual data
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_basic_cdata.Rdata"), verbose = TRUE)

# Data for comparing miTAGs, ASVs and functional unknown fraction
load(url("http://osd2014.metagenomics.eu/osd2014_shotgun/data/osd2014_eggnog_physeq_filt_objects.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_mitag_asv_physeq_filt_objects.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_shotgun/data/osd2014_unks_reads_physeq_filt_objects.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_shotgun/data/osd2014_unks_assm_physeq_filt_objects.Rdata"), verbose = TRUE)
osd2014_simka_bc <- read_delim("http://osd2014.metagenomics.eu/osd2014_shotgun/data/allVSall_paired_k21/mat_abundance_braycurtis.csv.gz", col_names = TRUE, delim = ";") %>%
  column_to_rownames("X1") %>%
  as.matrix()
# END: Load necessary data -----------------------------------------------------


glom_prop_tidy <- function(X, Y){
  mycols <- c("Sample", "OTU", Y, "Abundance")
  mycols1 <- c("Sample", Y)
  com <- paste("!is.na(", Y, ")", sep = "")

  tax_X <- as(tax_table(X), "matrix") %>%
    as_tibble(rownames = "OTU")

  if (taxa_are_rows(X)){
    X <- t(as(otu_table(X),"matrix"))
  }else{
    X <- as(otu_table(X),"matrix")
  }

  X <- X %>%
    as_tibble(rownames = "Sample") %>%
    gather(OTU, Abundance, -Sample) %>%
    #filter(Abundance > 0) %>%
    inner_join(tax_X)
  X_glom <- X %>% as_tibble() %>% select(mycols) %>% filter_(.dots = com) %>% group_by_(.dots = mycols1) %>% summarise(n = sum(Abundance)) %>% ungroup
  X_glom_prop <- X_glom %>% group_by(Sample) %>% mutate(N= sum(n), prop = n/N)
  X_glom_prop
}

osd2014_dada2_phyloseq_emp_genus <- glom_prop_tidy(osd2014_dada2_phyloseq_emp, "Genus")
osd2014_dada2_phyloseq_alma_genus <- glom_prop_tidy(osd2014_dada2_phyloseq_alma, "Genus")
osd2014_dada2_phyloseq_mitag_qiime97_genus <- glom_prop_tidy(osd2014_phyloseq_mitag_qiime97, "Genus")

osd2014_dada2_phyloseq_emp_family <- glom_prop_tidy(osd2014_dada2_phyloseq_emp, "Family")
osd2014_dada2_phyloseq_alma_family <- glom_prop_tidy(osd2014_dada2_phyloseq_alma, "Family")
osd2014_dada2_phyloseq_mitag_qiime97_family <- glom_prop_tidy(osd2014_phyloseq_mitag_qiime97, "Family")

osd2014_dada2_phyloseq_emp_order <- glom_prop_tidy(osd2014_dada2_phyloseq_emp, "Order")
osd2014_dada2_phyloseq_alma_order <- glom_prop_tidy(osd2014_dada2_phyloseq_alma, "Order")
osd2014_dada2_phyloseq_mitag_qiime97_order <- glom_prop_tidy(osd2014_phyloseq_mitag_qiime97, "Order")

osd2014_dada2_phyloseq_emp_class <- glom_prop_tidy(osd2014_dada2_phyloseq_emp, "Class")
osd2014_dada2_phyloseq_alma_class <- glom_prop_tidy(osd2014_dada2_phyloseq_alma, "Class")
osd2014_dada2_phyloseq_mitag_qiime97_class <- glom_prop_tidy(osd2014_phyloseq_mitag_qiime97, "Class")

osd2014_dada2_phyloseq_emp_phylum <- glom_prop_tidy(osd2014_dada2_phyloseq_emp, "Phylum")
osd2014_dada2_phyloseq_alma_phylum <- glom_prop_tidy(osd2014_dada2_phyloseq_alma, "Phylum")
osd2014_dada2_phyloseq_mitag_qiime97_phylum <- glom_prop_tidy(osd2014_phyloseq_mitag_qiime97, "Phylum")

# save(osd2014_dada2_phyloseq_emp_genus,
#      osd2014_dada2_phyloseq_alma_genus,
#      osd2014_dada2_phyloseq_mitag_qiime97_genus,
#
#      osd2014_dada2_phyloseq_emp_family,
#      osd2014_dada2_phyloseq_alma_family,
#      osd2014_dada2_phyloseq_mitag_qiime97_family,
#
#      osd2014_dada2_phyloseq_emp_order,
#      osd2014_dada2_phyloseq_alma_order,
#      osd2014_dada2_phyloseq_mitag_qiime97_order,
#
#      osd2014_dada2_phyloseq_emp_class,
#      osd2014_dada2_phyloseq_alma_class,
#      osd2014_dada2_phyloseq_mitag_qiime97_class,
#
#      osd2014_dada2_phyloseq_emp_phylum,
#      osd2014_dada2_phyloseq_alma_phylum,
#      osd2014_dada2_phyloseq_mitag_qiime97_phylum, file = "osd2014_16S_asv/data/osd2014_16S_mitag_qiime97_emp_alma_comp_objects.Rdata")


pair_cor <- function(X,Y){
  X %>% ungroup() %>% select(-n, -N,) %>% dplyr::rename(prop_x = prop) %>%
    inner_join(Y %>% ungroup() %>% select(-n, -N,) %>% dplyr::rename(prop_y = prop)) %>%
    group_by(Sample) %>%
    do(Hmisc::rcorr(.$prop_x, .$prop_y, "spearman") %>% broom::tidy()) %>% ungroup()
}

osd2014_emp_vs_alma_genus <- pair_cor(osd2014_dada2_phyloseq_emp_genus, osd2014_dada2_phyloseq_alma_genus) %>% mutate(comp = "ALMA vs EMP", rank = "Genus")
osd2014_emp_vs_mitag_qiime97_genus <- pair_cor(osd2014_dada2_phyloseq_emp_genus, osd2014_dada2_phyloseq_mitag_qiime97_genus) %>% mutate(comp = "EMP vs miTAG", rank = "Genus")
osd2014_alma_vs_mitag_qiime97_genus <- pair_cor(osd2014_dada2_phyloseq_alma_genus, osd2014_dada2_phyloseq_mitag_qiime97_genus) %>% mutate(comp = "ALMA vs miTAG", rank = "Genus")

colors <- c("#EA7580", "#F8CD9C", "#088BBE")
p_genus <- bind_rows(osd2014_emp_vs_alma_genus, osd2014_emp_vs_mitag_qiime97_genus, osd2014_alma_vs_mitag_qiime97_genus) %>%
  mutate(comp = fct_relevel(comp, c("ALMA vs EMP", "ALMA vs miTAG", "EMP vs miTAG"))) %>%
  ggplot(aes(comp, y = estimate, fill = comp), color = "black") +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = "black", width = 0.5,
                 jitter.height = 0.05, jitter.width = 0.075, errorbar.draw = TRUE, jitter.colour = "black", alpha = 0.7) +
  theme(legend.position = "none") +
  ylab(expression(Spearman~rho)) +
  xlab("") +
  theme_bw() + ylim(c(0.25,1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colors) +
  ggtitle("Genus")

osd2014_emp_vs_alma_family <- pair_cor(osd2014_dada2_phyloseq_emp_family, osd2014_dada2_phyloseq_alma_family) %>% mutate(comp = "ALMA vs EMP", rank = "Family")
osd2014_emp_vs_mitag_qiime97_family <- pair_cor(osd2014_dada2_phyloseq_emp_family, osd2014_dada2_phyloseq_mitag_qiime97_family) %>% mutate(comp = "EMP vs miTAG", rank = "Family")
osd2014_alma_vs_mitag_qiime97_family <- pair_cor(osd2014_dada2_phyloseq_alma_family, osd2014_dada2_phyloseq_mitag_qiime97_family) %>% mutate(comp = "ALMA vs miTAG", rank = "Family")

p_family <- bind_rows(osd2014_emp_vs_alma_family, osd2014_emp_vs_mitag_qiime97_family, osd2014_alma_vs_mitag_qiime97_family) %>%
  mutate(comp = fct_relevel(comp, c("ALMA vs EMP", "ALMA vs miTAG", "EMP vs miTAG"))) %>%
  ggplot(aes(comp, y = estimate, fill = comp), color = "black") +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = "black", width = 0.5,
                 jitter.height = 0.05, jitter.width = 0.075, errorbar.draw = TRUE, jitter.colour = "black", alpha = 0.7) +
  theme(legend.position = "none") +
  ylab(expression(Spearman~rho)) +
  xlab("") +
  theme_bw() + ylim(c(0.25,1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colors) +
  ggtitle("Family")

osd2014_emp_vs_alma_order <- pair_cor(osd2014_dada2_phyloseq_emp_order, osd2014_dada2_phyloseq_alma_order) %>% mutate(comp = "ALMA vs EMP", rank = "Order")
osd2014_emp_vs_mitag_qiime97_order <- pair_cor(osd2014_dada2_phyloseq_emp_order, osd2014_dada2_phyloseq_mitag_qiime97_order) %>% mutate(comp = "EMP vs miTAG", rank = "Order")
osd2014_alma_vs_mitag_qiime97_order <- pair_cor(osd2014_dada2_phyloseq_alma_order, osd2014_dada2_phyloseq_mitag_qiime97_order) %>% mutate(comp = "ALMA vs miTAG", rank = "Order")

p_order <- bind_rows(osd2014_emp_vs_alma_order, osd2014_emp_vs_mitag_qiime97_order, osd2014_alma_vs_mitag_qiime97_order) %>%
  mutate(comp = fct_relevel(comp, c("ALMA vs EMP", "ALMA vs miTAG", "EMP vs miTAG"))) %>%
  ggplot(aes(comp, y = estimate, fill = comp), color = "black") +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = "black", width = 0.5,
                 jitter.height = 0.05, jitter.width = 0.075, errorbar.draw = TRUE, jitter.colour = "black", alpha = 0.7) +
  theme(legend.position = "none") +
  ylab(expression(Spearman~rho)) +
  xlab("") +
  theme_bw() + ylim(c(0.25,1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colors) +
  ggtitle("Order")


osd2014_emp_vs_alma_class <- pair_cor(osd2014_dada2_phyloseq_emp_class, osd2014_dada2_phyloseq_alma_class) %>% mutate(comp = "ALMA vs EMP", rank = "Class")
osd2014_emp_vs_mitag_qiime97_class <- pair_cor(osd2014_dada2_phyloseq_emp_class, osd2014_dada2_phyloseq_mitag_qiime97_class) %>% mutate(comp = "EMP vs miTAG", rank = "Class")
osd2014_alma_vs_mitag_qiime97_class <- pair_cor(osd2014_dada2_phyloseq_alma_class, osd2014_dada2_phyloseq_mitag_qiime97_class) %>% mutate(comp = "ALMA vs miTAG", rank = "Class")

p_class <- bind_rows(osd2014_emp_vs_alma_class, osd2014_emp_vs_mitag_qiime97_class, osd2014_alma_vs_mitag_qiime97_class) %>%
  mutate(comp = fct_relevel(comp, c("ALMA vs EMP", "ALMA vs miTAG", "EMP vs miTAG"))) %>%
  ggplot(aes(comp, y = estimate, fill = comp), color = "black") +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = "black", width = 0.5,
                 jitter.height = 0.05, jitter.width = 0.075, errorbar.draw = TRUE, jitter.colour = "black", alpha = 0.7) +
  theme(legend.position = "none") +
  ylab(expression(Spearman~rho)) +
  xlab("") +
  theme_bw() + ylim(c(0.25,1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colors) +
  ggtitle("Class")


osd2014_emp_vs_alma_phylum <- pair_cor(osd2014_dada2_phyloseq_emp_phylum, osd2014_dada2_phyloseq_alma_phylum) %>% mutate(comp = "ALMA vs EMP", rank = "Phylum")
osd2014_emp_vs_mitag_qiime97_phylum <- pair_cor(osd2014_dada2_phyloseq_emp_phylum, osd2014_dada2_phyloseq_mitag_qiime97_phylum) %>% mutate(comp = "EMP vs miTAG", rank = "Phylum")
osd2014_alma_vs_mitag_qiime97_phylum <- pair_cor(osd2014_dada2_phyloseq_alma_phylum, osd2014_dada2_phyloseq_mitag_qiime97_phylum) %>% mutate(comp = "ALMA vs miTAG", rank = "Phylum")

p_phylum <- bind_rows(osd2014_emp_vs_alma_phylum, osd2014_emp_vs_mitag_qiime97_phylum, osd2014_alma_vs_mitag_qiime97_phylum) %>%
  mutate(comp = fct_relevel(comp, c("ALMA vs EMP", "ALMA vs miTAG", "EMP vs miTAG"))) %>%
  ggplot(aes(comp, y = estimate, fill = comp), color = "black") +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = "black", width = 0.5,
                 jitter.height = 0.05, jitter.width = 0.075, errorbar.draw = TRUE, jitter.colour = "black", alpha = 0.7) +
  theme(legend.position = "none") +
  ylab(expression(Spearman~rho)) +
  xlab("") +
  theme_bw() + ylim(c(0.25,1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colors) +
  ggtitle("Phylum")


ggpubr::ggarrange(p_genus, p_family, p_order, p_class, p_phylum, ncol = 5, legend = "top", common.legend = TRUE, labels = "AUTO")
ggsave(plot = last_plot(), filename = "osd2014_16S_asv/figures/osd2014_mitag_amplicon_comp_cor.pdf", width = 12, height = 3)


#
# osd2014_dada2_phyloseq_emp_prop_genus <- tax_glom(osd2014_dada2_phyloseq_emp_prop, taxrank = "Genus")
# osd2014_dada2_phyloseq_alma_prop_genus <- tax_glom(osd2014_dada2_phyloseq_alma_prop, taxrank = "Genus")
#
# osd2014_dada2_phyloseq_emp_prop_genus_long <- psmelt(osd2014_dada2_phyloseq_emp_prop_genus) %>%
#   as_tibble() %>%
#   select(Sample, Genus, Abundance) %>%
#   rename(prop_emp = Abundance)
#
# osd2014_dada2_phyloseq_alma_prop_genus_long <- psmelt(osd2014_dada2_phyloseq_alma_prop_genus) %>%
#   as_tibble() %>%
#   select(Sample, Genus, Abundance) %>%
#   rename(prop_alma = Abundance)
#
# osd2014_dada2_phyloseq_emp_prop_genus_long %>%
#   inner_join(osd2014_dada2_phyloseq_alma_prop_genus_long) %>%
#   ungroup() %>%
#   group_by(Sample) %>%
#   do(rcorr(.$prop_alma, .$prop_emp, "spearman") %>% broom::tidy()) %>%
#   ungroup() %>%
#   ggplot(aes("EMP", y = estimate)) +
#   geom_flat_violin(scale = "count", trim = FALSE, width = 0.5) +
#   stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
#                geom = "pointrange", position = position_nudge(0.05)) +
#   geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "down",
#                position = position_nudge(-0.025)) +
#   theme(legend.position = "none") +
#   ylab(expression(Spearman~rho)) +
#   xlab("") +
#   theme_light() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
#
# osd2014_dada2_phyloseq_emp_prop_genus_long %>%
#   inner_join(osd2014_dada2_phyloseq_alma_prop_genus_long) %>%
#   ungroup() %>%
#   group_by(Sample) %>%
#   do(rcorr(.$prop_alma, .$prop_emp, "spearman") %>% broom::tidy()) %>%
#   ungroup() %>%
#   ggplot(aes("EMP", y = estimate)) +
#   geom_boxjitter(outlier.color = NA, jitter.shape = 21, errorbar.draw = TRUE, width = 0.2, jitter.colour = "black", fill = "#103C54", alpha = 0.7) +
#   theme(legend.position = "none") +
#   ylab(expression(Spearman~rho)) +
#   xlab("") +
#   theme_light()

# osd_mitag <- read_tsv("~/Downloads/osd2014_mitag_qiime97_tax.tsv", col_names = T)
# osd_mitag$label <- plyr::mapvalues(osd_mitag$label,
#                   c('OSD114_2014-06-21_50m_NPL022', 'OSD118_2014-06-24_0.2m_NPL022', 'OSD72_2014-07-21_0.8m_NPL022', 'OSD159_2014-06-21_2m_NPL022'),
#                   c('OSD114_2014-06-21_1m_NPL022', 'OSD118_2014-06-21_0.2m_NPL022', 'OSD72_2014-06-21_0.8m_NPL022', 'OSD159_2014-06-19_2m_NPL022'))
#
# osd_mitag <- osd_mitag %>% filter(label %in% sample_names(osd2014_dada2_phyloseq_alma))
#
# osd2014_mitag <- osd_mitag %>% select(seq_id, label) %>% group_by(seq_id, label) %>% count %>%
#   spread(label, value = n, fill = 0) %>%
#   column_to_rownames("seq_id") %>%
#   as.matrix
#
# osd2014_mitag_qiime97_tax <- osd_mitag %>%
#   filter(seq_id %in% rownames(osd2014_mitag)) %>%
#   select(-asv, -label) %>%
#   unique() %>%
#   column_to_rownames("seq_id") %>%
#   as.matrix
#
# osd2014_mitag_qiime97_phyloseq <- phyloseq(otu_table(osd2014_mitag, taxa_are_rows = TRUE), tax_table(osd2014_mitag_qiime97_tax))
#
# save(osd2014_mitag_qiime97_phyloseq, file = "osd2014_16S_asv/data/osd2014_mitag_qiime97_phyloseq.Rdata")



# Which ASVs are common or different between EMP and ALMA -----------------
# library(RPostgreSQL)
# my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
# osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
#   collect(n = Inf)
# st_100_order_terrestrial <- tbl(my_db, "osd2014_st_order_terrestrial") %>%
#   collect(n = Inf) %>%
#   filter(label %in% osd2014_amp_mg_intersect$label)
#
# osd2014_silva_dada2_emp <- data.frame(as(tax_table(osd2014_dada2_phyloseq_emp), "matrix")) %>% rownames_to_column("asv") %>% as_tibble()
# osd2014_silva_dada2_alma <- data.frame(as(tax_table(osd2014_dada2_phyloseq_alma), "matrix")) %>% rownames_to_column("asv") %>% as_tibble()
#
# write_tsv(osd2014_silva_dada2_emp %>% select(asv_name, asv), path = "osd2014_16S_asv/data/osd2014_16S_asv_sequences_emp.tsv", col_names = FALSE)
# write_tsv(osd2014_silva_dada2_alma %>% select(asv_name, asv), path = "osd2014_16S_asv/data/osd2014_16S_asv_sequences_alma.tsv", col_names = FALSE)
#
# system("awk '{print $1\"\t\"$2}' osd2014_16S_asv/data/osd2014_16S_asv_sequences_emp.tsv | seqkit tab2fx > osd2014_16S_asv/data/osd2014_16S_asv_sequences_emp.fasta")
# system("awk '{print $1\"\t\"$2}' osd2014_16S_asv/data/osd2014_16S_asv_sequences_alma.tsv | seqkit tab2fx > osd2014_16S_asv/data/osd2014_16S_asv_sequences_alma.fasta")
#
# system("vsearch --usearch_global osd2014_16S_asv/data/osd2014_16S_asv_sequences_emp.fasta --db osd2014_16S_asv/data/osd2014_16S_asv_sequences_alma.fasta --blast6out osd2014_16S_asv/data/osd2014_16S_asv_sequences_emp_vs_alma.tsv --id 0.1")
# system("vsearch --usearch_global osd2014_16S_asv/data/osd2014_16S_asv_sequences_alma.fasta --db osd2014_16S_asv/data/osd2014_16S_asv_sequences_emp.fasta --blast6out osd2014_16S_asv/data/osd2014_16S_asv_sequences_alma_vs_emp.tsv --id 0.1")
#
# emp_vs_alma <- read_tsv("osd2014_16S_asv/data/osd2014_16S_asv_sequences_emp_vs_alma.tsv", col_names = FALSE)
# alma_vs_emp <- read_tsv("osd2014_16S_asv/data/osd2014_16S_asv_sequences_alma_vs_emp.tsv", col_names = FALSE)
#
# c(emp_vs_alma %>% filter(X3 == 100) %>% .$X1,
#   alma_vs_emp %>% filter(X3 == 100) %>% .$X2) %>% unique() %>% length()
#
# c(emp_vs_alma %>% filter(X3 == 100) %>% .$X2,
#   alma_vs_emp %>% filter(X3 == 100) %>% .$X1) %>% unique() %>% length()
#
# osd2014_dada2_phyloseq_emp_prop <- transform_sample_counts(osd2014_dada2_phyloseq_emp, function (x) x/sum(x))
# osd2014_dada2_phyloseq_alma_prop <- transform_sample_counts(osd2014_dada2_phyloseq_alma, function (x) x/sum(x))
#
#
# osd2014_dada2_phyloseq_emp_spec <- subset_taxa(osd2014_dada2_phyloseq_emp_prop, !(asv_name %in% (emp_vs_alma %>% filter(X3 == 100) %>% .$X1)))
# osd2014_dada2_phyloseq_emp_spec_l <- psmelt(osd2014_dada2_phyloseq_emp_spec) %>%
#   as_tibble() %>%
#   mutate(Sample = fct_relevel(Sample, st_100_order_terrestrial$label))
#
# osd2014_dada2_phyloseq_alma_spec <- subset_taxa(osd2014_dada2_phyloseq_alma_prop, !(asv_name %in% (alma_vs_emp %>% filter(X3 == 100) %>% .$X1)))
# osd2014_dada2_phyloseq_alma_spec_l <- psmelt(osd2014_dada2_phyloseq_alma_spec) %>%
#   as_tibble() %>%
#   mutate(Sample = fct_relevel(Sample, st_100_order_terrestrial$label))
#
#
# emp_excl_asv_family <- osd2014_dada2_phyloseq_emp_spec_l %>%
#   select(Kingdom, Phylum, Class, Order, Family, Genus, Species, asv_name) %>%
#   unique() %>%
#   count(Family, sort = TRUE) %>%
#   mutate(Family = fct_explicit_na(Family, "Not assigned"), Family = fct_reorder(Family, n)) %>%
#   head(25) %>%
#   ggplot(aes(Family, n)) +
#   geom_col() +
#   ggpubr::rotate() +
#   theme_light()+
#   ylab("Number of ASVs")
#
# alma_excl_asv_family <- osd2014_dada2_phyloseq_alma_spec_l %>%
#   select(Kingdom, Phylum, Class, Order, Family, Genus, Species, asv_name) %>%
#   unique() %>%
#   count(Family, sort = TRUE) %>%
#   mutate(Family = fct_explicit_na(Family, "Not assigned"), Family = fct_reorder(Family, n)) %>%
#   head(25) %>%
#   ggplot(aes(Family, n)) +
#   geom_col() +
#   ggpubr::rotate() +
#   theme_light()+
#   ylab("Number of ASVs")
#
# emp_alma_excl_asv_family <- ggpubr::ggarrange(emp_excl_asv_family, alma_excl_asv_family, labels = "AUTO")
# ggsave(emp_alma_excl_asv_family, filename = "osd2014_16S_asv/figures/emp_alma_excl_asv_family.pdf", width = 11.69, height = 8.27)
#
#
# emp_excl_asv_order <- osd2014_dada2_phyloseq_emp_spec_l %>%
#   select(Kingdom, Phylum, Class, Order, Family, Genus, Species, asv_name) %>%
#   unique() %>%
#   count(Order, sort = TRUE) %>%
#   mutate(Order = fct_explicit_na(Order, "Not assigned"), Order = fct_reorder(Order, n)) %>%
#   head(25) %>%
#   ggplot(aes(Order, n)) +
#   geom_col() +
#   ggpubr::rotate() +
#   theme_light()+
#   ylab("Number of ASVs")
#
# alma_excl_asv_order <- osd2014_dada2_phyloseq_alma_spec_l %>%
#   select(Kingdom, Phylum, Class, Order, Family, Genus, Species, asv_name) %>%
#   unique() %>%
#   count(Order, sort = TRUE) %>%
#   mutate(Order = fct_explicit_na(Order, "Not assigned"), Order = fct_reorder(Order, n)) %>%
#   head(25) %>%
#   ggplot(aes(Order, n)) +
#   geom_col() +
#   ggpubr::rotate() +
#   theme_light()+
#   ylab("Number of ASVs")
#
# emp_alma_excl_asv_order <- ggpubr::ggarrange(emp_excl_asv_order, alma_excl_asv_order, labels = "AUTO")
# ggsave(emp_alma_excl_asv_order, filename = "osd2014_16S_asv/figures/emp_alma_excl_asv_order.pdf", width = 11.69, height = 8.27)
#
#
#
# osd2014_dada2_phyloseq_emp_shar <- subset_taxa(osd2014_dada2_phyloseq_emp_prop, (asv_name %in% (emp_vs_alma %>% filter(X3 == 100) %>% .$X1)))
# osd2014_dada2_phyloseq_emp_shar_l <- psmelt(osd2014_dada2_phyloseq_emp_shar) %>%
#   as_tibble() %>%
#   mutate(Sample = fct_relevel(Sample, st_100_order_terrestrial$label))
#
# osd2014_dada2_phyloseq_alma_shar <- subset_taxa(osd2014_dada2_phyloseq_alma_prop, (asv_name %in% (alma_vs_emp %>% filter(X3 == 100) %>% .$X1)))
# osd2014_dada2_phyloseq_alma_shar_l <- psmelt(osd2014_dada2_phyloseq_alma_shar) %>%
#   as_tibble() %>%
#   mutate(Sample = fct_relevel(Sample, st_100_order_terrestrial$label))
#
#
# #modified ggplot
# ggplot(osd2014_dada2_phyloseq_emp_spec_l, aes(x = Sample, y = Abundance)) +
#   geom_col() +
#   facet_wrap(~Phylum) +
#   scale_y_continuous(labels = scales::percent)
#
# ggplot(osd2014_dada2_phyloseq_alma_spec_l, aes(x = Sample, y = Abundance)) +
#   geom_col() +
#   facet_wrap(~Phylum) +
#   scale_y_continuous(labels = scales::percent)
#
# ggplot(osd2014_dada2_phyloseq_emp_shar_l, aes(x = Sample, y = Abundance)) +
#   geom_col() +
#   facet_wrap(~Phylum) +
#   scale_y_continuous(labels = scales::percent)
#
# ggplot(osd2014_dada2_phyloseq_alma_shar_l, aes(x = Sample, y = Abundance)) +
#   geom_col() +
#   facet_wrap(~Phylum) +
#   scale_y_continuous(labels = scales::percent)
#
#
# gheatmap <- function(X, Y){
#   X <- X %>%
#     select(Sample, Phylum, Abundance) %>%
#     filter(Abundance > 0) %>%
#     group_by(Sample, Phylum) %>%
#     summarise(Abundance = sum(Abundance), n = n()) %>%
#     ungroup()
#
#   X1_euc <- X %>%
#     select(Sample, Phylum, Abundance) %>%
#     spread(Sample, Abundance, fill = 0) %>%
#     as.data.frame() %>%
#     column_to_rownames(var = "Phylum") %>%
#     dist()
#   X1_clu <- hclust(X1_euc, method = "ward.D2")
#
#   X2_euc <- X %>%
#     select(Sample, Phylum, Abundance) %>%
#     spread(Phylum, Abundance, fill = 0) %>%
#     as.data.frame() %>%
#     column_to_rownames(var = "Sample") %>%
#     dist()
#   X2_clu <- hclust(X2_euc, method = "ward.D2")
#
#
#   X <- X %>% mutate(Phylum = fct_relevel(Phylum, X1_clu$labels[X1_clu$order]), Sample = fct_relevel(Sample, X2_clu$labels[X2_clu$order]))
#   base_breaks <- function(n = 10){
#     function(x) {
#       axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
#     }
#   }
#   textcol <- "grey40"
#   #modified ggplot
#   ggplot(X, aes(x = Sample, y = Phylum, fill = 100 * (Abundance + 0.0000001), size = n))+
#     #geom_point(shape = 22, fill = "white", color = "white", size = 4) +
#     geom_point(shape = 22) +
#     #redrawing tiles to remove cross lines from legend
#     #geom_tile(colour = "white", size = 0.25, show.legend = FALSE)+
#     #remove axis labels, add title
#     labs(x="",y="",title="")+
#     #remove extra space
#     #scale_y_discrete(expand=c(0,0))+
#     #custom breaks on x-axis
#     #scale_x_discrete(expand=c(0,0))+
#     scale_size_continuous( trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum, range = c(2,4)) +
#     #custom colours for cut levels and na values
#     #scale_fill_gradientn(colours =rev(c("#d53e4f","#f46d43","#fdae61",
#     # "#fee08b","#e6f598","#abdda4","#ddf1da")), na.value="grey90", trans = "sqrt", labels = scales::percent) +
#     viridis::scale_fill_viridis(option = "D", trans = scales::sqrt_trans(), breaks = base_breaks(), labels = prettyNum) +
#     #scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(7,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
#     #mark year of vaccination
#     #geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
#     #equal aspect ratio x and y axis
#     coord_fixed() +
#     #set base size for all font elements
#     theme_grey(base_size=10)+
#     #theme options
#     theme(
#       legend.position = "bottom",
#       #remove legend title
#       legend.title=element_blank(),
#       #remove legend margin
#       legend.spacing = grid::unit(0,"cm"),
#       #change legend text properties
#       legend.text=element_text(colour=textcol,size=7,face="bold"),
#       #change legend key height
#       legend.key.height=grid::unit(0.2,"cm"),
#       #set a slim legend
#       legend.key.width=grid::unit(0.8,"cm"),
#       #set x axis text size and colour
#       axis.text.x=element_text(hjust = 1, vjust = 0.5, colour=textcol, angle = 90, size = 6),
#       #axis.ticks.y = element_blank(),
#       #axis.text.y = element_blank(),
#       #set y axis text colour and adjust vertical justification
#       #axis.text.y=element_text(vjust = 0.2,colour=textcol, size = 6),
#       #change axis ticks thickness
#       axis.ticks=element_line(size=0.4),
#       #change title font, size, colour and justification
#       plot.title=element_blank(),
#       #remove plot background
#       plot.background=element_blank(),
#       #remove plot border
#       panel.border=element_blank())
# }
# gheatmap(osd2014_dada2_phyloseq_emp_spec_l)
# gheatmap(osd2014_dada2_phyloseq_alma_spec_l)
# gheatmap(osd2014_dada2_phyloseq_emp_shar_l)
# gheatmap(osd2014_dada2_phyloseq_alma_shar_l)


# miTAGs vs simka dissimilarities -----------------------------------------
#
# bc_mt <- read.table(file = "~/Downloads/mat_abundance_braycurtis-norm.csv", header = T, row.names = 1, sep = ";", check.names = F)
#
# # fix names
#
# colnames(bc_mt) <- plyr::mapvalues(colnames(bc_mt),
#                                    c('OSD114_2014-06-21_50m_NPL022', 'OSD118_2014-06-24_0.2m_NPL022', 'OSD72_2014-07-21_0.8m_NPL022', 'OSD159_2014-06-21_2m_NPL022'),
#                                    c('OSD114_2014-06-21_1m_NPL022', 'OSD118_2014-06-21_0.2m_NPL022', 'OSD72_2014-06-21_0.8m_NPL022', 'OSD159_2014-06-19_2m_NPL022'))
# rownames(bc_mt) <- plyr::mapvalues(rownames(bc_mt),
#                                    c('OSD114_2014-06-21_50m_NPL022', 'OSD118_2014-06-24_0.2m_NPL022', 'OSD72_2014-07-21_0.8m_NPL022', 'OSD159_2014-06-21_2m_NPL022'),
#                                    c('OSD114_2014-06-21_1m_NPL022', 'OSD118_2014-06-21_0.2m_NPL022', 'OSD72_2014-06-21_0.8m_NPL022', 'OSD159_2014-06-19_2m_NPL022'))
#
#
# bc_mt <- bc_mt[osd2014_amp_mg_intersect$label, osd2014_amp_mg_intersect$label]
#
#
# bc_mt_long <- broom::tidy(as.dist(bc_mt)) %>%
#   rename()

#library(RPostgreSQL)  # loads the PostgreSQL driver
#drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
#con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
#dbWriteTable(con, c("osd_analysis", "osd2014_simka_bc"), value=bc_mt_long,overwrite=TRUE,row.names=FALSE)



# Comparing miTAGs, ASVs and unknowns -------------------------------------


osd2014_dada2_phyloseq_alpha_genus <-  tax_glom(osd2014_dada2_phyloseq_alpha, taxrank = "Genus")

osd2014_dada2_phyloseq_alpha_genus_deseq <- phyloseq_to_deseq2(osd2014_dada2_phyloseq_alpha_genus, ~1)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(osd2014_dada2_phyloseq_alpha_genus_deseq), 1, gm_mean)
diagdds <- estimateSizeFactors(osd2014_dada2_phyloseq_alpha_genus_deseq, type = "poscounts")
diagdds <- estimateDispersions(diagdds)
diagvst <- varianceStabilizingTransformation(diagdds, blind = FALSE)
diagvst <- assay(diagvst)
diagvst[diagvst < 0] <- 0
osd2014_dada2_phyloseq_alpha_genus_vst <- osd2014_dada2_phyloseq_alpha_genus
osd2014_dada2_phyloseq_alpha_genus_norm <- osd2014_dada2_phyloseq_alpha_genus
otu_table(osd2014_dada2_phyloseq_alpha_genus_vst) <- otu_table(diagvst, taxa_are_rows = TRUE)
otu_table(osd2014_dada2_phyloseq_alpha_genus_norm) <- otu_table(counts(diagdds, normalized = TRUE), taxa_are_rows = TRUE)

osd2014_dada2_phyloseq_alpha_genus_prop <- transform_sample_counts(osd2014_dada2_phyloseq_alpha_genus, function(x) x/sum(x))
mean_prop <- taxa_sums(osd2014_dada2_phyloseq_alpha_genus_prop)/nsamples(osd2014_dada2_phyloseq_alpha_genus)
osd2014_dada2_phyloseq_beta_vst_genus <- prune_taxa(mean_prop > 1e-5, osd2014_dada2_phyloseq_alpha_genus_vst)

osd2014_mitag_phyloseq_beta_vst_bc <- phyloseq::distance(osd2014_mitag_phyloseq_beta_css, "bray")
osd2014_dada2_phyloseq_beta_vst_bc <- phyloseq::distance(osd2014_dada2_phyloseq_beta_css, "bray")
osd2014_dada2_phyloseq_beta_vst_genus_bc <- phyloseq::distance(osd2014_dada2_phyloseq_beta_vst_genus, "bray")
osd2014_eggnog_phyloseq_beta_vst_bc <- phyloseq::distance(osd2014_eggnog_phyloseq_beta_css, "bray")
osd2014_unks_reads_phyloseq_beta_vst_bc <- phyloseq::distance(osd2014_unks_reads_phyloseq_beta_css, "bray")
osd2014_unks_assm_phyloseq_beta_vst_bc <- phyloseq::distance(osd2014_unks_assm_phyloseq_beta_css, "bray")


mitag_bc <- broom::tidy(as.dist(as.matrix(osd2014_mitag_phyloseq_beta_vst_bc)[osd2014_amp_mg_intersect$label,osd2014_amp_mg_intersect$label])) %>% as_tibble()
asv_bc <- broom::tidy(as.dist(as.matrix(osd2014_dada2_phyloseq_beta_vst_bc)[osd2014_amp_mg_intersect$label,osd2014_amp_mg_intersect$label])) %>% as_tibble()
asv_genus_bc <- broom::tidy(as.dist(as.matrix(osd2014_dada2_phyloseq_beta_vst_genus_bc)[osd2014_amp_mg_intersect$label,osd2014_amp_mg_intersect$label])) %>% as_tibble()
eggnog_bc <- broom::tidy(as.dist(as.matrix(osd2014_eggnog_phyloseq_beta_vst_bc)[osd2014_amp_mg_intersect$label,osd2014_amp_mg_intersect$label])) %>% as_tibble()
unk_reads_bc <- broom::tidy(as.dist(as.matrix(osd2014_unks_reads_phyloseq_beta_vst_bc)[osd2014_amp_mg_intersect$label,osd2014_amp_mg_intersect$label])) %>% as_tibble()
unk_assm_bc <- broom::tidy(as.dist(as.matrix(osd2014_unks_assm_phyloseq_beta_vst_bc)[osd2014_amp_mg_intersect$label,osd2014_amp_mg_intersect$label])) %>% as_tibble()

simka_bc <- broom::tidy(as.dist(osd2014_simka_bc[osd2014_amp_mg_intersect$label,osd2014_amp_mg_intersect$label]))



mitag_bc %>%
  left_join(simka_bc, by = c("item1", "item2")) %>%
  ggplot(aes(1-distance.y, 1-distance.x)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  ylab("Similarity in community composition") +
  xlab("Similarity in community functional attributes") +
  theme_light()

eggnog_bc %>%
  left_join(simka_bc, by = c("item1", "item2")) %>%
  ggplot(aes(1-distance.y, 1-distance.x)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  ylab("Similarity in community composition") +
  xlab("Similarity in community functional attributes") +
  theme_light()

asv_bc %>%
  left_join(simka_bc, by = c("item1", "item2")) %>%
  ggplot(aes(1-distance.y, 1-distance.x)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  ylab("Similarity in community composition") +
  xlab("Similarity in community functional attributes") +
  theme_light()

asv_genus_bc %>%
  left_join(simka_bc, by = c("item1", "item2")) %>%
  ggplot(aes(1-distance.y, 1-distance.x)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  ylab("Similarity in community composition") +
  xlab("Similarity in community functional attributes") +
  theme_light()

asv_bc %>%
  left_join(eggnog_bc, by = c("item1", "item2")) %>%
  ggplot(aes(1-distance.y, 1-distance.x)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  ylab("Similarity in community composition") +
  xlab("Similarity in community functional attributes") +
  theme_light()

asv_bc %>%
  left_join(eggnog_bc, by = c("item1", "item2")) %>%
  ggplot(aes(1-distance.y, 1-distance.x)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  ylab("Similarity in community composition") +
  xlab("Similarity in community functional attributes") +
  theme_light()


 simka_bc %>%
  left_join(unk_reads_bc, by = c("item1", "item2")) %>%
  dplyr::rename(unk_assm = distance.x, simka = distance.y) %>%
  ggplot(aes(1-unk_assm, 1-simka)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Similarity in community functional attributes (simka)") +
  ylab("Similarity in community functional attributes (unks)") +
  theme_light()


unk_reads_bc %>%
  left_join(unk_assm_bc, by = c("item1", "item2")) %>%
  dplyr::rename(unk_assm = distance.x, simka = distance.y) %>%
  ggplot(aes(1-unk_assm, 1-simka)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Similarity in community functional attributes (unks)") +
  ylab("Similarity in community composition") +
  theme_light()

unk_reads_bc %>%
  left_join(mitag_bc, by = c("item1", "item2")) %>%
  dplyr::rename(unk_assm = distance.x, simka = distance.y) %>%
  ggplot(aes(1-unk_assm, 1-simka)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Similarity in community functional attributes (simka)") +
  ylab("Similarity in community composition") +
  theme_light()


eggnog_bc %>%
  left_join(mitag_bc, by = c("item1", "item2")) %>%
  dplyr::rename(unk_assm = distance.x, simka = distance.y) %>%
  ggplot(aes(1-unk_assm, 1-simka)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Similarity in community functional attributes (eggnog)") +
  ylab("Similarity in community composition") +
  theme_light()


simka_bc %>%
  left_join(mitag_bc, by = c("item1", "item2")) %>%
  dplyr::rename(unk_assm = distance.x, simka = distance.y) %>%
  ggplot(aes(1-unk_assm, 1-simka)) +
  geom_point(shape = 21, fill = "grey", alpha = 0.7) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Similarity in community functional attributes (simka)") +
  ylab("Similarity in community composition") +
  theme_light()

all <- bind_rows(
  simka_bc %>%
    left_join(mitag_bc, by = c("item1", "item2")) %>% mutate(class = "simka vs mitag"),
  simka_bc %>%
    left_join(eggnog_bc, by = c("item1", "item2")) %>% mutate(class = "simka vs eggnog"),
  simka_bc %>%
    left_join(unk_reads_bc, by = c("item1", "item2")) %>% mutate(class = "simka vs unks"),
  simka_bc %>%
    left_join(asv_bc, by = c("item1", "item2")) %>% mutate(class = "simka vs ASV"),
  unk_reads_bc %>%
    left_join(asv_bc, by = c("item1", "item2")) %>% mutate(class = "unks vs ASV"),
  unk_reads_bc %>%
    left_join(mitag_bc, by = c("item1", "item2")) %>% mutate(class = "unks vs mitag"),
  unk_reads_bc %>%
    left_join(eggnog_bc, by = c("item1", "item2")) %>% mutate(class = "unks vs eggnog"))

all <- all %>%
  left_join(osd2014_cdata %>% select(label, contains("meow")) %>% left_join(osd2014_meow_regions) %>% dplyr::rename(item1 = label)) %>%
  left_join(osd2014_cdata %>% select(label, contains("meow")) %>% left_join(osd2014_meow_regions) %>% dplyr::rename(item2 = label), by = "item2") %>%
  mutate(same_meow = ifelse(meow_province.x == meow_province.y, TRUE, FALSE)) %>% as_tibble()


all <- all %>%
  mutate(region.x = ifelse(is.na(meow_province.x), "Other", meow_province.x),
         region.y = ifelse(is.na(meow_province.y), "Other", meow_province.y),
         same_region = ifelse(meow_province.x == meow_province.y, meow_province.y, "Other")) %>%
  filter(same_region %in% osd2014_meow_regions$meow_province)



all_filt_s <- all %>% filter(same_region != "Other", grepl("simka", class)) %>% select(distance.x, distance.y, class)
all_filt_u <- all %>% filter(same_region != "Other", !grepl("simka", class)) %>% select(distance.x, distance.y, class)



p1 <- ggplot(all %>% filter(same_region != "Other", grepl("simka", class)) %>% mutate(meow_region.y = fct_relevel(meow_region.y, osd2014_meow_regions$meow_region %>% unique %>% sort(decreasing = T))), aes(1 - distance.x, 1 - distance.y)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(data = all_filt_s, colour = "black", alpha = .2, shape = 21, fill = "grey") +
  #geom_smooth(aes(color = class), method = "gam", formula = y ~ s(x, bs = "cs"), size = 0.5) +
  geom_point(shape = 21, color = "black", aes(fill = class)) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  ylab("Bray-curtis similarity") +
  xlab("Bray-curtis similarity") +
  theme_bw() +
  facet_wrap(.~meow_region.y, scales = "free")

p2 <- ggplot(all %>% filter(same_region != "Other", !grepl("simka", class)) %>%
         mutate(meow_region.y = fct_relevel(meow_region.y, osd2014_meow_regions$meow_region %>% unique %>% sort(decreasing = T))),
       aes(1 - distance.x, 1 - distance.y)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(data = all_filt_u, colour = "black", alpha = .2, shape = 21, fill = "grey") +
  #geom_smooth(aes(color = class), method = "gam", formula = y ~ s(x, bs = "cs"), size = 0.5) +
  geom_point(shape = 21, color = "black", aes(fill = class)) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  ylab("Bray-curtis similarity") +
  xlab("Bray-curtis similarity") +
  theme_bw() +
  facet_wrap(.~meow_region.y, scales = "free")

ggpubr::ggarrange(p1, p2, nrow = 2)



