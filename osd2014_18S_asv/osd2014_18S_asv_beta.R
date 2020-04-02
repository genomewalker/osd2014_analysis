library(phyloseq)
library(vegan)
library(tidyverse)
library(RPostgreSQL)
library(adespatial)
library(microbiomeSeq)
library(GUniFrac)
library(ggpubr)
library(DESeq2)

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

# load("osd2014_18S_asv/data/osd2014_18S_asv_beta.Rdata", verbose = TRUE)
# load(url("http://osd2014.metagenomics.eu/osd2014_18S_asv/data/osd2014_18S_asv_beta.Rdata"), verbose = TRUE)

# END: WARNING!! ---------------------------------------------------------------




# BEGIN: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------

# Load necessary data -----------------------------------------------------
# Use if you have the postgres DB in place
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)
osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)
st_100_order_terrestrial <- tbl(my_db, "osd2014_st_order_coastal") %>%
  collect(n = Inf)
osd2014_cdata <- osd2014_cdata %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% left_join(osd2014_meow_regions)
osd2014_haversine_distance <- tbl(my_db, "osd2014_haversine_distance") %>%
  collect(n = Inf) %>%
  dplyr::rename(haversine = distance) %>%
  inner_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  inner_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label")) %>%
  mutate(same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE))
osd2014_selected_meow_provinces <-tbl(my_db, "osd2014_selected_meow_provinces") %>%
  collect(n = Inf)
osd2014_simka_k31_bc <- tbl(my_db, "osd2014_simka_k31_bc") %>%
  collect(n = Inf)
osd2014_simka_k21_bc <- tbl(my_db, "osd2014_simka_k21_bc") %>%
  collect(n = Inf)


# If downloaded file at osd2014_18S_asv/data/ use:
load("osd2014_18S_asv/data/osd2014_18S_asv_physeq_filt_objects_with_phylo.Rdata", verbose = TRUE)
load("osd2014_18S_asv/data/osd2014_18S_asv_pina_tina_results.Rdata", verbose = TRUE)
load("osd2014_18S_asv/data/osd2014_18S_bMNTD_bNTI_results.Rdata", verbose = TRUE)
load("osd2014_18S_asv/data/osd2014_rc-analysis_9999.Rda", verbose = TRUE)

# Basic contextual data
load("osd2014_18S_asv/data/osd2014_basic_cdata.Rdata", verbose = TRUE)

# If remote use
load(url("http://osd2014.metagenomics.eu/osd2014_18S_asv/data/osd2014_18S_asv_physeq_filt_objects_with_phylo.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_18S_asv/data/osd2014_18S_asv_pina_tina_results.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_18S_asv/data/osd2014_mitag_qiime97_phyloseq.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_18S_asv/data/osd2014_bMNTD_bNTI_results.Rda"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_18S_asv/data/osd2014_rc-analysis_9999.Rda"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_18S_asv/data/osd2014_simka_bc.Rdata"), verbose = TRUE)

# Basic contextual data
load(url("http://osd2014.metagenomics.eu/osd2014_18S_asv/data/osd2014_basic_cdata.Rdata"), verbose = TRUE)
# Load necessary data -----------------------------------------------------

# END: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------







# Rarefy samples ----------------------------------------------------------

osd2014_dada2_phyloseq_beta_rar <- rarefy_even_depth(osd2014_dada2_phyloseq_beta)
osd2014_dada2_phyloseq_beta_rar_prop <- transform_sample_counts(osd2014_dada2_phyloseq_beta_rar, function(X)X/sum(X))

osd2014_dada2_phyloseq_beta_rar_deseq <- phyloseq_to_deseq2(osd2014_dada2_phyloseq_beta_rar, ~1)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(osd2014_dada2_phyloseq_beta_rar_deseq), 1, gm_mean)
diagdds <- estimateSizeFactors(osd2014_dada2_phyloseq_beta_rar_deseq, type = "poscounts")
diagdds <- estimateDispersions(diagdds)
diagvst <- varianceStabilizingTransformation(diagdds, blind = FALSE)
diagvst <- assay(diagvst)
diagvst[diagvst < 0] <- 0
osd2014_dada2_phyloseq_beta_rar_vst <- osd2014_dada2_phyloseq_beta_rar
osd2014_dada2_phyloseq_beta_rar_norm <- osd2014_dada2_phyloseq_beta_rar
otu_table(osd2014_dada2_phyloseq_beta_rar_vst) <- otu_table(diagvst, taxa_are_rows = TRUE)
otu_table(osd2014_dada2_phyloseq_beta_rar_norm) <- otu_table(counts(diagdds, normalized = TRUE), taxa_are_rows = TRUE)


# Explore relationships between distances ---------------------------------

osd2014_18S_asv_tina_uw <- osd2014_pina_tina_results[[1]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[1]]$name) %>%
  as_tibble()
osd2014_18S_asv_tina_w <- osd2014_pina_tina_results[[2]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[2]]$name) %>%
  as_tibble()
osd2014_18S_asv_tina_uw_filt <- osd2014_pina_tina_results[[3]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[3]]$name) %>%
  as_tibble()
osd2014_18S_asv_tina_w_filt <- osd2014_pina_tina_results[[4]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[4]]$name) %>%
  as_tibble()
osd2014_18S_asv_pina_uw <- osd2014_pina_tina_results[[5]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[5]]$name) %>%
  as_tibble()
osd2014_18S_asv_pina_w <- osd2014_pina_tina_results[[6]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[6]]$name) %>%
  as_tibble()
osd2014_18S_asv_bc <- phyloseq::distance(osd2014_dada2_phyloseq_beta_vst, "bray") %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = "Bray-Curtis") %>%
  as_tibble()

osd2014_18S_asv_bc_rar <- phyloseq::distance(osd2014_dada2_phyloseq_beta_rar_prop, "bray") %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = "Bray-Curtis rar") %>%
  as_tibble()

osd2014_18S_asv_jc <- phyloseq::distance(osd2014_dada2_phyloseq_beta_vst, method = "jaccard", binary = FALSE) %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = "Jaccard") %>%
  as_tibble()

# osd2014_mitag_qiime97_phyloseq_t <- phyloseq:::veganifyOTU(osd2014_mitag_qiime97_phyloseq)
# osd2014_mitag_qiime97_phyloseq_t <- osd2014_mitag_qiime97_phyloseq_t[rownames(osd2014_pina_tina_results[[4]]$cs), ]
# osd2014_mitag_bc <- vegan::vegdist(osd2014_mitag_qiime97_phyloseq_t) %>%
#   as.dist() %>%
#   broom::tidy() %>%
#   mutate(class = "Bray-Curtis") %>%
#   as_tibble()


osd2014_18S_asv_wu <- phyloseq::distance(osd2014_dada2_phyloseq_beta_vst, method = "wunifrac") %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = "Weighted Unifrac") %>%
  as_tibble()

osd2014_18S_asv_wu_rar <- phyloseq::distance(osd2014_dada2_phyloseq_beta_rar_prop, method = "wunifrac") %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = "Weighted Unifrac rar") %>%
  as_tibble()


osd2014_18S_asv_gu <- GUniFrac(otu.tab = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst), phy_tree(osd2014_dada2_phyloseq_beta_vst))$unifracs[, , "d_0.5"]%>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = "Generalised Unifrac") %>%
  as_tibble()

osd2014_18S_asv_gu_rar <- GUniFrac(otu.tab = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_rar_prop), phy_tree(osd2014_dada2_phyloseq_beta_rar))$unifracs[, , "d_0.5"]%>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = "Generalised Unifrac rar") %>%
  as_tibble()


osd2014_18S_asv_tina_uw %>%
  left_join(osd2014_18S_asv_tina_w, by = c("item1", "item2")) %>% mutate(class = "TINA, unweighted vs TINA, weighted") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("TINA, unweighted") +
  ylab("TINA, weighted")

osd2014_18S_asv_pina_uw %>%
  left_join(osd2014_18S_asv_pina_w, by = c("item1", "item2")) %>% mutate(class = "PINA, unweighted vs PINA, weighted") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("PINA, unweighted") +
  ylab("PINA, weighted") +
  theme_bw() +
  stat_cor()

osd2014_18S_asv_tina_w %>%
  left_join(osd2014_18S_asv_pina_w, by = c("item1", "item2")) %>% mutate(class = "TINA, weighted vs PINA, weighted") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("TINA, weighted") +
  ylab("PINA, weighted") +
  theme_bw() +
  stat_cor()

osd2014_18S_asv_tina_w %>%
  left_join(osd2014_18S_asv_bc, by = c("item1", "item2")) %>% mutate(class = "TINA, weighted vs ASV BC") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("TINA, weighted") +
  ylab("ASV BC")  +
  theme_bw() +
  stat_cor()

# osd2014_18S_asv_tina_w %>%
#   left_join(osd2014_mitag_bc, by = c("item1", "item2")) %>% mutate(class = "TINA, weighted vs miTAG BC") %>%
#   ggplot(aes(distance.x, distance.y)) +
#   geom_point() +
#   geom_smooth() +
#   geom_abline() +
#   xlim(c(0,1)) +
#   ylim(c(0,1)) +
#   xlab("TINA, weighted") +
#   ylab("miTAG BC")

osd2014_18S_asv_tina_w %>%
  left_join(osd2014_simka_k21_bc, by = c("item1", "item2")) %>% mutate(class = "TINA, weighted vs SIMKA BC k=21") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("TINA, weighted") +
  ylab("SIMKA BC k=21")

osd2014_18S_asv_tina_w %>%
  left_join(osd2014_simka_k31_bc, by = c("item1", "item2")) %>% mutate(class = "TINA, weighted vs SIMKA BC k=31") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("TINA, weighted") +
  ylab("SIMKA BC k=31")

osd2014_18S_asv_bc %>%
  left_join(osd2014_simka_k21_bc, by = c("item1", "item2")) %>% mutate(class = "ASV BC vs SIMKA BC k=21") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("ASV BC") +
  ylab("SIMKA BC k=21")

osd2014_18S_asv_bc %>%
  left_join(osd2014_simka_k31_bc, by = c("item1", "item2")) %>% mutate(class = "ASV BC vs SIMKA BC k=31") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("ASV BC") +
  ylab("SIMKA BC k=31")

osd2014_18S_asv_bc %>%
  left_join(osd2014_18S_asv_bc_rar, by = c("item1", "item2")) %>% mutate(class = "ASV BC vs ASV BC rar") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("ASV BC") +
  ylab("ASV BC rar") +
  theme_bw() +
  stat_cor()

# osd2014_18S_asv_bc %>%
#   left_join(osd2014_mitag_bc, by = c("item1", "item2")) %>% mutate(class = "ASV BC vs miTAG BC") %>%
#   ggplot(aes(distance.x, distance.y)) +
#   geom_point() +
#   geom_smooth() +
#   geom_abline() +
#   xlim(c(0,1)) +
#   ylim(c(0,1)) +
#   xlab("ASV BC") +
#   ylab("miTAG BC") +
#   theme_bw() +
#   stat_cor()

osd2014_18S_asv_wu %>%
  left_join(osd2014_18S_asv_wu_rar, by = c("item1", "item2")) %>% mutate(class = "ASV WU vs ASV WU rar") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlab("ASV WU") +
  ylab("ASV WU rar") +
  theme_bw() +
  stat_cor()

osd2014_18S_asv_gu %>%
  left_join(osd2014_18S_asv_gu_rar, by = c("item1", "item2")) %>% mutate(class = "ASV GU vs ASV GU rar") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlab("ASV GU") +
  ylab("ASV GU rar") +
  theme_bw() +
  stat_cor()


osd2014_18S_asv_tina_w_meow <- osd2014_18S_asv_tina_w_filt %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label"))

osd2014_18S_asv_pina_w_meow <- osd2014_18S_asv_pina_w %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label"))

osd2014_18S_asv_bc_meow <- osd2014_18S_asv_bc%>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label"))

osd2014_simka_k21_bc_meow <- osd2014_simka_k21_bc%>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label"))

osd2014_18S_asv_gu_meow <- osd2014_18S_asv_gu %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label"))

bind_rows(osd2014_18S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            mutate(class = "TINA weighted"),
          osd2014_18S_asv_pina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            mutate(class = "PINA unweighted"),
          osd2014_18S_asv_bc_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            mutate(class = "ASV BC"),
          osd2014_18S_asv_gu_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            mutate(class = "ASV GU")
          # osd2014_simka_k21_bc_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
          #   filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
          #   filter(same_region ==  TRUE) %>% mutate(class = "SIMKA BC k = 21")
) %>%
  ggplot(aes(x = class, y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5) +
  facet_wrap(~meow_region.x, scales = "free")


bind_rows(osd2014_18S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            filter(same_region ==  TRUE) %>% mutate(class = "TINA weighted"),
          osd2014_18S_asv_pina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            filter(same_region ==  TRUE) %>% mutate(class = "PINA unweighted"),
          osd2014_18S_asv_bc_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            filter(same_region ==  TRUE) %>% mutate(class = "ASV BC"),
          osd2014_18S_asv_gu_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            filter(same_region ==  TRUE) %>% mutate(class = "ASV GU")
          # osd2014_simka_k21_bc_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
          #   filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
          #   filter(same_region ==  TRUE) %>% mutate(class = "SIMKA BC k = 21")
) %>%
  ggplot(aes(x = meow_province.x, y = distance, color = class, fill = class)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5) +
  facet_wrap(~meow_region.x, scales = "free")



# We will select wTINA, gUNIFRAC and BC -----------------------------------
osd2014_cdata_df <- osd2014_cdata %>% filter(label %in% sample_names(osd2014_dada2_phyloseq_beta), !is.na(meow_region)) %>% as.data.frame() %>% column_to_rownames("label")

selected_samples <- osd2014_cdata %>% filter(label %in% sample_names(osd2014_dada2_phyloseq_beta), !is.na(meow_region)) %>% .$label
osd2014_dada2_phyloseq_beta_vst_sel <- prune_samples(selected_samples, osd2014_dada2_phyloseq_beta_vst)
osd2014_dada2_phyloseq_beta_vst_sel <-  prune_taxa(taxa_sums(osd2014_dada2_phyloseq_beta_vst_sel) > 0, osd2014_dada2_phyloseq_beta_vst_sel)


# 1: "TINA, unweighted", 2: "TINA, weighted", 3: "TINA, unweighted filtered", 4: "TINA, weighted filtered", 5: "PINA, unweighted", 6: "PINA, weighted"
m_tina <- metaMDS(as.dist(osd2014_pina_tina_results[[2]]$cs[selected_samples,selected_samples]),  autotransform = TRUE,  trymax = 1000)
m_tina_filt <- metaMDS(as.dist(osd2014_pina_tina_results[[4]]$cs[selected_samples,selected_samples]),  autotransform = TRUE, trymax = 1000)
m_pina <- metaMDS(as.dist(osd2014_pina_tina_results[[6]]$cs[selected_samples,selected_samples]),  autotransform = TRUE, trymax = 1000)
m_bc <- metaMDS(vegdist(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), method = "bray"),  autotransform = TRUE, trymax = 1000)
m_uf <- metaMDS(GUniFrac(otu.tab = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), phy_tree(osd2014_dada2_phyloseq_beta_vst_sel))$unifracs[, , "d_0.5"],  autotransform = TRUE, trymax = 1000, noshare = 0.01)



p_tina <- scores(m_tina) %>%
  as_tibble(rownames = "label") %>%
  left_join(osd2014_cdata) %>%
  mutate(is_region = ifelse(meow_province %in% osd2014_meow_regions$meow_province, meow_province, "Other")) %>%
  filter(is_region != "Other") %>%
  mutate(is_region = fct_relevel(is_region, c("Tropical Northwestern Atlantic", "Warm Temperate Northwest Atlantic",
                                              "Cold Temperate Northwest Atlantic", "Lusitanian", "Mediterranean Sea", "Northern European Seas")))

colors <- (c("#C059D1", "#D2D0B7", "#B99ED0", "#ABE172", "#7AD4CB",  "#DB8763"))


surf <- function(X, Y) {
  ordi <- vegan::ordisurf(X, Y, plot = FALSE, bs="ds")
  ordi.grid <- ordi$grid #extracts the ordisurf object
  #str(ordi.grid) #it's a list though - cannot be plotted as is
  ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
  ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
  ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
}

p1 <- ggplot() +
  ggConvexHull::geom_convexhull(data = p_tina, aes(NMDS1, NMDS2, color = is_region), fill = "grey80", alpha = 0.2, size = 0.5) +
  #stat_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, color = ..level..)) +
  geom_point(data = p_tina, aes(NMDS1, NMDS2, fill = is_region), color = "grey20", shape = 21, size = 1) +
  coord_fixed() +
  theme_bw() +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) + theme(legend.position = "none")
p1

p2 <- ggplot() +
  geom_rug(data = p_tina, aes(NMDS1, NMDS2, color = water_temperature), sides = "b") +
  coord_fixed() +
  theme_bw() +
  scale_color_gradient2(low = "blue", mid = "blue", high = "red", midpoint = 15) +
  theme(legend.position = "right")
p2

ggarrange(p1, p2, legend = "top")
ggsave(plot = last_plot(), filename = "osd2014_18S_asv/figures/osd2014_dada2_phyloseq_beta_nmds.pdf", width = 6, height = 4)


# p1 <- ggplot() +
#   #stat_contour(data = ordi.mite.na, aes(x = x, y = y, z = z, color = ..level..)) +
#   ggConvexHull::geom_convexhull(data = p_tina %>% rename(is_region1 = is_region), aes(NMDS1, NMDS2, group = is_region1), alpha = 0.2, size = 0.1, fill = NA, color = "black", linetype = "dashed") +
#   geom_point(data = p_tina %>% select(-is_region), aes(NMDS1, NMDS2), color = "grey20", shape = 21, fill = "grey80") +
#   ggConvexHull::geom_convexhull(data = p_tina, aes(NMDS1, NMDS2), fill = "grey80", alpha = 0.2, size = 1, color = "black") +
#   geom_point(data = p_tina, aes(NMDS1, NMDS2, fill = water_temperature), color = "grey20", shape = 21) +
#   coord_fixed() +
#   theme_light() +
#   scale_fill_continuous() +
#   scale_color_manual(values = colors) + theme(legend.position = "none") +
#   facet_wrap(~is_region)
# p1

ggarrange(p1, p2, legend = "none")

p_tina <- scores(m_tina) %>%
  as_tibble(rownames = "label") %>%
  left_join(osd2014_cdata) %>%
  mutate(is_region = ifelse(meow_province %in% osd2014_meow_regions$meow_province, meow_province, "Other")) %>%
  filter(is_region != "Other") %>%
  ggplot(aes(NMDS1, NMDS2, color = water_temperature, fill = is_region, shape = is_region)) +
  ggConvexHull::geom_convexhull(alpha = 0.3, color = "grey50") +
  geom_point() +
  coord_fixed() +
  theme_light() +
  viridis::scale_color_viridis(option = "B") +
  ggtitle("TINA")

p_tina_filt <- scores(m_tina_filt) %>%
  as_tibble(rownames = "label") %>%
  left_join(osd2014_cdata) %>%
  mutate(is_region = ifelse(meow_province %in% osd2014_meow_regions$meow_province, meow_province, "Other")) %>%
  filter(is_region != "Other") %>%
  ggplot(aes(NMDS1, NMDS2, color = water_temperature, fill = is_region, shape = is_region)) +
  ggConvexHull::geom_convexhull(alpha = 0.3, color = "grey50") +
  geom_point() +
  coord_fixed() +
  theme_light() +
  viridis::scale_color_viridis(option = "B") +
  ggtitle("TINA filt")

p_pina <- scores(m_pina) %>%
  as_tibble(rownames = "label") %>%
  left_join(osd2014_cdata) %>%
  mutate(is_region = ifelse(meow_province %in% osd2014_meow_regions$meow_province, meow_province, "Other")) %>%
  filter(is_region != "Other") %>%
  ggplot(aes(NMDS1, NMDS2, color = water_temperature, fill = is_region, shape = is_region)) +
  ggConvexHull::geom_convexhull(alpha = 0.3, color = "grey50") +
  geom_point() +
  coord_fixed() +
  theme_light() +
  viridis::scale_color_viridis(option = "B") +
  ggtitle("PINA")

p_bc <- scores(m_bc) %>% as_tibble(rownames = "label") %>%
  left_join(osd2014_cdata) %>%
  mutate(is_region = ifelse(meow_province %in% osd2014_meow_regions$meow_province, meow_province, "Other")) %>%
  filter(is_region != "Other") %>%
  ggplot(aes(NMDS1, NMDS2, color = water_temperature, fill = is_region, shape = is_region)) +
  ggConvexHull::geom_convexhull(alpha = 0.3, color = "grey50") +
  geom_point() +
  coord_fixed() +
  theme_light() +
  viridis::scale_color_viridis(option = "B") +
  ggtitle("Bray-Curtis")


p_uf <- scores(m_uf) %>% as_tibble(rownames = "label") %>%
  left_join(osd2014_cdata) %>%
  mutate(is_region = ifelse(meow_province %in% osd2014_meow_regions$meow_province, meow_province, "Other")) %>%
  filter(is_region != "Other") %>%
  ggplot(aes(NMDS1, NMDS2, color = water_temperature, fill = is_region, shape = is_region)) +
  ggConvexHull::geom_convexhull(alpha = 0.3, color = "grey50") +
  geom_point() +
  coord_fixed() +
  theme_light() +
  viridis::scale_color_viridis(option = "B") +
  ggtitle("UNIFRAC")

ggpubr::ggarrange(p_tina, p_tina_filt, p_bc, p_uf, common.legend = TRUE, align = "hv", nrow = 1, ncol = 4)


perm_tina_mp <- adonis(as.dist(osd2014_pina_tina_results[[2]]$cs[selected_samples,selected_samples]) ~ meow_province, data = osd2014_cdata_df)
perm_tina_wt <- adonis(as.dist(osd2014_pina_tina_results[[2]]$cs[selected_samples,selected_samples]) ~ water_temperature, data = osd2014_cdata_df)
perm_tina_iso <- adonis(as.dist(osd2014_pina_tina_results[[2]]$cs[selected_samples,selected_samples]) ~ dist_coast_iso3_code, data = osd2014_cdata_df)
perm_tina_mp_wt <- adonis(as.dist(osd2014_pina_tina_results[[2]]$cs[selected_samples,selected_samples]) ~ meow_province + water_temperature, data = osd2014_cdata_df)
betadisper_tina_mp <- betadisper(as.dist(osd2014_pina_tina_results[[2]]$cs[selected_samples,selected_samples]), osd2014_cdata_df$meow_province, type = "centroid")
p_betadisper_tina_mp <- permutest(betadisper_tina_mp)
betadisper_tina_iso <- betadisper(as.dist(osd2014_pina_tina_results[[2]]$cs[selected_samples,selected_samples]), osd2014_cdata_df$dist_coast_iso3_code, type = "centroid")
p_betadisper_tina_iso <- permutest(betadisper_tina_iso)


perm_tina_filt_mp <- adonis(as.dist(osd2014_pina_tina_results[[4]]$cs[selected_samples,selected_samples]) ~ meow_province, data = osd2014_cdata_df)
perm_tina_filt_wt <- adonis(as.dist(osd2014_pina_tina_results[[4]]$cs[selected_samples,selected_samples]) ~ water_temperature, data = osd2014_cdata_df)
perm_tina_filt_iso <- adonis(as.dist(osd2014_pina_tina_results[[4]]$cs[selected_samples,selected_samples]) ~ dist_coast_iso3_code, data = osd2014_cdata_df)
perm_tina_filt_mp_wt <- adonis(as.dist(osd2014_pina_tina_results[[4]]$cs[selected_samples,selected_samples]) ~ meow_province + water_temperature, data = osd2014_cdata_df)
betadisper_tina_filt_mp <- betadisper(as.dist(osd2014_pina_tina_results[[4]]$cs[selected_samples,selected_samples]), osd2014_cdata_df$meow_province, type = "centroid")
p_betadisper_tina_filt_mp <- permutest(betadisper_tina_filt_mp)
betadisper_tina_filt_iso <- betadisper(as.dist(osd2014_pina_tina_results[[4]]$cs[selected_samples,selected_samples]), osd2014_cdata_df$dist_coast_iso3_code, type = "centroid")
p_betadisper_tina_filt_iso <- permutest(betadisper_tina_filt_iso)



perm_pina_mp <- adonis(as.dist(osd2014_pina_tina_results[[6]]$cs[selected_samples,selected_samples]) ~ meow_province, data = osd2014_cdata_df)
perm_pina_wt <- adonis(as.dist(osd2014_pina_tina_results[[6]]$cs[selected_samples,selected_samples]) ~ water_temperature, data = osd2014_cdata_df)
perm_pina_iso <- adonis(as.dist(osd2014_pina_tina_results[[6]]$cs[selected_samples,selected_samples]) ~ dist_coast_iso3_code, data = osd2014_cdata_df)
perm_pina_mp_wt <- adonis(as.dist(osd2014_pina_tina_results[[6]]$cs[selected_samples,selected_samples]) ~ meow_province + water_temperature, data = osd2014_cdata_df)
betadisper_pina_mp <- betadisper(as.dist(osd2014_pina_tina_results[[6]]$cs[selected_samples,selected_samples]), osd2014_cdata_df$meow_province, type = "centroid")
p_betadisper_pina_mp <- permutest(betadisper_pina_mp)
betadisper_pina_iso <- betadisper(as.dist(osd2014_pina_tina_results[[6]]$cs[selected_samples,selected_samples]), osd2014_cdata_df$dist_coast_iso3_code, type = "centroid")
p_betadisper_pina_iso <- permutest(betadisper_pina_iso)


perm_bc_mp <- adonis(vegdist(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), method = "bray") ~ meow_province, data = osd2014_cdata_df)
perm_bc_wt <- adonis(vegdist(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), method = "bray") ~ water_temperature, data = osd2014_cdata_df)
perm_bc_iso <- adonis(vegdist(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), method = "bray") ~ dist_coast_iso3_code, data = osd2014_cdata_df)
perm_bc_mp_wt <- adonis(vegdist(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), method = "bray") ~ meow_province + water_temperature, data = osd2014_cdata_df)
betadisper_bc_mp <- betadisper(vegdist(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), method = "bray"), osd2014_cdata_df$meow_province, type = "centroid")
p_betadisper_bc_mp <- permutest(betadisper_bc_mp)
betadisper_bc_iso <- betadisper(vegdist(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), method = "bray"), osd2014_cdata_df$dist_coast_iso3_code, type = "centroid")
p_betadisper_bc_iso <- permutest(betadisper_bc_iso)


perm_uf_mp <- adonis(as.dist(GUniFrac(otu.tab = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), phy_tree(osd2014_dada2_phyloseq_beta_vst_sel))$unifracs[, , "d_0.5"]) ~ meow_province, data = osd2014_cdata_df)
perm_uf_wt <- adonis(as.dist(GUniFrac(otu.tab = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), phy_tree(osd2014_dada2_phyloseq_beta_vst_sel))$unifracs[, , "d_0.5"]) ~ water_temperature, data = osd2014_cdata_df)
perm_uf_iso <- adonis(as.dist(GUniFrac(otu.tab = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), phy_tree(osd2014_dada2_phyloseq_beta_vst_sel))$unifracs[, , "d_0.5"]) ~ dist_coast_iso3_code, data = osd2014_cdata_df)
perm_uf_mp_wt <- adonis(as.dist(GUniFrac(otu.tab = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), phy_tree(osd2014_dada2_phyloseq_beta_vst_sel))$unifracs[, , "d_0.5"]) ~ meow_province + water_temperature, data = osd2014_cdata_df)
betadisper_uf_mp <- betadisper(as.dist(GUniFrac(otu.tab = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), phy_tree(osd2014_dada2_phyloseq_beta_vst_sel))$unifracs[, , "d_0.5"]), osd2014_cdata_df$meow_province, type = "centroid")
p_betadisper_uf_mp <- permutest(betadisper_uf_mp)
betadisper_uf_iso <- betadisper(as.dist(GUniFrac(otu.tab = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst_sel), phy_tree(osd2014_dada2_phyloseq_beta_vst_sel))$unifracs[, , "d_0.5"]), osd2014_cdata_df$dist_coast_iso3_code, type = "centroid")
p_betadisper_uf_iso <- permutest(betadisper_uf_iso)


osd2014_permanova_results <- bind_rows(
  tibble(class = "wTINA", group = "MEOW province", r2 = perm_tina_mp$aov.tab$R2[[1]], F_Model = perm_tina_mp$aov.tab$F.Model[[1]], pvalue = perm_tina_mp$aov.tab$`Pr(>F)`[[1]]),
  tibble(class = "wTINA", group = "Temperature", r2 = perm_tina_wt$aov.tab$R2[[1]], F_Model = perm_tina_wt$aov.tab$F.Model[[1]], pvalue = perm_tina_wt$aov.tab$`Pr(>F)`[[1]]),
  #tibble(class = "wTINA", group = "Closest country", r2 = perm_tina_iso$aov.tab$R2[[1]], F_Model = perm_tina_iso$aov.tab$F.Model[[1]], pvalue = perm_tina_iso$aov.tab$`Pr(>F)`[[1]]),
  # tibble(class = "wTINA filtered", group = "MEOW province", r2 = perm_tina_filt_mp$aov.tab$R2[[1]], F_Model = perm_tina_filt_mp$aov.tab$F.Model[[1]], pvalue = perm_tina_filt_mp$aov.tab$`Pr(>F)`[[1]]),
  tibble(class = "Bray-Curtis", group = "MEOW province", r2 = perm_bc_mp$aov.tab$R2[[1]], F_Model = perm_bc_mp$aov.tab$F.Model[[1]], pvalue = perm_bc_mp$aov.tab$`Pr(>F)`[[1]]),
  tibble(class = "Bray-Curtis", group = "Temperature", r2 = perm_bc_iso$aov.tab$R2[[1]], F_Model = perm_bc_iso$aov.tab$F.Model[[1]], pvalue = perm_bc_iso$aov.tab$`Pr(>F)`[[1]]),
  # tibble(class = "wTINA filtered", group = "Temperature", r2 = perm_tina_filt_wt$aov.tab$R2[[1]], F_Model = perm_tina_filt_wt$aov.tab$F.Model[[1]], pvalue = perm_tina_filt_wt$aov.tab$`Pr(>F)`[[1]]),
  # #tibble(class = "wTINA filtered", group = "Closest country", r2 = perm_tina_filt_iso$aov.tab$R2[[1]], F_Model = perm_tina_filt_iso$aov.tab$F.Model[[1]], pvalue = perm_tina_filt_iso$aov.tab$`Pr(>F)`[[1]]),
  tibble(class = "wPINA", group = "MEOW province", r2 = perm_pina_mp$aov.tab$R2[[1]], F_Model = perm_pina_mp$aov.tab$F.Model[[1]], pvalue = perm_pina_mp$aov.tab$`Pr(>F)`[[1]]),
  tibble(class = "wPINA", group = "Temperature", r2 = perm_pina_wt$aov.tab$R2[[1]], F_Model = perm_pina_wt$aov.tab$F.Model[[1]], pvalue = perm_pina_wt$aov.tab$`Pr(>F)`[[1]]),
  #tibble(class = "wPINA", group = "Closest country", r2 = perm_pina_iso$aov.tab$R2[[1]], F_Model = perm_pina_iso$aov.tab$F.Model[[1]], pvalue = perm_pina_iso$aov.tab$`Pr(>F)`[[1]]),
  #tibble(class = "Bray-Curtis", group = "Closest country", r2 = perm_bc_wt$aov.tab$R2[[1]], F_Model = perm_bc_wt$aov.tab$F.Model[[1]], pvalue = perm_bc_wt$aov.tab$`Pr(>F)`[[1]]),
  tibble(class = "gUnifrac", group = "MEOW province", r2 = perm_uf_mp$aov.tab$R2[[1]], F_Model = perm_uf_mp$aov.tab$F.Model[[1]], pvalue = perm_uf_mp$aov.tab$`Pr(>F)`[[1]]),
  tibble(class = "gUnifrac", group = "Temperature", r2 = perm_uf_wt$aov.tab$R2[[1]], F_Model = perm_uf_wt$aov.tab$F.Model[[1]], pvalue = perm_uf_wt$aov.tab$`Pr(>F)`[[1]]),
  #tibble(class = "gUnifrac", group = "Closest country", r2 = perm_uf_iso$aov.tab$R2[[1]], F_Model = perm_uf_iso$aov.tab$F.Model[[1]], pvalue = perm_uf_iso$aov.tab$`Pr(>F)`[[1]]),
) %>%
  mutate(signif = ifelse(pvalue <= 0.001, "Significant", "No significant"),
         class = fct_relevel(class, rev(c("wTINA", "Bray-Curtis", "wPINA", "gUnifrac"))))


ggplot(osd2014_permanova_results, aes(x = class, y = r2, fill = signif)) +
  geom_point(shape = 22, size = 2) +
  facet_wrap(~group) +
  ggpubr::rotate() +
  theme_bw() +
  scale_y_continuous() +
  scale_fill_manual(values = c("black", "red")) +
  theme_cleveland() +
  ylab("Permanova R2") +
  theme(legend.position = "top")

ggsave(plot = last_plot(), filename = "osd2014_18S_asv/figures/osd2014_dada2_phyloseq_beta_permanova.pdf", width = 4, height = 3)



# Local Contribution to Beta Diversity (LCBD) -----------------------------


# bNTI & bMNTD analysis -----------------------------------------------------------

# bMNTD values higher than expected by chance indicate that
# communities are under heterogeneous selection (18). In contrast, bMNTD values which
# are lower than expected by chance indicate that communities are experiencing
# homogeneous selection.
# Selection can act in two opposite directions, it can constrain (homogeneous selection)
# or promote (heterogeneous selection) the divergence of communities



# LCBD --------------------------------------------------------------------
# Local contributions to beta diversity (LCBD indices) represent the degree of uniqueness of
# the sites in terms of their species compositions. They can be computed in all cases: raw
# (not recommended) or transformed data, as well as dissimilarity matrices. See Legendre and
# De Cáceres (2013) for details. LCBD indices are tested for significance by random, independent
# permutations within the columns of Y. This permutation method tests H0 that the species are
# distributed at random, independently of one another, among the sites, while preserving the
# species abundance distributions in the observed data. See Legendre and De Cáceres (2013) for discussion.

# RC bray values between -0.95 and +0.95 point to a community assembly governed by drift.
# On the contrary, RC bray values > +0.95 or < -0.95 indicate that community turnover is driven by dispersal
# limitation or homogenizing dispersal respectively

beta_div <- beta.div(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta), method = "hellinger", sqrt.D = F,
                     samp = T, nperm = 999)
df_LCBD <- tibble(label = names(beta_div$LCBD),
                  LCBD = beta_div$LCBD,
                  p.LCBD = beta_div$p.LCBD,
                  padj.LCBD = beta_div$p.adj)



# bNTI: distribution are denoted as b-Nearest Taxon Index (bNTI), with | bNTI | > 2 being
# considered as significant departures from random phylogenetic turnover, pointing to the
# action of selection.

selection_bNTI <- broom::tidy(as.dist(weighted.bNTI)) %>%
  as_tibble() %>%
  filter(!is.na(distance)) %>%
  rename(bNTI = distance) %>%
  inner_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  inner_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label")) %>%
  mutate(same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE), selection = ifelse(abs(bNTI) > 2, TRUE, FALSE)) %>%
  filter(item1 %in% selected_samples, item2 %in% selected_samples, same_region == TRUE, selection == TRUE)

no_selection_bNTI <- broom::tidy(as.dist(weighted.bNTI)) %>%
  as_tibble() %>%
  filter(!is.na(distance)) %>%
  rename(bNTI = distance) %>%
  inner_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  inner_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label")) %>%
  mutate(same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE), selection = ifelse(abs(bNTI) > 2, TRUE, FALSE)) %>%
  filter(item1 %in% selected_samples, item2 %in% selected_samples, same_region == TRUE, selection == FALSE) %>%
  inner_join(osd2014_18S_rc_results_9999) %>%
  mutate(rc_type = case_when(rc > 0.95 ~ "Dispersal limitation",
                             rc < -0.95 ~ "Homogenizing dispersal",
                             TRUE ~ "Undominated" ))


library(picante)
otu <- t(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst))
phylo <- phy_tree(osd2014_dada2_phyloseq_beta_vst)
match.phylo.otu <- match.phylo.data(phylo, otu)

bMNTD_rand_mean <- apply(rand.weighted.bMNTD.comp, c(1,2), mean)
rownames(bMNTD_rand_mean) <- colnames(match.phylo.otu$data)
colnames(bMNTD_rand_mean) <- colnames(match.phylo.otu$data)

bMNTD_rand_lowCI <- apply(rand.weighted.bMNTD.comp, c(1,2), quantile, probs = c(0.025))
rownames(bMNTD_rand_lowCI) <- colnames(match.phylo.otu$data)
colnames(bMNTD_rand_lowCI) <- colnames(match.phylo.otu$data)

bMNTD_rand_uppCI <- apply(rand.weighted.bMNTD.comp, c(1,2), quantile, probs = c(0.975))
rownames(bMNTD_rand_uppCI) <- colnames(match.phylo.otu$data)
colnames(bMNTD_rand_uppCI) <- colnames(match.phylo.otu$data)

bMNTD_rand_mean <- as.dist(bMNTD_rand_mean) %>%
  broom::tidy() %>%
  as_tibble() %>%
  filter(!is.na(distance)) %>%
  rename(bMNTD_rand_mean = distance)

bMNTD_rand_lowCI <- as.dist(bMNTD_rand_lowCI) %>%
  broom::tidy() %>%
  as_tibble() %>%
  filter(!is.na(distance)) %>%
  rename(bMNTD_rand_lowCI = distance)

bMNTD_rand_uppCI <- as.dist(bMNTD_rand_uppCI) %>%
  broom::tidy() %>%
  as_tibble() %>%
  filter(!is.na(distance)) %>%
  rename(bMNTD_rand_uppCI = distance)

bMNTD_rand <- bMNTD_rand_mean %>%
  inner_join(bMNTD_rand_lowCI) %>%
  inner_join(bMNTD_rand_uppCI)

bMNTD <- broom::tidy(as.dist(beta.mntd.weighted)) %>%
  as_tibble() %>%
  filter(!is.na(distance)) %>%
  rename(bMNTD = distance) %>%
  inner_join(bMNTD_rand, by = c("item1" = "item1", "item2" = "item2")) %>%
  mutate(type_selection = ifelse(bMNTD > bMNTD_rand_mean, "Heterogeneous", "Homogeneous")) %>%
  inner_join(selection_bNTI)

# bind_rows(bMNTD %>% select(item1, item2, selection, type_selection, meow_province.x) %>% mutate(rc_type = NA),
#           no_selection_bNTI %>% select(item1, item2, selection, rc_type,meow_province.x) %>% mutate(type_selection = NA)
#           ) %>%
#   group_by(meow_province.x, selection, type_selection) %>%
#   dplyr::count() %>%
#   group_by(meow_province.x) %>%
#   mutate(N = sum(n), prop = n/N) %>%
#   filter(selection == TRUE) %>%
#   ggplot(aes(meow_province.x, prop, fill = type_selection)) +
#   geom_col() +
#   ggpubr::rotate()

bind_rows(bMNTD %>% select(item1, item2, selection, type_selection, meow_province.x) %>% mutate(rc_type = NA),
          no_selection_bNTI %>% select(item1, item2, selection, rc_type, meow_province.x) %>% mutate(type_selection = NA)
) %>%
  group_by(meow_province.x, selection, rc_type, type_selection) %>%
  dplyr::count() %>%
  group_by(meow_province.x) %>%
  mutate(N = sum(n), prop = n/N) %>%
  ungroup() %>%
  mutate(meow_province.x = fct_relevel(meow_province.x, rev(c("Tropical Northwestern Atlantic", "Warm Temperate Northwest Atlantic", "Cold Temperate Northwest Atlantic",
                                                              "Lusitanian", "Mediterranean Sea", "Northern European Seas"))),
         rc_type = fct_relevel(rc_type, (c("Undominated", "Dispersal limitation", "Homogenizing dispersal")))) %>%
  #filter(selection == FALSE) %>%
  ggplot(aes(meow_province.x, prop, fill = rc_type, color = type_selection)) +
  geom_col() +
  ggpubr::rotate()+
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "top")

ggsave(plot = last_plot(), filename = "osd2014_18S_asv/figures/osd2014_dada2_phyloseq_beta_selection.pdf", width = 6, height = 8)

# Exploring bray-curtis/TINA divergence between transects -----------------
seq_bc <- function(X){
  l <- st_100_order_terrestrial %>%
    filter(label %in% (osd2014_cdata %>% filter(eval(parse(text = X))) %>% .$label))  %>%
    filter(label %in% selected_samples) %>%
    .$label %>% combn(, m = 2) %>% t %>%
    as_tibble() %>%
    inner_join(st_100_order_terrestrial, by = c("V1" = "label"))  %>%
    inner_join(st_100_order_terrestrial, by = c("V2" = "label")) %>%
    mutate(posdif = abs(position.x - position.y)) %>%
    filter(posdif == 1)
  s_bc <- bind_rows(lapply(1:nrow(l), function (x){
    physeq <- prune_samples(c(l$V1[[x]], l$V2[[x]]), osd2014_dada2_phyloseq_beta_vst)
    physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
    d <- phyloseq::distance(physeq = physeq, method = "bray")
    tibble(V1 = l$V1[[x]], V2 = l$V2[[x]], distance = d)
  })) %>%
    inner_join(l)
  physeq <- prune_samples(c(l$V1, l$V2) %>% unique(), osd2014_dada2_phyloseq_beta)
  physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  beta_div <- beta.div(phyloseq:::veganifyOTU(physeq), method = "hellinger", sqrt.D = F,
                       samp = T, nperm = 999)
  df_LCBD <- tibble(label = names(beta_div$LCBD),
                    LCBD = beta_div$LCBD,
                    p.LCBD = beta_div$p.LCBD,
                    padj.LCBD = beta_div$p.adj) %>%
    mutate(signif = ifelse(padj.LCBD < 0.05, "Significant", "No significant")) %>%
    left_join(st_100_order_terrestrial %>% select(label, position)) %>%
    left_join(s_bc %>% select(distance, position.x), by = c("position" = "position.x"))

  list(seq_bc = s_bc, LCB = df_LCBD, sample_order = l)
}


seq_tina <- function(X){
  l <- st_100_order_terrestrial %>%
    filter(label %in% (osd2014_cdata %>% filter(eval(parse(text = X))) %>% .$label))  %>%
    filter(label %in% selected_samples) %>%
    .$label %>% combn(, m = 2) %>% t %>%
    as_tibble() %>%
    inner_join(st_100_order_terrestrial, by = c("V1" = "label"))  %>%
    inner_join(st_100_order_terrestrial, by = c("V2" = "label")) %>%
    mutate(posdif = abs(position.x - position.y)) %>%
    filter(posdif == 1)
  s_tina <- bind_rows(lapply(1:nrow(l), function (x){
    d <- osd2014_pina_tina_results[[2]]$cs[l$V1[[x]], l$V2[[x]]]
    tibble(V1 = l$V1[[x]], V2 = l$V2[[x]], distance = d)
  })) %>%
    inner_join(l)

  physeq <- prune_samples(c(l$V1, l$V2) %>% unique(), osd2014_dada2_phyloseq_beta)
  physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  beta_div <- beta.div(phyloseq:::veganifyOTU(physeq), method = "hellinger", sqrt.D = F,
                       samp = T, nperm = 999)
  df_LCBD <- tibble(label = names(beta_div$LCBD),
                    LCBD = beta_div$LCBD,
                    p.LCBD = beta_div$p.LCBD,
                    padj.LCBD = beta_div$p.adj) %>%
    mutate(signif = ifelse(padj.LCBD < 0.05, "Significant", "No significant")) %>%
    left_join(st_100_order_terrestrial %>% select(label, position)) %>%
    left_join(s_tina %>% select(distance, position.x), by = c("position" = "position.x"))
  list(seq_tina = s_tina, LCB = df_LCBD, sample_order = l)
}


flt <- paste('meow_province == "Warm Temperate Northwest Atlantic"', 'meow_province == "Tropical Northwestern Atlantic"', 'meow_province == "Cold Temperate Northwest Atlantic"', sep = "|")
seq_tina_nwa <- seq_tina(X = flt)
flt <- paste('meow_province == "Mediterranean Sea"', sep = "|")
seq_tina_med <- seq_tina(X = flt)
flt <- paste('meow_province == "Northern European Seas"', sep = "|")
seq_tina_nes <- seq_tina(X = flt)
flt <- paste('meow_province == "Lusitanian"', sep = "|")
seq_tina_lus <- seq_tina(X = flt)

flt <- paste('meow_province == "Warm Temperate Northwest Atlantic"', 'meow_province == "Tropical Northwestern Atlantic"', 'meow_province == "Cold Temperate Northwest Atlantic"', sep = "|")
seq_bc_nwa <- seq_bc(X = flt)
flt <- paste('meow_province == "Mediterranean Sea"', sep = "|")
seq_bc_med <- seq_bc(X = flt)
flt <- paste('meow_province == "Northern European Seas"', sep = "|")
seq_bc_nes <- seq_bc(X = flt)
flt <- paste('meow_province == "Lusitanian"', sep = "|")
seq_bc_lus <- seq_bc(X = flt)


nwa_seq <- ggplot() +
  geom_point(data = seq_tina_nwa$seq_tina, aes(position.x, 1- distance)) +
  geom_line(data = seq_tina_nwa$seq_tina, aes(position.x, 1- distance)) +
  #geom_point(data = seq_bc_nwa$seq_bc, aes(position.x, 1- distance)) +
  geom_line(data = seq_bc_nwa$seq_bc, aes(position.x, 1- distance), linetype = 2)  +
  geom_point(data = seq_bc_nwa$LCB, aes(position, 1 - distance, size = LCBD, fill = signif), shape = 21) +
  theme_bw() +
  scale_fill_manual(values = c("#282C34", "red")) +
  ylim(c(0,1))+
  scale_size_continuous(range = c(1,3))

med_seq <- ggplot() +
  geom_point(data = seq_tina_med$seq_tina, aes(position.x, 1- distance)) +
  geom_line(data = seq_tina_med$seq_tina, aes(position.x, 1- distance)) +
  geom_line(data = seq_bc_med$seq_bc, aes(position.x, 1- distance), linetype = 2)  +
  geom_point(data = seq_bc_med$LCB, aes(position, 1 - distance, size = LCBD, fill = signif), shape = 21) +
  theme_bw() +
  scale_fill_manual(values = c("#282C34", "red")) +
  ylim(c(0,1))+
  scale_size_continuous(range = c(1,3))

nes_seq <- ggplot() +
  geom_point(data = seq_tina_nes$seq_tina, aes(position.x, 1- distance)) +
  geom_line(data = seq_tina_nes$seq_tina, aes(position.x, 1- distance)) +
  geom_line(data = seq_bc_nes$seq_bc, aes(position.x, 1- distance), linetype = 2)  +
  geom_point(data = seq_bc_nes$LCB, aes(position, 1 - distance, size = LCBD, fill = signif), shape = 21) +
  theme_bw() +
  scale_fill_manual(values = c("#282C34", "red")) +
  ylim(c(0,1))+
  scale_size_continuous(range = c(1,3))

lus_seq <- ggplot() +
  geom_point(data = seq_tina_lus$seq_tina, aes(position.x, 1- distance)) +
  geom_line(data = seq_tina_lus$seq_tina, aes(position.x, 1- distance)) +
  geom_line(data = seq_bc_lus$seq_bc, aes(position.x, 1- distance), linetype = 2)  +
  geom_point(data = seq_bc_lus$LCB, aes(position, 1 - distance, size = LCBD, fill = signif), shape = 21) +
  theme_bw() +
  scale_fill_manual(values = c("#282C34", "red")) +
  ylim(c(0,1))+
  scale_size_continuous(range = c(1,3))


ggpubr::ggarrange(nwa_seq, lus_seq, med_seq, nes_seq, common.legend = TRUE, labels = c("a", "b", "c", "d"))

ggsave(plot = last_plot(), filename = "osd2014_18S_asv/figures/osd2014_dada2_phyloseq_beta_seqBC.pdf", width = 7, height = 4)

# Distance decay ----------------------------------------------------------



osd2014_haversine_distance %>%
  inner_join(osd2014_18S_asv_bc) %>%
  filter(same_region == TRUE, meow_province.x %in% osd2014_selected_meow_provinces$meow_province) %>%
  mutate(meow_region = ifelse(grepl("east", meow_region.x), meow_province.x, meow_region.x)) %>%
  select(haversine, distance, meow_region) %>%
  ggplot(aes(haversine/1000, distance)) +
  geom_point(shape = 21, fill = "grey30", size = 2) +
  geom_smooth(size = 1) +
  facet_wrap(~meow_region, scales = "free_x") +
  theme_bw()

osd2014_haversine_distance %>%
  inner_join(osd2014_18S_asv_tina_w) %>%
  filter(same_region == TRUE, meow_province.x %in% osd2014_selected_meow_provinces$meow_province) %>%
  ggplot(aes(haversine/1000, distance)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~meow_region.x, scales = "free_x")



# Mantel & p-Mantel tests -------------------------------------------------


osd2014_18S_asv_bc_w <- phyloseq::distance(osd2014_dada2_phyloseq_beta_vst, "bray")

osd2014_coords <- osd2014_cdata %>%
  collect(n = Inf) %>%
  dplyr::select(label, start_lon, start_lat) %>%
  arrange(label) %>%
  as.data.frame() %>%
  column_to_rownames("label")

osd2014_haversine <- geosphere::distm(osd2014_coords, fun=geosphere::distHaversine)
rownames(osd2014_haversine) <- rownames(osd2014_coords)
colnames(osd2014_haversine) <- rownames(osd2014_coords)
osd2014_haversine <- as.dist(osd2014_haversine)

osd2014_18S_asv_gu_w <- as.dist(GUniFrac(otu.tab = phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst), phy_tree(osd2014_dada2_phyloseq_beta_vst))$unifracs[, , "d_0.5"])

osd2014_delta_temp <- with(osd2014_cdata %>% select(label, water_temperature), abs(outer(water_temperature, water_temperature, "-")))
dimnames(osd2014_delta_temp) <- list(osd2014_cdata$label, osd2014_cdata$label)
osd2014_delta_temp <- as.dist(osd2014_delta_temp)


calculate_mantel_bc <- function(X, Y) {
  labels <- osd2014_cdata %>% filter(eval(parse(text = X))) %>%
    filter(label %in% selected_samples) %>% .$label
  osd2014_haversine_sel <- as.matrix(osd2014_haversine)[labels, labels] %>% as.dist()
  osd2014_delta_temp_sel <- as.matrix(osd2014_delta_temp)[labels, labels] %>% as.dist()
  osd2014_bc_sel <- as.matrix(osd2014_18S_asv_bc_w)[labels, labels] %>% as.dist()
  m <- mantel(osd2014_bc_sel, osd2014_haversine_sel)
  m_p <- mantel.partial(osd2014_bc_sel, osd2014_haversine_sel, osd2014_delta_temp_sel)
  l <- Y
  if(is.null(Y)){
    m_c = mantel.correlog(D.eco = osd2014_bc_sel, D.geo = osd2014_haversine_sel/1000, r.type = "spearman")
  }else{
    osd2014_haversine_sel_r <- osd2014_haversine_sel %>%
      broom::tidy() %>%
      mutate(distance = distance/1000) %>%
      arrange((distance)) %>%
      mutate(d_class =( l * trunc(distance/l, 1))+ l ) %>%
      .$d_class %>% unique
    m_c = mantel.correlog(D.eco = osd2014_bc_sel, D.geo = osd2014_haversine_sel/1000, r.type = "spearman", break.pts = osd2014_haversine_sel_r)
  }
  list(mantel = m, mantel.partial = m_p, mantel.corr = m_c)
}


flt <- paste('meow_province == "Warm Temperate Northwest Atlantic"', 'meow_province == "Tropical Northwestern Atlantic"', 'meow_province == "Cold Temperate Northwest Atlantic"', sep = "|")
mantel_nwa <- calculate_mantel_bc(X = flt, Y = NULL)
mantel_nwa$mantel.corr %>% plot
flt <- paste('meow_province == "Mediterranean Sea"', sep = "|")
mantel_med <- calculate_mantel_bc(X = flt, Y = NULL)
mantel_med$mantel.corr %>% plot
flt <- paste('meow_province == "Northern European Seas"', sep = "|")
mantel_nes <- calculate_mantel_bc(X = flt, Y = NULL)
mantel_nes$mantel.corr %>% plot
flt <- paste('meow_province == "Lusitanian"', sep = "|")
mantel_lus <- calculate_mantel_bc(X = flt, Y = NULL)
mantel_lus$mantel.corr %>% plot



calculate_mantel_phylo <- function(X, Y) {
  labels <- osd2014_cdata %>% filter(eval(parse(text = X))) %>%
    filter(label %in% selected_samples) %>%
    .$label
  osd2014_18S_asv_gu_w_sel <- (beta.mntd.weighted)[labels, labels] %>% as.dist()
  osd2014_delta_temp_sel <- as.matrix(osd2014_delta_temp)[labels, labels] %>% as.dist()
  osd2014_bc_sel <- as.matrix(osd2014_18S_asv_bc_w)[labels, labels] %>% as.dist()
  m <- mantel(osd2014_bc_sel, osd2014_18S_asv_gu_w_sel)
  m_p <- mantel.partial(osd2014_bc_sel, osd2014_18S_asv_gu_w_sel, osd2014_delta_temp_sel)
  l <- Y
  if(is.null(Y)){
    m_c = mantel.correlog(D.eco = osd2014_bc_sel, D.geo = osd2014_18S_asv_gu_w_sel, r.type = "spearman")
  }else{
    osd2014_phylo_sel_r <- osd2014_18S_asv_gu_w_sel %>%
      broom::tidy() %>%
      arrange((distance)) %>%
      mutate(d_class =( l * trunc(distance/l, 1))+ l ) %>%
      .$d_class %>% unique
    m_c = mantel.correlog(D.eco = osd2014_bc_sel, D.geo = osd2014_18S_asv_gu_w_sel, r.type = "spearman", break.pts = osd2014_phylo_sel_r)
  }
  list(mantel = m, mantel.partial = m_p, mantel.corr = m_c)
}


flt <- paste('meow_province == "Warm Temperate Northwest Atlantic"', 'meow_province == "Tropical Northwestern Atlantic"', 'meow_province == "Cold Temperate Northwest Atlantic"', sep = "|")
mantel_nwa_phylo <- calculate_mantel_phylo(X = flt, Y = NULL)
mantel_nwa_phylo$mantel.corr %>% plot
flt <- paste('meow_province == "Mediterranean Sea"', sep = "|")
mantel_med_phylo <- calculate_mantel_phylo(X = flt, Y = NULL)
mantel_med_phylo$mantel.corr %>% plot
flt <- paste('meow_province == "Northern European Seas"', sep = "|")
mantel_nes_phylo <- calculate_mantel_phylo(X = flt, Y = NULL)
mantel_nes_phylo$mantel.corr %>% plot
flt <- paste('meow_province == "Lusitanian"', sep = "|")
mantel_lus <- calculate_mantel_bc(X = flt, Y = NULL)
mantel_lus$mantel.corr %>% plot

# BEGIN: Save objects ------------------------------------------------------------
# WARNING!!! You might not want to run this code --------------------------
save.image("osd2014_18S_asv/data/osd2014_18S_asv_beta.Rdata")
# END: Save objects ------------------------------------------------------------
