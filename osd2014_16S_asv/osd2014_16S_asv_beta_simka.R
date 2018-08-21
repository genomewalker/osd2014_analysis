library(RPostgreSQL)  # loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_simka_k21_bc"), value=simka_k21,overwrite=TRUE,row.names=FALSE)
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)
osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)
osd2014_cdata <- osd2014_cdata %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% left_join(osd2014_meow_regions)
osd2014_16S_asv_tina_uw <- pina_tina_results[[1]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = pina_tina_results[[1]]$name) %>%
  as_tibble()
osd2014_16S_asv_tina_w <- pina_tina_results[[1]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = pina_tina_results[[2]]$name) %>%
  as_tibble()
osd2014_16S_asv_tina_w <- pina_tina_results[[2]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = pina_tina_results[[2]]$name) %>%
  as_tibble()
osd2014_16S_asv_tina_uw_filt <- pina_tina_results[[3]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = pina_tina_results[[3]]$name) %>%
  as_tibble()
osd2014_16S_asv_tina_uw_filt
osd2014_16S_asv_tina_uw <- osd2014_pina_tina_results[[1]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = pina_tina_results[[1]]$name) %>%
  as_tibble()
simka_k31_m <- simka_k31_m[rownames(osd2014_pina_tina_results[[1]]$cs), rownames(osd2014_pina_tina_results[[1]]$cs)]
simka_k31 <- simka_k31_m %>%
  as.dist() %>%
  broom::tidy() %>%
  as_tibble()
simka_k21_m <- read.table(file = "~/Downloads/allVSall_paired_k21/mat_abundance_braycurtis.csv.gz", header = T, row.names = 1, sep = ";", check.names = F)
simka_k21_m <- simka_k21_m[rownames(osd2014_pina_tina_results[[1]]$cs), rownames(osd2014_pina_tina_results[[1]]$cs)]
simka_k21 <- simka_k21_m %>%
  as.dist() %>%
  broom::tidy() %>%
  as_tibble()
library(RPostgreSQL)  # loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_simka_k21_bc"), value=simka_k21,overwrite=TRUE,row.names=FALSE)
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)
osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)
osd2014_cdata <- osd2014_cdata %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% left_join(osd2014_meow_regions)
osd2014_16S_asv_tina_uw <- osd2014_pina_tina_results[[1]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = pina_tina_results[[1]]$name) %>%
  as_tibble()
osd2014_16S_asv_tina_uw_filt <- osd2014_pina_tina_results[[3]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[3]]$name) %>%
  as_tibble()
osd2014_16S_asv_tina_uw_filt
osd2014_16S_asv_pina_w <- osd2014_pina_tina_results[[6]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[6]]$name) %>%
  as_tibble()
osd2014_16S_asv_pina_w
osd2014_16S_asv_tina_uw <- osd2014_pina_tina_results[[1]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[1]]$name) %>%
  as_tibble()
osd2014_16S_asv_tina_w <- osd2014_pina_tina_results[[2]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[2]]$name) %>%
  as_tibble()
osd2014_16S_asv_tina_uw_filt <- osd2014_pina_tina_results[[3]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[3]]$name) %>%
  as_tibble()
osd2014_16S_asv_tina_w_filt <- osd2014_pina_tina_results[[4]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[4]]$name) %>%
  as_tibble()
osd2014_16S_asv_pina_uw <- osd2014_pina_tina_results[[5]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[5]]$name) %>%
  as_tibble()
osd2014_16S_asv_pina_w <- osd2014_pina_tina_results[[6]]$cs %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = osd2014_pina_tina_results[[6]]$name) %>%
  as_tibble()
osd2014_16S_asv_bc <- distance(osd2014_dada2_phyloseq_beta_vst, "bray") %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = "Bray-Curtis") %>%
  as_tibble()
osd2014_16S_asv_jc <- distance(osd2014_dada2_phyloseq_beta_vst, method = "jaccard", binary = FALSE) %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = "Jaccard") %>%
  as_tibble()
osd2014_mitag_qiime97_phyloseq_t <- phyloseq:::veganifyOTU(osd2014_mitag_qiime97_phyloseq)
osd2014_mitag_qiime97_phyloseq_t <- osd2014_mitag_qiime97_phyloseq_t[rownames(pina_tina_results[[4]]$cs), ]
osd2014_16S_mitag_bc <- vegan::vegdist(osd2014_mitag_qiime97_phyloseq_t) %>%
  as.dist() %>%
  broom::tidy() %>%
  mutate(class = "Bray-Curtis") %>%
  as_tibble()
osd2014_simka_k31_bc <- tbl(my_db, "osd2014_simka_k31_bc") %>%
  collect(n = Inf)
osd2014_simka_k21_bc <- tbl(my_db, "osd2014_simka_k21_bc") %>%
  collect(n = Inf)
osd2014_16S_asv_tina_uw %>%
  left_join(osd2014_16S_asv_tina_w, by = c("item1", "item2")) %>% mutate(class = "TINA, unweighted vs TINA, weighted") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("TINA, unweighted") +
  ylab("TINA, weighted")
osd2014_16S_asv_tina_uw %>%
  left_join(osd2014_16S_asv_tina_uw_filt, by = c("item1", "item2")) %>% mutate(class = "TINA, unweighted vs TINA, weighted") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("TINA, unweighted") +
  ylab("TINA, weighted")
osd2014_16S_asv_tina_uw_filt %>%
  left_join(osd2014_16S_asv_tina_w_filt, by = c("item1", "item2")) %>% mutate(class = "TINA, unweighted vs TINA, weighted") %>%
  ggplot(aes(distance.x, distance.y)) +
  geom_point() +
  geom_smooth() +
  geom_abline() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("TINA, unweighted") +
  ylab("TINA, weighted")
osd2014_16S_asv_tina_w_meow <- osd2014_16S_asv_tina_w_filt %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label"))
osd2014_16S_asv_pina_w_meow <- osd2014_16S_asv_pina_w %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label"))
osd2014_16S_asv_bc_meow <- osd2014_16S_asv_bc%>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label"))
osd2014_simka_k21_bc_meow <- osd2014_simka_k21_bc%>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item1" = "label")) %>%
  left_join(osd2014_cdata %>% select(label, meow_province, meow_region), by = c("item2" = "label"))
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x ="tina", y = distance, color = same_region)) +
  geom_boxplot() +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5) +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5) +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5, errorbar.draw = TRUE, width = 0.2) +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5, errorbar.draw = TRUE, width = 0.4) +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5, errorbar.draw = TRUE, width = 0.4, jitter.width = 0.6) +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5, errorbar.draw = TRUE) +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5, errorbar.draw = TRUE, outlier.colour = NULL) +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5, errorbar.draw = TRUE, outlier.colour = NULL, outlier.fill = NULL) +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5, errorbar.draw = TRUE, outlier.size = 0) +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5, errorbar.draw = TRUE, outlier.size = NULL) +
  ggpubr::stat_compare_means()
osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x = "tina", y = distance, color = same_region, fill = same_region)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.5, color = "black", alpha = 0.5, errorbar.draw = TRUE, outlier.shape = NA ) +
  ggpubr::stat_compare_means()
bind_rows(osd2014_16S_asv_tina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            filter(same_region ==  TRUE) %>% mutate(class = "TINA weighted"),
          osd2014_16S_asv_pina_w_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            filter(same_region ==  TRUE) %>% mutate(class = "PINA unweighted"),
          osd2014_16S_asv_bc_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
            filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
            filter(same_region ==  TRUE) %>% mutate(class = "ASV BC"),
          # osd2014_simka_k21_bc_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
          #   filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
          #   filter(same_region ==  TRUE) %>% mutate(class = "SIMKA BC k = 21")
)
osd2014_16S_asv_bc_meow %>% mutate(same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x ="tina", y = distance, color = same_region)) +
  ggpol::geom_boxjitter() +
  ggpubr::stat_compare_means()
osd2014_16S_asv_bc_meow %>% mutate( same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  filter(same_region ==  TRUE) %>%
  ggplot(aes(x = meow_province.x, y = distance )) +
  geom_boxplot() +
  ggpubr::stat_compare_means()
osd2014_16S_asv_bc_meow %>% mutate(same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x ="tina", y = distance, color = same_region)) +
  ggpol::geom_boxjitter() +
  ggpubr::stat_compare_means()
osd2014_16S_asv_bc_meow %>% mutate(same_region = ifelse(meow_province.x ==  meow_province.y, TRUE, FALSE)) %>%
  filter(!(is.na(meow_region.x) | is.na(meow_region.y))) %>%
  ggplot(aes(x ="tina", y = distance, color = same_region)) +
  ggpol::geom_boxjitter() +
  ggpubr::stat_compare_means()
