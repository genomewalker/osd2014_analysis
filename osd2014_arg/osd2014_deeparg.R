library(tidyverse)
library(DESeq2)
library(tidyverse)
library(phyloseq)
library(caret)
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
load(file = "osd2014_16S_asv/data/osd2014_mitag_asv_physeq_filt_objects.Rdata", verbose = TRUE)

mitag_counts <- sample_sums(osd2014_mitag_phyloseq_alpha) %>% as.data.frame() %>% as_tibble(rownames = "label")
names(mitag_counts) <- c("label", "counts")
osd2014_deeparg <- read_tsv("~/Downloads/osd2014_deeparg_1.tsv", col_names = F) %>%
  select(X1, X6, X8)
names(osd2014_deeparg) <- c("label", "class", "prob")

osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)

osd2014_cardb <- tbl(my_db, "osd2014_cardb_mmseqs_reads") %>%
  collect(n = Inf) %>%
  inner_join(tbl(my_db, "cardb_metadata") %>%
  collect(n = Inf))

osd2014_cardb_filt <- osd2014_cardb %>%
  mutate(qcov = (qend-qstart+1)/qlen) %>%
  select(label, id, qcov, amr_gene_family, drug_class, model_name, aro_name) %>%
  filter(id > 0.8)

osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)

osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf) %>%
  filter(meow_province %in% osd2014_meow_regions$meow_province, label %in% osd2014_amp_mg_intersect$label)

osd2014_reads_num_orfs <- tbl(my_db, "osd2014_reads_num_orfs") %>%
  collect(n = Inf)
osd2014_cardb_filt %>% View
  inner_join(osd2014_cdata %>% select(label, meow_province)) %>%
  inner_join(osd2014_reads_num_orfs) %>%
  filter(label %in% osd2014_cdata$label) %>%
  inner_join(mitag_counts) %>%
  group_by(label, n_orfs, counts, meow_province) %>%
  dplyr::summarise(N=n()) %>% mutate(prop = N/n_orfs) %>%
  ggplot(aes(meow_province, prop)) +
  ggpol::geom_boxjitter() +scale_y_continuous(labels = scales::percent)

osd2014_deeparg %>%
  filter(prob >= 0.8) %>%
  group_by(label) %>%
  summarise(N = n()) %>%
  inner_join(osd2014_reads_num_orfs) %>%
  mutate(prop = N/n_orfs) %>%
  filter(label %in% osd2014_cdata$label) %>%
  inner_join(osd2014_cdata %>% select(label, meow_province)) %>%
  inner_join(mitag_counts) %>%
  ungroup() %>%
  ggplot(aes(meow_province, prop)) +
  ggpol::geom_boxjitter() +scale_y_continuous(labels = scales::percent)

