library(tidyverse)

sparcc_correlation <- read_tsv("osd2014_18S_asv/data/median_correlation.tsv.gz", col_names = TRUE, comment = "|") %>%
  dplyr::rename(otu_id = `#OTU ID`)

sparcc_pvalue <- read_tsv("osd2014_18S_asv/data/pvalues.tsv.gz", col_names = TRUE, comment = "|") %>%
  dplyr::rename(otu_id = `#OTU ID`)

osd2014_sparcc <- sparcc_correlation %>%
  gather(item2, correlation, -otu_id) %>%
  rename(item1 = otu_id) %>%
  filter(item1 != item2) %>%
  inner_join(sparcc_pvalue %>%
  gather(item2, pvalue, -otu_id) %>%
  rename(item1 = otu_id) %>%
  filter(item1 != item2))

save(osd2014_sparcc, file =  "osd2014_18S_asv/data/osd2014_18S_asv_networks.Rdata")
