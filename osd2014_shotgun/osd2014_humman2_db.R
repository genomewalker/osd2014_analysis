library(tidyverse)

humann2_mod2label <- tbl(my_db, "humann2_mod2label") %>%
  collect(n = Inf)

path_abun <- read_tsv("~/ownCloud/OSD_paper/BTEX/humann2_pathabundance.tsv", col_names = TRUE) %>%
  dplyr::rename(pathway = `# Pathway`) %>%
  gather(label, abundance, -pathway) %>%
  mutate(label = gsub("_Abundance", "", label)) %>%
  filter(pathway != "UNINTEGRATED", pathway != "UNMAPPED") %>%
  mutate(label = plyr::mapvalues(label, humann2_mod2label$mod, humann2_mod2label$orig))


path_cov <- read_tsv("~/ownCloud/OSD_paper/BTEX/humann2_pathcoverage.tsv", col_names = TRUE) %>%
  dplyr::rename(pathway = `# Pathway`) %>%
  gather(label, coverage, -pathway) %>%
  mutate(label = gsub("_Coverage", "", label)) %>%
  filter(pathway != "UNINTEGRATED", pathway != "UNMAPPED") %>%
  mutate(label = plyr::mapvalues(label, humann2_mod2label$mod, humann2_mod2label$orig))


path_relab <- read_tsv("~/ownCloud/OSD_paper/BTEX/humann2_pathabundance_relab.tsv", col_names = TRUE) %>%
  dplyr::rename(pathway = `# Pathway`) %>%
  gather(label, relabundance, -pathway) %>%
  mutate(label = gsub("_Abundance", "", label)) %>%
  filter(pathway != "UNINTEGRATED", pathway != "UNMAPPED") %>%
  mutate(label = plyr::mapvalues(label, humann2_mod2label$mod, humann2_mod2label$orig))


mod_relab <- read_tsv("~/ownCloud/OSD_paper/BTEX/humann2_modabundance_relab.tsv", col_names = TRUE) %>%
  dplyr::rename(module = `# Pathway`) %>%
  gather(label, relabundance, -module) %>%
  mutate(label = gsub("_Abundance", "", label)) %>%
  filter(module != "UNINTEGRATED", module != "UNMAPPED") %>%
  mutate(label = plyr::mapvalues(label, humann2_mod2label$mod, humann2_mod2label$orig))



library(RPostgreSQL)  # loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_humman2_mod_relab"), value=mod_relab,overwrite=TRUE,row.names=FALSE)
