library(geosphere)
library(tidyverse)

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

# load("osd2014_ancillary_data/data/osd2014_haversine_distance_get.Rdata", verbose = TRUE)
# load(url("http://osd2014.metagenomics.eu/osd2014_ancillary_data/data/osd2014_haversine_distance_get.Rdata"), verbose = TRUE)

# END: WARNING!! ---------------------------------------------------------------

# Load necessary data -----------------------------------------------------
# Use if you have the postgres DB in place
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

osd2014_coords <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf) %>%
  dplyr::select(label, start_lon, start_lat) %>%
  arrange(label) %>%
  as.data.frame() %>%
  column_to_rownames("label")

# If downloaded file at osd2014_ancillary_data/data/ use:
# Basic contextual data
load("osd2014_16S_asv/data/osd2014_basic_cdata.Rdata", verbose = TRUE)
osd2014_coords <- osd2014_cdata %>%
  collect(n = Inf) %>%
  dplyr::select(label, start_lon, start_lat) %>%
  arrange(label) %>%
  as.data.frame() %>%
  column_to_rownames("label")

# If remote use
# Basic contextual data
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_basic_cdata.Rdata"), verbose = TRUE)
osd2014_coords <- osd2014_cdata %>%
  collect(n = Inf) %>%
  dplyr::select(label, start_lon, start_lat) %>%
  arrange(label) %>%
  as.data.frame() %>%
  column_to_rownames("label")


geo_dm <- geosphere::distm(osd2014_coords, fun=geosphere::distHaversine)

rownames(geo_dm) <- rownames(osd2014_coords)
colnames(geo_dm) <- rownames(osd2014_coords)
geo_dm_df <- broom::tidy(as.dist(geo_dm)) %>% as.tibble()

# BEGIN: Save objects ------------------------------------------------------------
# WARNING!!! You might not want to run this code --------------------------
save.image(file = "osd2014_ancillary_data/data/osd2014_haversine_distance_get.Rdata")
save(geo_dm_df, file = "osd2014_ancillary_data/data/osd2014_haversine_distance.Rdata")
# END: Save objects ------------------------------------------------------------


# drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
# con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
# dbWriteTable(con, c("osd_analysis", "osd2014_haversine_distance"), value=geo_dm_df,overwrite=TRUE,row.names=FALSE)
