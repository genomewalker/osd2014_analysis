library(geosphere)
library(tidyverse)

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

osd2014_coords <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf) %>%
  dplyr::select(label, start_lon, start_lat) %>%
  arrange(label) %>%
  as.data.frame() %>%
  column_to_rownames("label")



geo_dm <- geosphere::distm(osd2014_coords, fun=geosphere::distHaversine)

rownames(geo_dm) <- rownames(osd2014_coords)
colnames(geo_dm) <- rownames(osd2014_coords)
geo_dm_df <- broom::tidy(as.dist(geo_dm)) %>% as.tibble()

drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_haversine_distance"), value=geo_dm_df,overwrite=TRUE,row.names=FALSE)
