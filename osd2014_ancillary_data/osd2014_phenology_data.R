library(stringi)
library(tidyverse)
library(TeachingDemos)

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

# load("osd2014_ancillary_data/data/osd2014_phenology_data_get.Rdata", verbose = TRUE)
# load(url("http://osd2014.metagenomics.eu/osd2014_ancillary_data/data/osd2014_phenology_data_get.Rdata"), verbose = TRUE)

# END: WARNING!! ---------------------------------------------------------------


# BEGIN: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------

# Load necessary data -----------------------------------------------------
# Use if you have the postgres DB in place
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  dplyr::select(label, start_lat, start_lon, local_date) %>%
  mutate(local_date = gsub("/", "-", local_date), local_date = gsub("-14", "-2014", local_date)) %>%
  separate(local_date, into = c("day", "month", "year"), remove = FALSE) %>%
  mutate(day = as.integer(day), month = as.integer(month), year = as.integer(year))


# If downloaded file at osd2014_ancillary_data/data/ use:
PP_june <- read.csv("osd2014_ancillary_data/data/vgpm.2014152.all.xyz", h=T, sep=" ", fill=T)
PP_july <- read.csv("osd2014_ancillary_data/data/vgpm.2014182.all.xyz", h=T, sep=" ",, fill=T)

PP_8d_161 <- read.table("osd2014_ancillary_data/data/vgpm.2014161.all.xyz", h=T, fill=T, sep=" ")
PP_8d_169 <- read.table("osd2014_ancillary_data/data/vgpm.2014169.all.xyz", h=T, fill=T, sep=" ")
PP_8d_177 <- read.table("osd2014_ancillary_data/data/vgpm.2014177.all.xyz", h=T, fill=T, sep=" ")
PP_8d_185 <- read.table("osd2014_ancillary_data/data//vgpm.2014185.all.xyz", h=T, fill=T, sep=" ")
PP_8d_193 <- read.table("osd2014_ancillary_data/data/vgpm.2014193.all.xyz", h=T, fill=T, sep=" ")
PP_8d_201 <- read.table("osd2014_ancillary_data/data/vgpm.2014201.all.xyz", h=T, fill=T, sep=" ")


# Basic contextual data
load("osd2014_16S_asv/data/osd2014_basic_cdata.Rdata", verbose = TRUE)
osd2014_cdata <- osd2014_cdata %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  dplyr::select(label, start_lat, start_lon, local_date) %>%
  mutate(local_date = gsub("/", "-", local_date), local_date = gsub("-14", "-2014", local_date)) %>%
  separate(local_date, into = c("day", "month", "year"), remove = FALSE) %>%
  mutate(day = as.integer(day), month = as.integer(month), year = as.integer(year))
# If remote use
PP_june <- read.csv("http://osd2014.metagenomics.eu/osd2014_ancillary_data/data/vgpm.2014152.all.xyz", h=T, sep=" ", fill=T)
PP_july <- read.csv("http://osd2014.metagenomics.eu/osd2014_ancillary_data/data/vgpm.2014182.all.xyz", h=T, sep=" ",, fill=T)

PP_8d_161 <- read.table("http://osd2014.metagenomics.eu/osd2014_ancillary_data/data/vgpm.2014161.all.xyz", h=T, fill=T, sep=" ")
PP_8d_169 <- read.table("http://osd2014.metagenomics.eu/osd2014_ancillary_data/data/vgpm.2014169.all.xyz", h=T, fill=T, sep=" ")
PP_8d_177 <- read.table("http://osd2014.metagenomics.eu/osd2014_ancillary_data/data/vgpm.2014177.all.xyz", h=T, fill=T, sep=" ")
PP_8d_185 <- read.table("http://osd2014.metagenomics.eu/osd2014_ancillary_data/data//vgpm.2014185.all.xyz", h=T, fill=T, sep=" ")
PP_8d_193 <- read.table("http://osd2014.metagenomics.eu/osd2014_ancillary_data/data/vgpm.2014193.all.xyz", h=T, fill=T, sep=" ")
PP_8d_201 <- read.table("http://osd2014.metagenomics.eu/osd2014_ancillary_data/data/vgpm.2014201.all.xyz", h=T, fill=T, sep=" ")

# Basic contextual data
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_basic_cdata.Rdata"), verbose = TRUE)
osd2014_cdata <- osd2014_cdata %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  dplyr::select(label, start_lat, start_lon, local_date) %>%
  mutate(local_date = gsub("/", "-", local_date), local_date = gsub("-14", "-2014", local_date)) %>%
  separate(local_date, into = c("day", "month", "year"), remove = FALSE) %>%
  mutate(day = as.integer(day), month = as.integer(month), year = as.integer(year))
# Load necessary data -----------------------------------------------------

# END: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------



my_date = data.frame(local_date = c("12-06-2014","13-06-2014","17-06-2014","18-06-2014",
                                    "19-06-2014","20-06-2014","21-06-2014","22-06-2014",
                                    "23-06-2014","24-06-2014","25-06-2014","02-07-2014",
                                    "13-07-2014","19-07-2014"),
                     num = c(162,163,167,168,169,170,171,172,173,174,175,182,193,199))


myCoord <- osd2014_cdata %>% left_join(my_date)
range(myCoord$start_lat)
range(myCoord$start_lon)


#-------------------------------------------------------------------------------
#1# PP_MONTHLY
#-------------------------------------------------------------------------------
PP_june = PP_june[, 1:3]
dim(PP_june); head(PP_june)

PP_july = PP_july[, 1:3]
dim(PP_july); head(PP_july)


out_pp_mon = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), pp_mon_0.5=as.numeric(), pp_mon_0.1=as.numeric())
for (j in 1:nrow(myCoord)){
    osd_id = myCoord[j,1]
    time   = as.numeric(myCoord[j,6] )          #june
    long   = as.numeric(myCoord[j,3])           #long
    lat    = as.numeric(myCoord[j,2])         #lat

  if(time==6){temp=PP_june}else{temp=PP_july}

    #kkk = which(abs(temp$lat-lat)==min(abs(temp$lat-lat)))
    #t1 = temp[kkk,]
    #lll = which(abs(t1$lon-long)==min(abs(t1$lon-long)))
    #t2 = t1[lll,]
    #out_pp_mon[j,]= c(osd_id, long, lat, t2[,1], t2[,2], t2[,3])

   #EXTRAPOLATION at 0.5 degrees (50 kms) radius
   longMin = as.numeric(long - 0.5)
   longMax = as.numeric(long + 0.5)
   latMin  = as.numeric(lat - 0.5)
   latMax  = as.numeric(lat + 0.5)

   temp_grid1 = data.frame(temp[temp$lat>latMin,])
   temp_grid2 = data.frame(temp_grid1[temp_grid1$lat<latMax,])
   temp_grid3 = data.frame(temp_grid2[temp_grid2$lon>longMin,])
   temp_grid4 = data.frame(temp_grid3[temp_grid3$lon<longMax,])

   longMin0.1 = as.numeric(long - 0.1)
   longMax0.1 = as.numeric(long + 0.1)
   latMin0.1  = as.numeric(lat - 0.1)
   latMax0.1 = as.numeric(lat + 0.1)

   temp_grid5 = data.frame(temp[temp$lat>latMin0.1,])
   temp_grid6 = data.frame(temp_grid5[temp_grid5$lat<latMax0.1,])
   temp_grid7 = data.frame(temp_grid6[temp_grid6$lon>longMin0.1,])
   temp_grid8 = data.frame(temp_grid7[temp_grid7$lon<longMax0.1,])

   out_pp_mon[j,]= c(osd_id, long, lat, max(temp_grid4$value), max(temp_grid8$value[!is.na(temp_grid8$value)]))
}
out_pp_mon
out_pp_mon[out_pp_mon== -9999] <-NA


#-------------------------------------------------------------------------------
#2# PP_8d
#-------------------------------------------------------------------------------
PP_8d_161 = PP_8d_161[, 1:3]
PP_8d_169 = PP_8d_169[, 1:3]
PP_8d_177 = PP_8d_177[, 1:3]
PP_8d_185 = PP_8d_185[, 1:3]
PP_8d_193 = PP_8d_193[, 1:3]
PP_8d_201 = PP_8d_201[, 1:3]

# 12 june = 162th day  :: #19th july = 199th day
# 11 june - 161
# 19 june - 169
# 27 june - 177
# 5  july - 185
# 13 july - 193
# 21 july - 201
head(PP_8d_169)

# previously the values were chosen based on the nearest coordainate.
#EXTRAPOLATION will be done and grid will be selected around the OSD site and maximal value will be listed.

out_pp_8d = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), pp_8d_0.5=as.numeric(), pp_8d_0.1=as.numeric())
for (j in 1:nrow(myCoord)){

    osd_id = myCoord[j,1]
    num   = as.numeric(myCoord[j,8])           #june
    long   = as.numeric(myCoord[j,3])           #long
    lat    = as.numeric(myCoord[j,2])           #lat


    if(161  %<% num %<% 169){temp=PP_8d_161}else{
      if(168 %<% num %<% 177){temp=PP_8d_169}else{
        if(176 %<% num %<% 185){temp=PP_8d_177}else{
          if(184 %<% num %<% 193){temp=PP_8d_185}else{
            if(192 %<% num %<% 201){temp=PP_8d_193}else{temp=PP_8d_201}}}}}

  #OLD METHOD::
    #kkk = which(abs(temp$lat-lat)==min(abs(temp$lat-lat)))
    #t1 = temp[kkk,]
    #lll = which(abs(t1$lon-long)==min(abs(t1$lon-long)))
    #t2 = t1[lll,]
    #out_pp_8d[j,]= c(osd_id, long, lat, t2[,1], t2[,2], t2[,3])

   #EXTRAPOLATION at 0.5 degrees (50 kms) radius
   longMin = as.numeric(long - 0.5)
   longMax = as.numeric(long + 0.5)
   latMin  = as.numeric(lat - 0.5)
   latMax  = as.numeric(lat + 0.5)

   temp_grid1 = data.frame(temp[temp$lat>latMin,])
   temp_grid2 = data.frame(temp_grid1[temp_grid1$lat<latMax,])
   temp_grid3 = data.frame(temp_grid2[temp_grid2$lon>longMin,])
   temp_grid4 = data.frame(temp_grid3[temp_grid3$lon<longMax,])

   longMin0.1 = as.numeric(long - 0.1)
   longMax0.1 = as.numeric(long + 0.1)
   latMin0.1  = as.numeric(lat - 0.1)
   latMax0.1 = as.numeric(lat + 0.1)

   temp_grid5 = data.frame(temp[temp$lat>latMin0.1,])
   temp_grid6 = data.frame(temp_grid5[temp_grid5$lat<latMax0.1,])
   temp_grid7 = data.frame(temp_grid6[temp_grid6$lon>longMin0.1,])
   temp_grid8 = data.frame(temp_grid7[temp_grid7$lon<longMax0.1,])

   out_pp_8d[j,]= c(osd_id, long, lat, max(temp_grid4$value), max(temp_grid8$value[!is.na(temp_grid8$value)]))

}
out_pp_8d
out_pp_8d[out_pp_8d== -9999] <-NA

priProOSD <- out_pp_8d %>% as_tibble() %>%
  left_join(out_pp_mon %>% as_tibble()) %>%
  dplyr::select(-lat, -long)
write.table(priProOSD, "osd2014_ancillary_data/data/priProOSD.csv", sep=",", quote=F, row.names=F)

# library(RPostgreSQL)  # loads the PostgreSQL driver
# drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
# con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
# dbWriteTable(con, c("osd_analysis", "osd2014_phenology_data"), value=priProOSD,overwrite=TRUE,row.names=FALSE)


# BEGIN: Save objects ------------------------------------------------------------
# WARNING!!! You might not want to run this code --------------------------
save.image(file = "osd2014_ancillary_data/data/osd2014_phenology_data_get.Rdata")
save(priProOSD, file = "osd2014_ancillary_data/data/osd2014_phenology_data.Rdata")
# END: Save objects ------------------------------------------------------------









