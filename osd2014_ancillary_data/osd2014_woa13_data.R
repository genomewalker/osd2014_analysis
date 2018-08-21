library(RNetCDF)
library(ncdf)
library(fields)
library(ncdf4)
library(RPostgreSQL)
library(tidyverse)

# READ COORDINATE FILE
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  dplyr::select(label, start_lat, start_lon, local_date, water_depth) %>%
  mutate(local_date = gsub("/", "-", local_date), local_date = gsub("-14", "-2014", local_date)) %>%
  separate(local_date, into = c("day", "month", "year"), remove = FALSE) %>%
  mutate(day = as.integer(day), month = as.integer(month), year = as.integer(year))

my_date = data.frame(local_date = c("12-06-2014","13-06-2014","17-06-2014","18-06-2014",
                                    "19-06-2014","20-06-2014","21-06-2014","22-06-2014",
                                    "23-06-2014","24-06-2014","25-06-2014","02-07-2014",
                                    "13-07-2014","19-07-2014"),
                     num = c(162,163,167,168,169,170,171,172,173,174,175,182,193,199))


myCoord <- osd2014_cdata %>% left_join(my_date)
range(myCoord$start_lat)
range(myCoord$start_lon)

temp_june = read.table("osd2014_ancillary_data/data/woa13_decav_t06mn01v2.csv", header = F, fill=T, sep=",", comment.char = '#')
names(temp_june) <- c("LATITUDE", "LONGITUDE", 0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,
                      300,325,350,375,400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500)
temp_july = read.table("osd2014_ancillary_data/data/woa13_decav_t07mn01v2.csv", header = F, fill=T, sep=",", comment.char = '#')
names(temp_july) <- c("LATITUDE", "LONGITUDE", 0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,
                      300,325,350,375,400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500)

salinity_june = read.table("osd2014_ancillary_data/data/woa13_decav_s06mn01v2.csv", header = FALSE, fill=T, sep=",", comment.char = '#')
names(salinity_june) <- c("LATITUDE", "LONGITUDE",0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,
                          325,350,375,400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500)

salinity_july = read.table("osd2014_ancillary_data/data/woa13_decav_s07mn01v2.csv", header = FALSE, fill=T, sep=",", comment.char = '#')
names(salinity_july) <- c("LATITUDE", "LONGITUDE",0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,
                          325,350,375,400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500)

phosphate_june = read.csv("osd2014_ancillary_data/data/woa13_all_p06mn01.csv", header = FALSE, fill=T, sep=",", comment.char = '#')
names(phosphate_june) <- c("LATITUDE", "LONGITUDE", 0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,
                           325,350,375,400,425,450,475,500)
phosphate_july = read.csv("osd2014_ancillary_data/data/woa13_all_p07mn01.csv", header = FALSE, fill=T, sep=",", comment.char = '#')
names(phosphate_july) <- c("LATITUDE", "LONGITUDE", 0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,
                           325,350,375,400,425,450,475,500)

oxy_june = read.table("osd2014_ancillary_data/data/woa13_all_o06mn01.csv", header = FALSE, fill=T, sep=",", comment.char = '#')
names(oxy_june) <- c("LATITUDE", "LONGITUDE", 0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,325,
                     350,375,400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500)
oxy_july = read.csv("osd2014_ancillary_data/data/woa13_all_o07mn01.csv", header = FALSE, fill=T, sep=",", comment.char = '#')
names(oxy_july) <- c("LATITUDE", "LONGITUDE", 0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,325,
                     350,375,400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500)

silicate_june = read.table("osd2014_ancillary_data/data/woa13_all_i06mn01.csv", header = FALSE, fill=T, sep=",", comment.char = '#')
names(silicate_june) <- c("LATITUDE", "LONGITUDE", 0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,325,
                          350,375,400,425,450,475,500)
silicate_july = read.csv("osd2014_ancillary_data/data/woa13_all_i07mn01.csv", header = FALSE, fill=T, sep=",", comment.char = '#')
names(silicate_july) <- c("LATITUDE", "LONGITUDE", 0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,325,
                          350,375,400,425,450,475,500)

nitrate_june = read.table("osd2014_ancillary_data/data/woa13_all_n06mn01.csv", header = FALSE, fill=T, sep=",", comment.char = '#')
names(nitrate_june) <- c("LATITUDE", "LONGITUDE", 0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,325,
                         350,375,400,425,450,475,500)
nitrate_july = read.csv("osd2014_ancillary_data/data/woa13_all_n07mn01.csv", header = FALSE, fill=T, sep=",", comment.char = '#')
names(nitrate_july) <- c("LATITUDE", "LONGITUDE", 0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,325,
                         350,375,400,425,450,475,500)



#1# TEMPERATURE
out_temp = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), long_woa13=as.numeric(), lat_woa13=as.numeric(), temp=as.numeric())
for (j in 1:nrow(myCoord)){
  osd_id = myCoord[j,1]
  time   = as.numeric(myCoord[j,6])          #june
  depth  = as.numeric(myCoord[j,8])
  long   = as.numeric(myCoord[j,3])           #long
  lat    = as.numeric(myCoord[j,2])           #lat

  if(time==6){temp=temp_june}else{temp=temp_july}

  #get the lat long and appropriate depth column

  x <- as.numeric(colnames(temp)[3:ncol(temp)])
  col_id <- which.min(abs(x - as.numeric(depth))) + 2

  temp_new = temp[,c(1,2,col_id)]

  kkk = which(abs(temp_new$LATITUDE-lat)==min(abs(temp_new$LATITUDE-lat)))
  t1 = temp_new[kkk,]
  lll = which(abs(t1$LONGITUDE-long)==min(abs(t1$LONGITUDE-long)))
  t2 = t1[lll,]
  out_temp[j,]= c(osd_id, long, lat, t2[,2], t2[,1], t2[,3])
}
out_temp



#2# SALINITY
out_sal = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), long_woa13=as.numeric(), lat_woa13=as.numeric(), sal=as.numeric())
for (j in 1:nrow(myCoord)){
  osd_id = myCoord[j,1]
  time   = as.numeric(myCoord[j,6])          #june
  depth  = as.numeric(myCoord[j,8])
  long   = as.numeric(myCoord[j,3])           #long
  lat    = as.numeric(myCoord[j,2])          #lat

  if(time==6){salinity=salinity_june}else{salinity=salinity_july}

  #get the lat long and appropriate depth column
  x <- as.numeric(colnames(salinity)[3:ncol(salinity)])
  col_id <- which.min(abs(x - as.numeric(depth))) + 2
  salinity_new = salinity[,c(1,2,col_id)]

  kkk = which(abs(salinity_new$LATITUDE-lat)==min(abs(salinity_new$LATITUDE-lat)))
  t1 = salinity_new[kkk,]
  lll = which(abs(t1$LONGITUDE-long)==min(abs(t1$LONGITUDE-long)))
  t2 = t1[lll,]
  out_sal[j,]= c(osd_id, long, lat, t2[,2], t2[,1], t2[,3])
}
out_sal



#3# PHOSPHATE
out_phos = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), long_woa13=as.numeric(), lat_woa13=as.numeric(), phos=as.numeric())
for (j in 1:nrow(myCoord)){
  osd_id = myCoord[j,1]
  time   = as.numeric(myCoord[j,6])          #june
  depth  = as.numeric(myCoord[j,8])
  long   = as.numeric(myCoord[j,3])           #long
  lat    = as.numeric(myCoord[j,2])

  if(time==6){phos=phosphate_june}else{phos=phosphate_july}

  #get the lat long and appropriate depth column
  x <- as.numeric(colnames(phos)[3:ncol(phos)])
  col_id <- which.min(abs(x - as.numeric(depth))) + 2
  phos_new = phos[,c(1,2,col_id)]

  kkk = which(abs(phos_new$LATITUDE-lat)==min(abs(phos_new$LATITUDE-lat)))
  t1 = phos_new[kkk,]
  lll = which(abs(t1$LONGITUDE-long)==min(abs(t1$LONGITUDE-long)))
  t2 = t1[lll,]
  out_phos[j,]= c(osd_id, long, lat, t2[,2], t2[,1], t2[,3])
}
out_phos



#4# DISSOLVED OXYGEN
out_oxy = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), long_woa13=as.numeric(), lat_woa13=as.numeric(), oxy=as.numeric())
for (j in 1:nrow(myCoord)){
  osd_id = myCoord[j,1]
  time   = as.numeric(myCoord[j,6])          #june
  depth  = as.numeric(myCoord[j,8])
  long   = as.numeric(myCoord[j,3])           #long
  lat    = as.numeric(myCoord[j,2])

  if(time==6){oxy=oxy_june}else{oxy=oxy_july}

  #get the lat long and appropriate depth column
  x <- as.numeric(colnames(oxy)[3:ncol(oxy)])
  col_id <- which.min(abs(x - as.numeric(depth))) + 2
  oxy_new = oxy[,c(1,2,col_id)]

  kkk = which(abs(oxy_new$LATITUDE-lat)==min(abs(oxy_new$LATITUDE-lat)))
  t1 = oxy_new[kkk,]
  lll = which(abs(t1$LONGITUDE-long)==min(abs(t1$LONGITUDE-long)))
  t2 = t1[lll,]
  out_oxy[j,]= c(osd_id, long, lat, t2[,2], t2[,1], t2[,3])
}
out_oxy



#5# NITRATE
out_nit = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), long_woa13=as.numeric(), lat_woa13=as.numeric(), nitrate=as.numeric())
for (j in 1:nrow(myCoord)){
  osd_id = myCoord[j,1]
  time   = as.numeric(myCoord[j,6])          #june
  depth  = as.numeric(myCoord[j,8])
  long   = as.numeric(myCoord[j,3])           #long
  lat    = as.numeric(myCoord[j,2])

  if(time==6){nitrate=nitrate_june}else{nitrate=nitrate_july}

  #get the lat long and appropriate depth column
  x <- as.numeric(colnames(nitrate)[3:ncol(nitrate)])
  col_id <- which.min(abs(x - as.numeric(depth))) + 2
  nitrate_new = nitrate[,c(1,2,col_id)]

  kkk = which(abs(nitrate_new$LATITUDE-lat)==min(abs(nitrate_new$LATITUDE-lat)))
  t1 = nitrate_new[kkk,]
  lll = which(abs(t1$LONGITUDE-long)==min(abs(t1$LONGITUDE-long)))
  t2 = t1[lll,]
  out_nit[j,]= c(osd_id, long, lat, t2[,2], t2[,1], t2[,3])
}
out_nit



#5# SILICATE
out_sili = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), long_woa13=as.numeric(), lat_woa13=as.numeric(), silicate=as.numeric())
for (j in 1:nrow(myCoord)){
  osd_id = myCoord[j,1]
  time   = as.numeric(myCoord[j,6])          #june
  depth  = as.numeric(myCoord[j,8])
  long   = as.numeric(myCoord[j,3])           #long
  lat    = as.numeric(myCoord[j,2])

  if(time==6){silicate=silicate_june}else{silicate=silicate_july}

  #get the lat long and appropriate depth column
  x <- as.numeric(colnames(silicate)[3:ncol(silicate)])
  col_id <- which.min(abs(x - as.numeric(depth))) + 2
  silicate_new = silicate[,c(1,2,col_id)]

  kkk = which(abs(silicate_new$LATITUDE-lat)==min(abs(silicate_new$LATITUDE-lat)))
  t1 = silicate_new[kkk,]
  lll = which(abs(t1$LONGITUDE-long)==min(abs(t1$LONGITUDE-long)))
  t2 = t1[lll,]
  out_sili[j,]= c(osd_id, long, lat, t2[,2], t2[,1], t2[,3])
}
out_sili

out=cbind(out_temp %>% dplyr::select(-long, -lat), out_sal %>% dplyr::select(sal), out_nit %>% dplyr::select(nitrate), out_phos%>% dplyr::select(phos),
          out_sili%>% dplyr::select(silicate), out_oxy%>% dplyr::select(oxy))
head(out)

library(RPostgreSQL)  # loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_woa13_data"), value=out,overwrite=TRUE,row.names=FALSE)

