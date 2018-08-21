library(RNetCDF)
library(ncdf)
library(fields)
library(ncdf4)
mycdf <- nc_open("osd2014_ancillary_data/data/ETOPO_PRE1_1m_fer_Y3000.nc", verbose = TRUE, write = FALSE)
Fer  <- ncvar_get(mycdf, "Fer")
latmat  <- ncvar_get(mycdf, "ETOPO60Y")
longmat <- ncvar_get(mycdf, "ETOPO60X")
longmat <- longmat-20.5
deptht  <- ncvar_get(mycdf, "deptht")
time_counter <- ncvar_get(mycdf, "time_counter")

summary(Fer)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's
#0.000e+00 0.000e+00 0.000e+00 3.568e+30 0.000e+00 1.041e+36  10098420

#units: mmol/m3"
#long_name: Dissolved Iron Concentration


mapLat = data.frame(lat=latmat, map=c(1:length(latmat)))
mapLong = data.frame(long=longmat, map=c(1:length(longmat)))
mapLong$long180 =  ifelse(mapLong$long>180, (((mapLong$long + 180) %% 360) - 180), mapLong$long)

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

my_date = data.frame(local_date = c("12-06-2014","13-06-2014","17-06-2014","18-06-2014",
                                    "19-06-2014","20-06-2014","21-06-2014","22-06-2014",
                                    "23-06-2014","24-06-2014","25-06-2014","02-07-2014",
                                    "13-07-2014","19-07-2014"),
                     num = c(162,163,167,168,169,170,171,172,173,174,175,182,193,199))


myCoord <- osd2014_cdata %>% left_join(my_date)
out = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(),
                 long_gd=as.numeric(), lat_gd=as.numeric(),iron=as.numeric())

for (i in 1:nrow(myCoord)){
    osd_id = myCoord[i,1]
    month  =  as.numeric(myCoord[i,6])
    long   = as.numeric(myCoord[i,3])          #long
    long1   = ifelse(round(long)%% 2 == 0, round(long), round(long)-1)
    lat    = as.numeric(myCoord[i,2])           #lat
    lat1   = ifelse(round(lat)%% 2 == 0, round(lat), round(lat)-1)
    if(lat>lat1){lat1 = lat1+0.5}else{lat1 = lat1-0.5}

    long_gd =  mapLong$map[match(long1, mapLong$long180)]
    lat_gd  =  mapLat$map[match(lat1, mapLat$lat)]
    iron = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,month), count=c(1,1,1,1))

    #iron1 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,1), count=c(1,1,5,1))
    #iron2 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,2), count=c(1,1,5,1))
    #iron3 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,3), count=c(1,1,5,1))
    #iron4 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,4), count=c(1,1,5,1))
    #iron5 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,5), count=c(1,1,5,1))
    #iron6 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,6), count=c(1,1,1,1))
    #iron7 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,7), count=c(1,1,5,1))
    #iron8 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,8), count=c(1,1,5,1))
    #iron9 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,9), count=c(1,1,5,1))
    #iron10 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,10), count=c(1,1,5,1))
    #iron11 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,11), count=c(1,1,5,1))
    #iron12 = ncvar_get(mycdf, "Fer", start=c(long_gd,lat_gd,1,12), count=c(1,1,5,1))

    out[i,]= c(osd_id, long, lat, long_gd, lat_gd, iron)
}
out
write.table(out[,-c(4,5)], "Iron_osd_sites.csv", sep=",", quote=F, row.names=F)
library(RPostgreSQL)  # loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_iron_data"), value=out %>% dplyr::select(-lat, long),overwrite=TRUE,row.names=FALSE)









