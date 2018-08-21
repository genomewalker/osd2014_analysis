#help
#https://www.image.ucar.edu/GSP/Software/Netcdf/


library(ncdf)
library(fields)
library(chron)
library(RNetCDF)

ex.nc = open.ncdf("osd2014_ancillary_data/data/mld_DReqDTm02_c1m_reg2.0.nc", write=FALSE, readunlim=FALSE)

print(ex.nc)
# # this prints out sensible information
summary(ex.nc)


lonmat = get.var.ncdf( ex.nc, "lon")          # coordinate variable
latmat = get.var.ncdf( ex.nc, "lat")          # coordinate variable
timearr  = get.var.ncdf( ex.nc, "time")         # time variable

a = get.var.ncdf( ex.nc, "mld_raw")      #
b = get.var.ncdf( ex.nc, "n_profiles")   #
c = get.var.ncdf( ex.nc, "med_dev")      #
d = get.var.ncdf( ex.nc, "mld_smth")     #
e = get.var.ncdf( ex.nc, "mld")          #
f = get.var.ncdf( ex.nc, "krig_std_dev") #
g = get.var.ncdf( ex.nc, "mask")         #


********************************************************************************
library(reshape)
myMldDF = data.frame(a[,,6])
rownames(myMldDF) = lonmat
colnames(myMldDF) = latmat

myMldDF.melt = melt(myMldDF)
colnames(myMldDF.melt) = c("latitude","MLD")
myMldDF.melt$longitude = rep(lonmat, length(latmat))
myMldDF.melt = myMldDF.melt[, c(1,3,2)]
myMldDF.melt$longitude180 =  ifelse(myMldDF.melt$longitude>180, (((myMldDF.melt$longitude + 180) %% 360) - 180), myMldDF.melt$longitude)
head(myMldDF.melt)
dim(myMldDF.melt)
write.table(myMldDF.melt, "MLDtype3_june", sep=",", quote=F, row.names=F)


myMldDF_july = data.frame(a[,,7])
rownames(myMldDF_july) = lonmat
colnames(myMldDF_july) = latmat
myMldDF_july.melt = melt(myMldDF_july)
colnames(myMldDF_july.melt) = c("latitude","MLD")
myMldDF_july.melt$longitude = rep(lonmat, length(latmat))
myMldDF_july.melt = myMldDF_july.melt[, c(1,3,2)]
myMldDF_july.melt$longitude180 =  ifelse(myMldDF_july.melt$longitude>180, (((myMldDF_july.melt$longitude + 180) %% 360) - 180), myMldDF.melt$longitude)
write.table(myMldDF.melt, "MLDtype3_july", sep=",", quote=F, row.names=F)



#z1 = get.var.ncdf( ex.nc, "mld", start=c(50,44,6), count=c(1,1,1))
#z1

mapLat = data.frame(lat=latmat, map=c(1:length(latmat)))
mapLong = data.frame(long=lonmat, map=c(1:length(lonmat)))
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

out = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), long_gd=as.numeric(), lat_gd=as.numeric(), mld=as.numeric())

for (i in 1:nrow(myCoord)){
    osd_id = myCoord[i,1]
    time   = as.numeric(myCoord[i,6])          #june

    long   = as.numeric(myCoord[i,3])          #long
    long1   = ifelse(round(long)%% 2 == 0, round(long), round(long)-1)
    lat    = as.numeric(myCoord[i,2])           #lat
    lat1   = ifelse(round(lat)%% 2 == 0, round(lat), round(lat)-1)

# osd_id; long; lat; long_gd; lat_gd

    long_gd =  mapLong$map[match(long1, mapLong$long180)]
    lat_gd  =  mapLat$map[match(lat1, mapLat$lat)]

    mld = get.var.ncdf(ex.nc, "mld", start=c(long_gd,lat_gd,time), count=c(1,1,1))
    out[i,]= c(osd_id, long, lat, long_gd, lat_gd, round(mld,2))
}
out
write.table(out, "mld_osd_sites", sep=",", quote=F, row.names=F)


library(RPostgreSQL)  # loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_mld_data"), value=out %>% dplyr::select(-lat, long),overwrite=TRUE,row.names=FALSE)




