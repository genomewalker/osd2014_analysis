library(RNetCDF)
library(ncdf)
library(fields)
library(ncdf4)


# READ COORDINATE FILE
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
range(myCoord$start_lat)
range(myCoord$start_lon)


####MONTHLY####

#1# READ CHL.NC FILES
mycdf_june <- nc_open("osd2014_ancillary_data/data/A20141522014181.L3m_MO_CHL_chlor_a_4km.nc.1", verbose = FALSE, write = FALSE)
mycdf_july <- nc_open("osd2014_ancillary_data/data/A20141822014212.L3m_MO_CHL_chlor_a_4km.nc", verbose = TRUE, write = FALSE)
print(mycdf_june)
lat <- ncvar_get(mycdf_june,'lat')
lon <- ncvar_get(mycdf_june,'lon')

library(raster)
tmpin6 <- raster("osd2014_ancillary_data/data/A20141522014181.L3m_MO_CHL_chlor_a_4km.nc")
tmpin7 <- raster("osd2014_ancillary_data/data/A20141822014212.L3m_MO_CHL_chlor_a_4km.nc")
out = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), my_var=as.numeric())    #VAR
colnames(out)[4] = "chlor_a_monthly"

for (i in 1:nrow(myCoord)){
  osd_id   = myCoord[i,1]
  mylat    = as.numeric(myCoord[i,2] )          #lat
  mylong   = as.numeric(myCoord[i,3])           #long
  month    = as.numeric(myCoord[i,6] )          #month

  longMin = as.numeric(mylong - 0.1)
  longMax = as.numeric(mylong + 0.1)
  latMin  = as.numeric(mylat - 0.1)
  latMax  = as.numeric(mylat + 0.1)

  temp_grid1 = lat[lat>latMin]
  temp_grid2 = temp_grid1[temp_grid1<latMax]
  temp_grid3 = lon[lon>longMin]
  temp_grid4 = temp_grid3[temp_grid3<longMax]

  latWindx=c()
  for(k in 1:length(temp_grid2)){
    t1 = which.min(abs(lat - temp_grid2[k]))
    latWindx = c(latWindx, t1)
  }
  lonWindx=c()
  for(k in 1:length(temp_grid4)){
    t1 = which.min(abs(lon - temp_grid4[k]))
    lonWindx = c(lonWindx, t1)
  }

  if(month == 6){
    vars    = tmpin6[latWindx,lonWindx]
    var_val = max(vars[!is.na(vars)])
    if(var_val == "-Inf"){var_val="NA"}
  }else{
    vars    = tmpin7[latWindx,lonWindx]
    var_val = max(vars[!is.na(vars)])
    if(var_val == "-Inf"){var_val="NA"}
  }

  out[i,]= c(osd_id, mylong, mylat, var_val)
  cat(paste("chlor_a_monthly", i, ":", var_val, "\n"))
  flush.console()
}
print(out)


#
#2# READ PAR.NC FILES
#
mycdf_june <- nc_open("osd2014_ancillary_data/data/A20141522014181.L3m_MO_PAR_par_4km.nc", verbose = TRUE, write = FALSE)
mycdf_july <- nc_open("osd2014_ancillary_data/data/A20141822014212.L3m_MO_PAR_par_4km.nc", verbose = TRUE, write = FALSE)
print(mycdf_june)
lat <- ncvar_get(mycdf_june,'lat')
lon <- ncvar_get(mycdf_june,'lon')

#library(raster)
tmpin6 <- raster("osd2014_ancillary_data/data/A20141522014181.L3m_MO_PAR_par_4km.nc")
tmpin7 <- raster("osd2014_ancillary_data/data/A20141822014212.L3m_MO_PAR_par_4km.nc")

for (i in 1:nrow(myCoord)){
  osd_id   = myCoord[i,1]
  mylat    = as.numeric(myCoord[i,2] )          #lat
  mylong   = as.numeric(myCoord[i,3])           #long
  month    = as.numeric(myCoord[i,6] )          #month

  longMin = as.numeric(mylong - 0.1)
  longMax = as.numeric(mylong + 0.1)
  latMin  = as.numeric(mylat - 0.1)
  latMax  = as.numeric(mylat + 0.1)

  temp_grid1 = lat[lat>latMin]
  temp_grid2 = temp_grid1[temp_grid1<latMax]
  temp_grid3 = lon[lon>longMin]
  temp_grid4 = temp_grid3[temp_grid3<longMax]

  latWindx=c()
  for(k in 1:length(temp_grid2)){
    t1 = which.min(abs(lat - temp_grid2[k]))
    latWindx = c(latWindx, t1)
  }
  lonWindx=c()
  for(k in 1:length(temp_grid4)){
    t1 = which.min(abs(lon - temp_grid4[k]))
    lonWindx = c(lonWindx, t1)
  }

  if(month == 6){
    vars    = tmpin6[latWindx,lonWindx]
    var_val = max(vars[!is.na(vars)])
    if(var_val == "-Inf"){var_val="NA"}
  }else{
    vars    = tmpin7[latWindx,lonWindx]
    var_val = max(vars[!is.na(vars)])
    if(var_val == "-Inf"){var_val="NA"}
  }

  out[i,5]= var_val
  print(paste("par_monthly:", i))
  flush.console()
}
colnames(out)[5] = "par_monthly"
print(out)

#
#3# READ POC.NC FILES
#
mycdf_june <- nc_open("osd2014_ancillary_data/data/A20141522014181.L3m_MO_POC_poc_4km.nc", verbose = TRUE, write = FALSE)
mycdf_july <- nc_open("osd2014_ancillary_data/data/A20141822014212.L3m_MO_POC_poc_4km.nc", verbose = TRUE, write = FALSE)
print(mycdf_june)
lat <- ncvar_get(mycdf_june,'lat')
lon <- ncvar_get(mycdf_june,'lon')

#library(raster)
tmpin6 <- raster("osd2014_ancillary_data/data/A20141522014181.L3m_MO_POC_poc_4km.nc")
tmpin7 <- raster("osd2014_ancillary_data/data/A20141822014212.L3m_MO_POC_poc_4km.nc")


for (i in 1:nrow(myCoord)){
  osd_id   = myCoord[i,1]
  mylat    = as.numeric(myCoord[i,2] )          #lat
  mylong   = as.numeric(myCoord[i,3])           #long
  month    = as.numeric(myCoord[i,6] )          #month

  longMin = as.numeric(mylong - 0.1)
  longMax = as.numeric(mylong + 0.1)
  latMin  = as.numeric(mylat - 0.1)
  latMax  = as.numeric(mylat + 0.1)

  temp_grid1 = lat[lat>latMin]
  temp_grid2 = temp_grid1[temp_grid1<latMax]
  temp_grid3 = lon[lon>longMin]
  temp_grid4 = temp_grid3[temp_grid3<longMax]

  latWindx=c()
  for(k in 1:length(temp_grid2)){
    t1 = which.min(abs(lat - temp_grid2[k]))
    latWindx = c(latWindx, t1)
  }
  lonWindx=c()
  for(k in 1:length(temp_grid4)){
    t1 = which.min(abs(lon - temp_grid4[k]))
    lonWindx = c(lonWindx, t1)
  }

  if(month == 6){
    vars    = tmpin6[latWindx,lonWindx]
    var_val = max(vars[!is.na(vars)])
    if(var_val == "-Inf"){var_val="NA"}
  }else{
    vars    = tmpin7[latWindx,lonWindx]
    var_val = max(vars[!is.na(vars)])
    if(var_val == "-Inf"){var_val="NA"}
  }

  out[i,6] = var_val
  print(paste("poc_monthly:", i))
  flush.console()
}
colnames(out)[6] = "poc_monthly"
print(out)


#
#4# READ KD490.NC FILES
#
mycdf_june <- nc_open("osd2014_ancillary_data/data/A20141522014181.L3m_MO_KD490_Kd_490_4km.nc", verbose = TRUE, write = FALSE)
mycdf_july <- nc_open("osd2014_ancillary_data/data/A20141822014212.L3m_MO_KD490_Kd_490_4km.nc", verbose = TRUE, write = FALSE)
print(mycdf_june)
lat <- ncvar_get(mycdf_june,'lat')
lon <- ncvar_get(mycdf_june,'lon')

#library(raster)
tmpin6 <- raster("osd2014_ancillary_data/data/A20141522014181.L3m_MO_KD490_Kd_490_4km.nc")
tmpin7 <- raster("osd2014_ancillary_data/data/A20141822014212.L3m_MO_KD490_Kd_490_4km.nc")

for (i in 1:nrow(myCoord)){
  osd_id   = myCoord[i,1]
  mylat    = as.numeric(myCoord[i,2] )          #lat
  mylong   = as.numeric(myCoord[i,3])           #long
  month    = as.numeric(myCoord[i,6] )          #month
  longMin = as.numeric(mylong - 0.1)
  longMax = as.numeric(mylong + 0.1)
  latMin  = as.numeric(mylat - 0.1)
  latMax  = as.numeric(mylat + 0.1)

  temp_grid1 = lat[lat>latMin]
  temp_grid2 = temp_grid1[temp_grid1<latMax]
  temp_grid3 = lon[lon>longMin]
  temp_grid4 = temp_grid3[temp_grid3<longMax]

  latWindx=c()
  for(k in 1:length(temp_grid2)){
    t1 = which.min(abs(lat - temp_grid2[k]))
    latWindx = c(latWindx, t1)
  }
  lonWindx=c()
  for(k in 1:length(temp_grid4)){
    t1 = which.min(abs(lon - temp_grid4[k]))
    lonWindx = c(lonWindx, t1)
  }

  if(month == 6){
    vars    = tmpin6[latWindx,lonWindx]
    var_val = max(vars[!is.na(vars)])
    if(var_val == "-Inf"){var_val="NA"}
  }else{
    vars    = tmpin7[latWindx,lonWindx]
    var_val = max(vars[!is.na(vars)])
    if(var_val == "-Inf"){var_val="NA"}
  }

  out[i,7] = var_val
  print(paste("KD490:", i))
  flush.console()
}
colnames(out)[7] = "KD490_monthly"
print(out)


#
#5# READ FLH.NC FILES
#
mycdf_june <- nc_open("osd2014_ancillary_data/data/A20141522014181.L3m_MO_FLH_ipar_4km.nc", verbose = TRUE, write = FALSE)
mycdf_july <- nc_open("osd2014_ancillary_data/data/A20141822014212.L3m_MO_FLH_ipar_4km.nc", verbose = TRUE, write = FALSE)
print(mycdf_june)
lat <- ncvar_get(mycdf_june,'lat')
lon <- ncvar_get(mycdf_june,'lon')

#library(raster)
tmpin6 <- raster("osd2014_ancillary_data/data/A20141522014181.L3m_MO_FLH_ipar_4km.nc")
tmpin7 <- raster("osd2014_ancillary_data/data/A20141822014212.L3m_MO_FLH_ipar_4km.nc")

for (i in 1:nrow(myCoord)){
  osd_id   = myCoord[i,1]
  mylat    = as.numeric(myCoord[i,2] )          #lat
  mylong   = as.numeric(myCoord[i,3])           #long
  month    = as.numeric(myCoord[i,6] )          #month
  longMin = as.numeric(mylong - 0.1)
  longMax = as.numeric(mylong + 0.1)
  latMin  = as.numeric(mylat - 0.1)
  latMax  = as.numeric(mylat + 0.1)

  temp_grid1 = lat[lat>latMin]
  temp_grid2 = temp_grid1[temp_grid1<latMax]
  temp_grid3 = lon[lon>longMin]
  temp_grid4 = temp_grid3[temp_grid3<longMax]

  latWindx=c()
  for(k in 1:length(temp_grid2)){
    t1 = which.min(abs(lat - temp_grid2[k]))
    latWindx = c(latWindx, t1)
  }
  lonWindx=c()
  for(k in 1:length(temp_grid4)){
    t1 = which.min(abs(lon - temp_grid4[k]))
    lonWindx = c(lonWindx, t1)
  }

  if(month == 6){
    vars    = tmpin6[latWindx,lonWindx]
    var_val = max(vars[!is.na(vars)])
    if(var_val == "-Inf"){var_val="NA"}
  }else{
    vars    = tmpin7[latWindx,lonWindx]
    var_val = max(vars[!is.na(vars)])
    if(var_val == "-Inf"){var_val="NA"}
  }

  out[i,8] = var_val
  print(paste("FLH:", i))
  flush.console()
}
colnames(out)[8] = "FLH_monthly"
print(out)



tmpin6 <- raster("osd2014_ancillary_data/data/A20141522014181.L3m_MO_CHL_chlor_a_4km.nc")
tmpin6
e <- extent(-18040095, -18040000, -9020047, -9020000)
r <- crop(tmpin6, e)

####8Day####
library(stringi)
library(TeachingDemos)
#
#1# READ CHL.NC FILES
#
mycdf_8d_1 <- nc_open("osd2014_ancillary_data/data/A20151612015168.L3m_8D_CHL_chlor_a_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_2 <- nc_open("osd2014_ancillary_data/data/A20151692015176.L3m_8D_CHL_chlor_a_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_3 <- nc_open("osd2014_ancillary_data/data/A20151772015184.L3m_8D_CHL_chlor_a_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_4 <- nc_open("osd2014_ancillary_data/data/A20151852015192.L3m_8D_CHL_chlor_a_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_5 <- nc_open("osd2014_ancillary_data/data/A20151932015200.L3m_8D_CHL_chlor_a_4km.nc", verbose = TRUE, write = FALSE)

lat <- ncvar_get(mycdf_8d_1,'lat')
lon <- ncvar_get(mycdf_8d_1,'lon')

#library(raster)
tmpin_8d_1 <- raster("osd2014_ancillary_data/data/A20151612015168.L3m_8D_CHL_chlor_a_4km.nc")
tmpin_8d_2 <- raster("osd2014_ancillary_data/data/A20151692015176.L3m_8D_CHL_chlor_a_4km.nc")
tmpin_8d_3 <- raster("osd2014_ancillary_data/data/A20151772015184.L3m_8D_CHL_chlor_a_4km.nc")
tmpin_8d_4 <- raster("osd2014_ancillary_data/data/A20151852015192.L3m_8D_CHL_chlor_a_4km.nc")
tmpin_8d_5 <- raster("osd2014_ancillary_data/data/A20151932015200.L3m_8D_CHL_chlor_a_4km.nc")

out_8d = data.frame(osd_id=as.numeric(), long=as.numeric(), lat=as.numeric(), my_var=as.numeric())    #VAR
colnames(out_8d)[4] = "chlor_a_8d"

for (i in 1:nrow(myCoord)){
  osd_id   = myCoord[i,1]
  mylat    = as.numeric(myCoord[i,2] )          #lat
  mylon   = as.numeric(myCoord[i,3])           #long
  num      = as.numeric(myCoord[i,8])           #june

  if(161  %<% num %<% 169){temp=tmpin_8d_1}else{
    if(168 %<% num %<% 177){temp=tmpin_8d_2}else{
      if(176 %<% num %<% 185){temp=tmpin_8d_3}else{
        if(184 %<% num %<% 193){temp=tmpin_8d_4}else{temp=tmpin_8d_5}}}}

  longMin = as.numeric(mylon - 0.1)
  longMax = as.numeric(mylon + 0.1)
  latMin  = as.numeric(mylat - 0.1)
  latMax  = as.numeric(mylat + 0.1)

  temp_grid1 = lat[lat>latMin]
  temp_grid2 = temp_grid1[temp_grid1<latMax]
  temp_grid3 = lon[lon>longMin]
  temp_grid4 = temp_grid3[temp_grid3<longMax]

  latWindx=c()
  for(k in 1:length(temp_grid2)){
    t1 = which.min(abs(lat - temp_grid2[k]))
    latWindx = c(latWindx, t1)
  }
  lonWindx=c()
  for(k in 1:length(temp_grid4)){
    t1 = which.min(abs(lon - temp_grid4[k]))
    lonWindx = c(lonWindx, t1)
  }

  vars    = temp[latWindx,lonWindx]
  var_val = max(vars[!is.na(vars)])
  if(var_val == "-Inf"){var_val="NA"}

  out_8d[i,]= c(osd_id, mylon, mylat, var_val)
  print(paste("chlor_a_8d:", i, var_val))
  flush.console()
}
print(out_8d)



#
#2# READ PAR.NC FILES
#
mycdf_8d_1 <- nc_open("osd2014_ancillary_data/data/A20151612015168.L3m_8D_PAR_par_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_2 <- nc_open("osd2014_ancillary_data/data/A20151692015176.L3m_8D_PAR_par_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_3 <- nc_open("osd2014_ancillary_data/data/A20151772015184.L3m_8D_PAR_par_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_4 <- nc_open("osd2014_ancillary_data/data/A20151852015192.L3m_8D_PAR_par_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_5 <- nc_open("osd2014_ancillary_data/data/A20151932015200.L3m_8D_PAR_par_4km.nc", verbose = TRUE, write = FALSE)

lat <- ncvar_get(mycdf_8d_1,'lat')
lon <- ncvar_get(mycdf_8d_1,'lon')

#library(raster)
tmpin_8d_1 <- raster("osd2014_ancillary_data/data/A20151612015168.L3m_8D_PAR_par_4km.nc")
tmpin_8d_2 <- raster("osd2014_ancillary_data/data/A20151692015176.L3m_8D_PAR_par_4km.nc")
tmpin_8d_3 <- raster("osd2014_ancillary_data/data/A20151772015184.L3m_8D_PAR_par_4km.nc")
tmpin_8d_4 <- raster("osd2014_ancillary_data/data/A20151852015192.L3m_8D_PAR_par_4km.nc")
tmpin_8d_5 <- raster("osd2014_ancillary_data/data/A20151932015200.L3m_8D_PAR_par_4km.nc")

for (i in 1:nrow(myCoord)){
  osd_id   = myCoord[i,1]
  mylat    = as.numeric(myCoord[i,2] )          #lat
  mylon   = as.numeric(myCoord[i,3])           #long
  num      = as.numeric(myCoord[i,8])           #june

  if(161  %<% num %<% 169){temp=tmpin_8d_1}else{
    if(168 %<% num %<% 177){temp=tmpin_8d_2}else{
      if(176 %<% num %<% 185){temp=tmpin_8d_3}else{
        if(184 %<% num %<% 193){temp=tmpin_8d_4}else{temp=tmpin_8d_5}}}}

  longMin = as.numeric(mylon - 0.1)
  longMax = as.numeric(mylon + 0.1)
  latMin  = as.numeric(mylat - 0.1)
  latMax  = as.numeric(mylat + 0.1)

  temp_grid1 = lat[lat>latMin]
  temp_grid2 = temp_grid1[temp_grid1<latMax]
  temp_grid3 = lon[lon>longMin]
  temp_grid4 = temp_grid3[temp_grid3<longMax]

  latWindx=c()
  for(k in 1:length(temp_grid2)){
    t1 = which.min(abs(lat - temp_grid2[k]))
    latWindx = c(latWindx, t1)
  }
  lonWindx=c()
  for(k in 1:length(temp_grid4)){
    t1 = which.min(abs(lon - temp_grid4[k]))
    lonWindx = c(lonWindx, t1)
  }

  vars    = temp[latWindx,lonWindx]
  var_val = max(vars[!is.na(vars)])
  if(var_val == "-Inf"){var_val="NA"}

  out_8d[i,5]= var_val
  print(paste("par_8d:", i))
  flush.console()
}
colnames(out_8d)[5] = "par_8d"
print(out_8d)


#
#3# READ POC.NC FILES
#
mycdf_8d_1 <- nc_open("osd2014_ancillary_data/data/A20151612015168.L3m_8D_POC_poc_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_2 <- nc_open("osd2014_ancillary_data/data/A20151692015176.L3m_8D_POC_poc_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_3 <- nc_open("osd2014_ancillary_data/data/A20151772015184.L3m_8D_POC_poc_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_4 <- nc_open("osd2014_ancillary_data/data/A20151852015192.L3m_8D_POC_poc_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_5 <- nc_open("osd2014_ancillary_data/data/A20151932015200.L3m_8D_POC_poc_4km.nc", verbose = TRUE, write = FALSE)

lat <- ncvar_get(mycdf_8d_1,'lat')
lon <- ncvar_get(mycdf_8d_1,'lon')

#library(raster)
tmpin_8d_1 <- raster("osd2014_ancillary_data/data/A20151612015168.L3m_8D_POC_poc_4km.nc")
tmpin_8d_2 <- raster("osd2014_ancillary_data/data/A20151692015176.L3m_8D_POC_poc_4km.nc")
tmpin_8d_3 <- raster("osd2014_ancillary_data/data/A20151772015184.L3m_8D_POC_poc_4km.nc")
tmpin_8d_4 <- raster("osd2014_ancillary_data/data/A20151852015192.L3m_8D_POC_poc_4km.nc")
tmpin_8d_5 <- raster("osd2014_ancillary_data/data/A20151932015200.L3m_8D_POC_poc_4km.nc")

for (i in 1:nrow(myCoord)){
  osd_id   = myCoord[i,1]
  mylat    = as.numeric(myCoord[i,2] )          #lat
  mylon   = as.numeric(myCoord[i,3])           #long
  num      = as.numeric(myCoord[i,8])           #june

  if(161  %<% num %<% 169){temp=tmpin_8d_1}else{
    if(168 %<% num %<% 177){temp=tmpin_8d_2}else{
      if(176 %<% num %<% 185){temp=tmpin_8d_3}else{
        if(184 %<% num %<% 193){temp=tmpin_8d_4}else{temp=tmpin_8d_5}}}}

  longMin = as.numeric(mylon - 0.1)
  longMax = as.numeric(mylon + 0.1)
  latMin  = as.numeric(mylat - 0.1)
  latMax  = as.numeric(mylat + 0.1)

  temp_grid1 = lat[lat>latMin]
  temp_grid2 = temp_grid1[temp_grid1<latMax]
  temp_grid3 = lon[lon>longMin]
  temp_grid4 = temp_grid3[temp_grid3<longMax]

  latWindx=c()
  for(k in 1:length(temp_grid2)){
    t1 = which.min(abs(lat - temp_grid2[k]))
    latWindx = c(latWindx, t1)
  }
  lonWindx=c()
  for(k in 1:length(temp_grid4)){
    t1 = which.min(abs(lon - temp_grid4[k]))
    lonWindx = c(lonWindx, t1)
  }

  vars    = temp[latWindx,lonWindx]
  var_val = max(vars[!is.na(vars)])
  if(var_val == "-Inf"){var_val="NA"}

  out_8d[i,6]= var_val
  print(paste("poc_8d:", i))
  flush.console()
}
colnames(out_8d)[6] = "poc_8d"
print(out_8d)


#
#4# READ POC.NC FILES
#
mycdf_8d_1 <- nc_open("osd2014_ancillary_data/data/A20151612015168.L3m_8D_KD490_Kd_490_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_2 <- nc_open("osd2014_ancillary_data/data/A20151692015176.L3m_8D_KD490_Kd_490_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_3 <- nc_open("osd2014_ancillary_data/data/A20151772015184.L3m_8D_KD490_Kd_490_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_4 <- nc_open("osd2014_ancillary_data/data/A20151852015192.L3m_8D_KD490_Kd_490_4km.nc", verbose = TRUE, write = FALSE)
mycdf_8d_5 <- nc_open("osd2014_ancillary_data/data/A20151932015200.L3m_8D_KD490_Kd_490_4km.nc", verbose = TRUE, write = FALSE)

lat <- ncvar_get(mycdf_8d_1,'lat')
lon <- ncvar_get(mycdf_8d_1,'lon')

#library(raster)
tmpin_8d_1 <- raster("osd2014_ancillary_data/data/A20151612015168.L3m_8D_KD490_Kd_490_4km.nc")
tmpin_8d_2 <- raster("osd2014_ancillary_data/data/A20151692015176.L3m_8D_KD490_Kd_490_4km.nc")
tmpin_8d_3 <- raster("osd2014_ancillary_data/data/A20151772015184.L3m_8D_KD490_Kd_490_4km.nc")
tmpin_8d_4 <- raster("osd2014_ancillary_data/data/A20151852015192.L3m_8D_KD490_Kd_490_4km.nc")
tmpin_8d_5 <- raster("osd2014_ancillary_data/data/A20151932015200.L3m_8D_KD490_Kd_490_4km.nc")

for (i in 1:nrow(myCoord)){
  osd_id   = myCoord[i,1]
  mylat    = as.numeric(myCoord[i,2] )          #lat
  mylon   = as.numeric(myCoord[i,3])           #long
  num      = as.numeric(myCoord[i,8])           #june

  if(161  %<% num %<% 169){temp=tmpin_8d_1}else{
    if(168 %<% num %<% 177){temp=tmpin_8d_2}else{
      if(176 %<% num %<% 185){temp=tmpin_8d_3}else{
        if(184 %<% num %<% 193){temp=tmpin_8d_4}else{temp=tmpin_8d_5}}}}

  longMin = as.numeric(mylon - 0.1)
  longMax = as.numeric(mylon + 0.1)
  latMin  = as.numeric(mylat - 0.1)
  latMax  = as.numeric(mylat + 0.1)

  temp_grid1 = lat[lat>latMin]
  temp_grid2 = temp_grid1[temp_grid1<latMax]
  temp_grid3 = lon[lon>longMin]
  temp_grid4 = temp_grid3[temp_grid3<longMax]

  latWindx=c()
  for(k in 1:length(temp_grid2)){
    t1 = which.min(abs(lat - temp_grid2[k]))
    latWindx = c(latWindx, t1)
  }
  lonWindx=c()
  for(k in 1:length(temp_grid4)){
    t1 = which.min(abs(lon - temp_grid4[k]))
    lonWindx = c(lonWindx, t1)
  }

  vars    = temp[latWindx,lonWindx]
  var_val = max(vars[!is.na(vars)])
  if(var_val == "-Inf"){var_val="NA"}

  out_8d[i,7]= var_val
  print(paste("KD490:", i))
  flush.console()
}
colnames(out_8d)[7] = "KD490_8d"
print(out_8d)


dt <- out %>% as_tibble() %>%
  left_join(out_8d %>% as_tibble()) %>%
  dplyr::select(-lat, -long)

library(RPostgreSQL)  # loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_satellite_data"), value=dt,overwrite=TRUE,row.names=FALSE)
