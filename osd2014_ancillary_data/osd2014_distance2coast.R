library(rgdal)   # for readOGR(...); loads package sp as well
library(rgeos)   # for gDistance(...)

setwd("~/Downloads/ne_10m_coastline/")
# WGS84 long/lat
wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# ESRI:54009 world mollweide projection, units = meters
# see http://www.spatialreference.org/ref/esri/54009/
mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

osd2014_registry_coords <- read_tsv("~/ownCloud/OSD_paper/osd2014_registry_site_coords.txt", col_names = T, trim_ws = T)

osd2014_registry_coords <- osd2014_registry_coords %>% right_join(osd2014_metadata %>%
  dplyr::select(osd_id,label, start_lat, start_lon)) %>%
  mutate(lat_diff = abs(site_lat - start_lat), lon_diff =  abs(site_lon - start_lon))




df <- osd2014_metadata %>%
  filter(label %in% osd2014_rgc_prop_long$label) %>%
  dplyr::select(label, start_lon, start_lat)

sp.points <- SpatialPoints(df[c("start_lon","start_lat")], proj4string=CRS(wgs.84))

coast  <- readOGR(dsn=".",layer="ne_10m_coastline", p4s=wgs.84)
coast.moll <- spTransform(coast,CRS(mollweide))
point.moll <- spTransform(sp.points,CRS(mollweide))

result <- sapply(1:length(sp.points), function(i)gDistance(point.moll[i],coast.moll))
result/1000   # distance in km
#  [1]   0.2185196   5.7132447   0.5302977  28.3381043 243.5410571 169.8712255   0.4182755  57.1516195 266.0498881 360.6789699

df$dist_shore_km <- result/1000

plot(coast)
points(sp.points[test],pch=20,col="red")

df %>% left_join(osd2014_metadata %>%
                   filter(label %in% osd2014_rgc_prop_long$label) %>%
                   dplyr::select(label, dist_coast_m) %>%
                   mutate(dist_coast_km = dist_coast_m/1000))
