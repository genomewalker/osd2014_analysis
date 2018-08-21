interpolate_alpha <- function(X){
  load("osd2014_ancillary_data/data/osd2014_meow.Rda")

  # graphics
  library(ggplot2)
  library(RColorBrewer)

  # data manipulation
  library(tidyverse)
  # spatial analysis
  # library(raster)
  # library(geoR)
  # library(gstat)
  # library(fields)
  # library(maps)
   library(sp)


  ########################################
  #Function to generate grids
  ########################################
  generate_grid <- function(ticks = 100){
    require(maps)
    ## Generate Grid Locations
    max.lat <- 70    # most Northern, Continental U.S.
    min.lat <- 0    # most Southern, Continental U.S.
    max.lon <- 50   # most Eastern, Continental U.S.
    min.lon <- -100  # most Western, Continental U.S.

    grid <- expand.grid(seq(min.lon, max.lon, length=ticks),
                        seq(min.lat, max.lat, length=ticks))
    colnames(grid) <- c("lon", "lat")
    rownames(grid) <- 1:nrow(grid)
    return(grid)
  }

  #############################
  ####  Data  #################
  #############################
  # dat is a data frame with three columns: c("Longitude", "Latitude", "MDS1") where MDS1 is the variable you want
  # to interpolate
  dat <- X %>% inner_join(osd2014_cdata) %>% dplyr::select(start_lon, start_lat, mean)
  names(dat) <- c("Longitude", "Latitude", "MDS1")
  s <- dat

  t <- as.data.frame(generate_grid(250))  # grid on which to make predictions

  sp::coordinates(dat) = ~Longitude + Latitude
  sp::coordinates(t) <- ~lon + lat
  sp::gridded(t) <- TRUE

  #############################
  ####  Interpolation (idw)  #################
  #############################
  #change name of variable to interpolate and define the level of interpolation with idp
  # idw.variable <- idw(formula = MDS1 ~ 1, locations = dat, newdata = t, idp=3)
  #
  # idp = seq(from = 1, to = 5, by = 0.5)
  # nmax = seq(from = 2, to = 10, by = 1)
  #
  # IDW.out <- vector(mode = "list")
  # o<-1
  # for (id in idp){
  #   for (m in idp){
  #   for (n in nmax) {
  #     IDW.out1 <- vector(length = length(dat))
  #     for (i in 1:length(dat)) {
  #       IDW.out1[i] <- idw(MDS1 ~ 1, dat[-i,], dat[i,], idp=id, nmax = n, maxdist = m)$var1.pred
  #     }
  #     e <- sqrt( sum((IDW.out1 - dat$MDS1)^2, na.rm = TRUE) / length(dat$MDS1))
  #     IDW.out[[o]] <- data_frame(idp = id, nmax = n, maxdist = m, error = e)
  #     o <- o + 1
  #   }
  #   }
  # }
  #
  # bind_rows(IDW.out) %>% arrange(error)
  # Plot the differences


  idw.variable <- gstat::idw(formula = MDS1 ~ 1, locations = dat, newdata = t, idp=2, nmax = 4, maxdist = 5)
  #detach("gstat::", unload=TRUE)
  #library(raster, quietly=TRUE)

  l <- c("Tropical Northwestern Atlantic", "Warm Temperate Northwest Atlantic", "Cold Temperate Northwest Atlantic",
         "Lusitanian", "Mediterranean Sea", "Northern European Seas")
  ne = provinces[match(toupper(l),toupper(provinces$PROVINCE)),]

  idw.output = raster::as.data.frame(idw.variable)
  idw_raster <- raster::rasterFromXYZ(idw.output)
  xyz_hull <- rgeos::gConvexHull(ne)
  wgs.84    <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  # ESRI:54009 world mollweide projection, units = meters
  # see http://www.spatialreference.org/ref/esri/54009/
  mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  newproj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"   #Mollweide
  oldproj <- "+proj=longlat +datum=WGS84 +no_defs"
  xyz_hull <- SpatialPoints(xyz_hull, proj4string=CRS(oldproj))
  xyz_hull <- spTransform(xyz_hull,CRS(newproj))
  xyz_buff <- rgeos::gBuffer(xyz_hull, width = 25000) # arbitrary 1km buffer
  idw_raster_crop <- raster::mask(idw_raster, ne)

  idw.output <- raster::as.data.frame(idw_raster_crop, xy = TRUE)

  names(idw.output)[1:3] <- c("long", "lat", "var1.pred")

  #############################
  ####  Map  #################
  #############################
  mapWorld <- borders("world", colour="white", fill="#0a0c2a") # create a layer of borders
  #mp <- ggplot() +   mapWorld
  theme_map <- function(...) {
    theme_minimal() +
      theme(
        #text = element_text(family = "Ubuntu Regular", color = "#22211d"),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
        panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#f5f5f2", color = NA),
        panel.background = element_rect(fill = "#f5f5f2", color = NA),
        legend.background = element_rect(fill = "#f5f5f2", color = NA),
        panel.border = element_blank(),
        ...
      )
  }
  global <- map_data("world")
  ggplot(data = idw.output, aes(x = long, y = lat)) +
    #geom_polygon(data = provinces.df_filt, aes(long,lat,group=group, fill = PROVINCE), alpha = 1) +
    geom_polygon(data = provinces.df_filt, aes(long,lat,group=group), color="black", fill = "grey") +
    geom_raster(data = idw.output, aes(fill = var1.pred), na.rm = TRUE) +
    geom_polygon(data = global, aes(x=long, y = lat, group = group)) +
    geom_path(data = provinces.df_filt, aes(long,lat,group=group), color="black") +
    geom_polygon(data=map_data("world"), aes(x=long, y=lat, group=group),
                 alpha=1) +
    geom_point(data=s, aes(x=Longitude, y=Latitude), size=2, alpha=0.7, shape=21, fill = "grey90") +
    coord_equal(xlim = c(-100, 50), ylim = c(0, 70)) +
    #scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4", mid = "ivory2", midpoint=mean(sqrt(idw.output$var1.pred), na.rm = TRUE), trans = "sqrt", na.value="transparent") +
    viridis::scale_fill_viridis(option = "D", trans = "sqrt", na.value="transparent") +
    # coord_equal() +
    theme_map(legend.position = "none")


}

