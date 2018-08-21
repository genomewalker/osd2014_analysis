library(rgdal)
library(maptools)
library(rgeos)
# From http://www.imachordata.com/meow-its-marine-ecoregions-in-r-2/

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)

osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)

ogrInfo("~/Downloads/MEOW/", "meow_ecos")
regions <- readOGR("~/Downloads/MEOW/", "meow_ecos")
provinces <- unionSpatialPolygons(regions, regions$PROVINCE)
realms <- unionSpatialPolygons(regions, regions$REALM)

prov_data <- regions@data %>%
  group_by(PROVINCE) %>%
  summarise(PROV_CODE = PROV_CODE[1], REALM = REALM[1], RLM_CODE=RLM_CODE[1], Lat_Zone=Lat_Zone[1])

provinces <- SpatialPolygonsDataFrame(provinces,
                                      data=data.frame(
                                        inner_join(data.frame(PROVINCE=names(provinces)),
                                             prov_data),
                                        row.names=row.names(provinces)))

regions@data$id = rownames(regions@data)
provinces@data$id = rownames(provinces@data)

regions.points = fortify(regions, ECOREGION="id")
provinces.points = fortify(provinces, PROVINCES="id")

regions.df = inner_join(regions.points, regions@data, by="id")
provinces.df = inner_join(provinces.points, provinces@data, by="id")


provinces.df_filt  <- bind_rows(provinces.df %>%
                                     filter(PROVINCE == "Mediterranean Sea") %>%
                                     mutate(region_new ="Mediterranean"),
                                provinces.df %>%
                                     filter(PROVINCE == "Cold Temperate Northwest Atlantic" | PROVINCE == "Lusitanian" | PROVINCE == "Northern European Seas" | PROVINCE == "Warm Temperate Northwest Atlantic") %>%
                                     mutate(region_new = "Eur North Atl"),
                                provinces.df %>%
                                     filter(grepl("Atlantic", PROVINCE) & grepl("North", PROVINCE)) %>%
                                     mutate(region_new = "USA North Atl"))






osd2014_cdata_old <- tbl(my_db, "osd2014_cdata_old") %>%
  collect(n = Inf)

osd2014_cdata_coord_paper <- osd2014_cdata %>%
  select(label, start_lat, start_lon) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  mutate(class = "a_paper")

osd2014_cdata_coord_all <- osd2014_cdata_old %>%
  select(label, start_lat, start_lon) %>%
  filter(!(label %in% osd2014_amp_mg_intersect$label)) %>%
  mutate(class = "b_all")




global <- map_data("world")

p1 <- ggplot() +
  geom_polygon(data = global, aes(x=long, y = lat, group = group)) +
  geom_point(data = bind_rows( osd2014_cdata_coord_all, osd2014_cdata_coord_paper), aes(x = start_lon, y = start_lat, shape = class, fill = class),
             size = 2.2, color = "grey") +
  xlab("") +
  ylab("") +
  coord_equal(xlim = c(-160, 170))  +
  scale_fill_manual(values = c("#C02942", "#2A67A0"), labels = c("NPL022 paper", "NPL022 all")) +
  theme_map(legend.position = "bottom") +
  #scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(21,22),  labels = c("NPL022 paper", "NPL022 all"))



#colors <- randomcoloR::distinctColorPalette(k = 6, runTsne = FALSE)
colors <- c("#B99ED0", "#ABE172", "#7AD4CB", "#DB8763", "#C059D1", "#D2D0B7")
p2 <- ggplot() +
  geom_polygon(data = provinces.df_filt, aes(long,lat,group=group, fill = PROVINCE), alpha = 1) +
  geom_path(data = provinces.df_filt, aes(long,lat,group=group), color="black") +
  geom_polygon(data = global, aes(x=long, y = lat, group = group)) +
  geom_point(data = bind_rows(osd2014_cdata_coord_paper), aes(x = start_lon, y = start_lat), fill = "#C02942", size = 2.2, color = "grey", shape = 21) +
  coord_equal(xlim = c(-100, 50), ylim = c(0, 70)) +
  # coord_equal() +
  theme_map(legend.position = "bottom") +
  scale_fill_manual(values = colors)

ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(0.65, 0.35))


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

save.image("osd2014_ancillary_data/data/osd2014_meow.Rda")
