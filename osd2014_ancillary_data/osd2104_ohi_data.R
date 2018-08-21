#================================================================================
#STEP 1: unzip files
#================================================================================
## set dir
setwd("~/Downloads/HALPERN")
## list all files for 2008

files_rescaled_2013 <- list.files(path = "2013_scaled/tif")
print(files_rescaled_2013)


#================================================================================
#STEP 2
#================================================================================
library(proj4)
library(tidyverse)
library(raster)
my_db<- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect") %>%
  collect(n = Inf)
osd2014_metadata <- read_tsv("~/ownCloud/OSD_paper/OSD-GC/data/osd2014_metadata_18-01-2017.tsv",
                             col_names = TRUE, trim_ws = TRUE,
                             col_types = list(mrgid = col_character(), biome_id = col_character(),
                                              feature_id = col_character(),material_id = col_character()))

osd2014_metadata_mg <- osd2014_metadata %>% filter(label %in% osd2014_amp_mg_intersect$label)
sc <- cbind(long=osd2014_metadata_mg$start_lon, lat=osd2014_metadata_mg$start_lat)                                #(long,lat)


r1 <- raster("~/Downloads/HALPERN/2013_scaled/tif_wgs84/global_cumul_impact_2013_all_layers_wgs84.tif")
r <- raster("~/Downloads/HALPERN/2013_scaled/tif/global_cumul_impact_2013_all_layers.tif")
#change WGS84 to Mollweide
oldproj <- "+proj=longlat +datum=WGS84 +no_defs"                                #WGS84
newproj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"   #Mollweide

wgs.84    <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# ESRI:54009 world mollweide projection, units = meters
# see http://www.spatialreference.org/ref/esri/54009/
mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


sp.points <- SpatialPoints(sc, proj4string=CRS(oldproj))
sc <- spTransform(sp.points,CRS(newproj))

#================================================================================
#STEP 3 : values for 2008
#================================================================================
require(rgdal)
library(raster)

#sink("GlobalMinMaxForVar.txt")
GlobalMinMaxForVar13 = data.frame(id=1:length(files_rescaled_2013), VarNames=gsub(".tif", "", files_rescaled_2013), minValue=rep(0,length(files_rescaled_2013)), maxValue=rep(0,length(files_rescaled_2013)))
head(GlobalMinMaxForVar13)

#OUTPUT FILE
out_rescaled_2013 = data.frame(osd_id=osd2014_metadata_mg$label, long=osd2014_metadata_mg$start_lon, lat=osd2014_metadata_mg$start_lat)
m=ncol(out_rescaled_2013)

#LOOP TO READ TIFF ONE BY ONE
for(i in 1:length(files_rescaled_2013)) {

  varName = gsub(".tif", "", files_rescaled_2013[i])
  print(varName); flush.console()

  ## Reading in a zip data file without unzipping it
  #r <- raster(list.files$Name[2])
  f <- paste("~/Downloads/HALPERN/2013_scaled/tif/", files_rescaled_2013[i], sep = "")
  r <- raster(f)
  print(r); flush.console()

  #Get min and max cell values from raster
  GlobalMinMaxForVar13$minValue[i] = cellStats(r, min)
  GlobalMinMaxForVar13$maxValue[i] = round(cellStats(r, max),3)

  #LOOP TO READ TIFF ONE BY ONE
  for(k in 1: nrow(sc@coords)){

    # if(k==16){
    # out_rescaled_2013[k, m+1] =  "NA"; colnames(out_rescaled_2013)[m+1] = (paste(varName, "50kms_max", sep="_"))
    # out_rescaled_2013[k, m+2] =  "NA"; colnames(out_rescaled_2013)[m+2] = (paste(varName, "50kms_min", sep="_"))
    # out_rescaled_2013[k, m+3] =  "NA"; colnames(out_rescaled_2013)[m+3] = (paste(varName, "10kms_max", sep="_"))
    # out_rescaled_2013[k, m+4] =  "NA"; colnames(out_rescaled_2013)[m+4] = (paste(varName, "10kms_min", sep="_"))
    #} else{
    ## MEHOD:1
    ## Extract values from Raster objects
    #equivalent to 0.5 degree or 50 kms approx radius
    e0 = extract(r, sc[k,], buffer = 100000, weights = TRUE, small = TRUE) %>% unlist

    e1 = extract(r, sc[k,], buffer = 50000, weights = TRUE, small = TRUE) %>% unlist

    #equivalent to 0.1 degree or 10 kms approx
    e2 = extract(r, sc[k,], buffer = 10000, weights = TRUE, small = TRUE) %>% unlist

    #equivalent to 0.1 degree or 10 kms approx
    e3 = extract(r, sc[k,], buffer = 5000, weights = TRUE, small = TRUE) %>% unlist

    e4 = extract(r, sc[k,], buffer = 1000, weights = TRUE, small = TRUE) %>% unlist

    out_rescaled_2013[k, m+1] =  round(max(e1, na.rm=TRUE),3) ; colnames(out_rescaled_2013)[m+1] = (paste(varName, "2013_50kms_max", sep="_"))
    out_rescaled_2013[k, m+2] =  round(median(e1, na.rm=TRUE),3)   ; colnames(out_rescaled_2013)[m+2] = (paste(varName, "2013_50kms_median", sep="_"))

    out_rescaled_2013[k, m+3] =  round(max(e2, na.rm=TRUE),3) ; colnames(out_rescaled_2013)[m+3] = (paste(varName, "2013_10kms_max", sep="_"))
    out_rescaled_2013[k, m+4] =  round(median(e2, na.rm=TRUE),3)   ; colnames(out_rescaled_2013)[m+4] = (paste(varName, "2013_10kms_median", sep="_"))

    out_rescaled_2013[k, m+5] =  round(max(e3, na.rm=TRUE),3) ; colnames(out_rescaled_2013)[m+5] = (paste(varName, "2013_5kms_max", sep="_"))
    out_rescaled_2013[k, m+6] =  round(median(e3, na.rm=TRUE),3)   ; colnames(out_rescaled_2013)[m+6] = (paste(varName, "2013_5kms_median", sep="_"))

    out_rescaled_2013[k, m+7] =  round(max(e4, na.rm=TRUE),3) ; colnames(out_rescaled_2013)[m+7] = (paste(varName, "2013_1kms_max", sep="_"))
    out_rescaled_2013[k, m+8] =  round(median(e4, na.rm=TRUE),3)   ; colnames(out_rescaled_2013)[m+8] = (paste(varName, "2013_1kms_median", sep="_"))

    out_rescaled_2013[k, m+9] =  round(max(e0, na.rm=TRUE),3) ; colnames(out_rescaled_2013)[m+9] = (paste(varName, "2013_100kms_max", sep="_"))
    out_rescaled_2013[k, m+10] =  round(median(e0, na.rm=TRUE),3)   ; colnames(out_rescaled_2013)[m+10] = (paste(varName, "2013_100kms_median", sep="_"))


    if(k %in% c(25,50,75,100,125,150,175,200)){cat(k, " | ")}
    flush.console()
    #}
  }
  cat("Next!!! \n")
  m=ncol(out_rescaled_2013)
}


out_rescaled_2013_median_long <- out_rescaled_2013 %>%
  dplyr::select(osd_id, contains("median")) %>%
  dplyr::rename(label = osd_id) %>%
  gather(ohi_variable, median, -label) %>% tbl_df %>%
  mutate(ohi_variable = gsub("_2013_all_layers", "", ohi_variable)) %>%
  mutate(ohi_variable = gsub("_2013_minus_", "_diff_", ohi_variable)) %>%
  separate("ohi_variable", c("ohi_variable", "buffer"), sep = "_2013_", remove = TRUE) %>%
  mutate(buffer = gsub("s_median", "", buffer))

out_rescaled_2013_max_long <- out_rescaled_2013 %>%
  dplyr::select(osd_id, contains("max")) %>%
  dplyr::rename(label = osd_id) %>%
  gather(ohi_variable, max, -label) %>% tbl_df %>%
  mutate(ohi_variable = gsub("_2013_all_layers", "", ohi_variable)) %>%
  mutate(ohi_variable = gsub("_2013_minus_", "_diff_", ohi_variable)) %>%
  separate("ohi_variable", c("ohi_variable", "buffer"), sep = "_2013_", remove = TRUE) %>%
  mutate(buffer = gsub("s_max", "", buffer))


library(RPostgreSQL)  # loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_halpern_scaled_median"), value=out_rescaled_2013_median_long,overwrite=TRUE,row.names=FALSE)



library(naturalsort)
out_rescaled_2013_median_long$ohi_variable <- factor(out_rescaled_2013_median_long$ohi_variable, levels = unique(out_rescaled_2013_median_long$ohi_variable) %>% naturalsort())

for (i in 1:13){
  fname <- paste("~/Downloads/halpern_median_", i, "_042017.pdf", sep = "")
  pdf(fname, width = 14, height = 3)
  p <- ggplot(out_rescaled_2013_median_long %>% filter(ohi_variable == "global_cumul_impact"), aes(median)) +
    geom_density(fill = "#333333", alpha = 0.8) +
    #ggforce::facet_wrap_paginate(~ohi_variable, scales = "free", ncol = 5, nrow = 1, page = i) +
    facet_grid(buffer~ohi_variable, scales = "free_y") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 6),
          axis.text = element_text(size = 8))
  print(p)
  dev.off()
}

out_rescaled_2013_median_long %>%
  filter(is.na(median)) %>%
  group_by(ohi_variable) %>%
  count() %>%
  separate("ohi_variable", c("variable", "buffer"), sep = "_2013_", remove = TRUE) %>%
  mutate(buffer = gsub("s_median", "", buffer)) %>%
  complete(variable, buffer, fill = list(n = 0)) %>%
  ggplot(aes(buffer, n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.2, size = 3) +
  facet_wrap(~variable, scales = "free_x") +
  scale_x_discrete(limits = c("1km", "5km", "10km", "50km", "100km")) +
  theme_bw() +
  xlab("Extraction buffer") +
  ylab("# of NAs")


library(maptools)
data(wrld_simpl)
wrldMoll <- spTransform(wrld_simpl, CRS(newproj))
for(i in 1:length(files_rescaled_2013)) {

  varName = gsub(".tif", "", files_rescaled_2013[i])
  print(varName); flush.console()

  ## Reading in a zip data file without unzipping it
  #r <- raster(list.files$Name[2])
  f <- paste("~/Downloads/HALPERN/2013_scaled/tif/", files_rescaled_2013[i], sep = "")
  r <- raster(f)

  fname <- paste("~/Downloads/halpern_", varName,".pdf", sep = "")
  pdf(fname, width = 7, height = 5)
  colors <- rev(colorRampPalette(c(brewer.pal(11, "Spectral"), "#FFFFFF"))(100))
  plot(r, col = colors, axes = F, box = F)
  plot(wrldMoll, col = "#E5E6E6", border = "#E5E6E6", add = TRUE, box = FALSE)
  points(sc, cex = 0.5)
  dev.off()
}


out_rescaled_2013_median_long_5km <- out_rescaled_2013_median_long %>%
  mutate(ohi_variable = gsub("_2013_all_layers", "", ohi_variable)) %>%
  mutate(ohi_variable = gsub("_2013_minus_2008", "diff", ohi_variable))


out_rescaled_2013_median_long_5km <- out_rescaled_2013_median_long_5km %>%
  separate("ohi_variable", c("variable", "buffer"), sep = "_2013_", remove = TRUE) %>%
  mutate(buffer = gsub("s_median", "", buffer))

out_rescaled_2013_median_long_5km$label <- factor(out_rescaled_2013_median_long_5km$label, levels = st_100_order_terrestrial)
out_rescaled_2013_median_long_5km$buffer <- factor(out_rescaled_2013_median_long_5km$buffer, levels = c("1km", "5km", "10km", "50km", "100km"))

for (i in 1:13){
  fname <- paste("~/Downloads/halpern_median_", i, "_by_sample.pdf", sep = "")
  pdf(fname, width = 14, height = 3)
  p <- ggplot(out_rescaled_2013_median_long, aes(label, median)) +
    geom_bar(stat = "identity") +
    ggforce::facet_grid_paginate(variable~buffer, scales = "free", ncol = 5, nrow = 1, page = i) +
    scale_x_discrete(limits = st_100_order_terrestrial) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("") +
    ylab("Average")
  print(p)
  dev.off()
}


out_rescaled_2013_median_long_lat_long <- out_rescaled_2013_median_long %>%
  filter(is.na(median)) %>%
  separate("ohi_variable", c("variable", "buffer"), sep = "_2013_", remove = TRUE) %>%
  mutate(buffer = gsub("s_median", "", buffer)) %>%
  filter(buffer == "10km") %>%
  left_join(osd2014_metadata)


wmap<-map_data("world")

pdf("~/Downloads/nas_map_10km.pdf", width = 21, height = 10)
p <- ggplot(wmap, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), fill = "#E5E6E6") +
  geom_path(aes(group = group), colour = "#E5E6E6") +
  geom_point(data = out_rescaled_2013_median_long_lat_long, aes(x = start_lon, y = start_lat),
             alpha = 0.8) +
  xlab("") +
  ylab("") +
  coord_equal(ratio = 1)  +
  facet_wrap(~variable) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")
print(p)
dev.off()






# filt <- osd2014_16s_prevalence %>%
#   filter(total_counts >= 10, prev > 2) %>% .$otu_name %>% droplevels()
#
# sites_rm <- c("OSD63_2014-06-20_0m_NPL022",
#               "OSD167_2014-06-21_0.1m_NPL022",
#               "OSD128_2014-06-21_1m_NPL022",
#               "OSD80_2014-06-21_2m_NPL022")
#
# sites_rm <- sample_names(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg)[!(sample_names(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg) %in% sites_rm)]
#
#
# osd2014_16s_otuXsample_physeq_filt_prev_beta_filt <- prune_taxa(filt %>% as.vector(), osd2014_16s_otuXsample_physeq_filt_prev_beta_mg)
# osd2014_16s_otuXsample_physeq_filt_prev_beta_filt <- prune_samples(sites_rm,osd2014_16s_otuXsample_physeq_filt_prev_beta_filt)
#
# osd2014_16s_otuXsample_physeq_filt_prev_beta_filt_css <- cssTrans(osd2014_16s_otuXsample_physeq_filt_prev_beta_filt, norm = T, log = F)
#
# osd2014_16s_otuXsample_physeq_filt_prev_beta_filt_css_m <- t(as(otu_table(osd2014_16s_otuXsample_physeq_filt_prev_beta_filt_css), "matrix"))
#
# plot(hclust(dist(osd2014_16s_otuXsample_physeq_filt_prev_beta_filt_css_m), "ward.D"))
#
#
# pca <- prcomp((osd2014_16s_otuXsample_physeq_filt_prev_beta_filt_css_m), center=F)
#
# osd2014_16s_otuXsample_physeq_filt_prev_beta_filt_css_m_hel <- decostand(osd2014_16s_otuXsample_physeq_filt_prev_beta_filt_css_m, "total")
#
# osd2014_16s_otuXsample_physeq_filt_prev_beta_filt_css_m_hel <- decostand(osd2014_16s_otuXsample_physeq_filt_prev_beta_filt_css_m_hel, "hellinger")
#
#
# mdsplot <- metaMDS(osd2014_16s_otuXsample_physeq_filt_prev_beta_filt_css_m_hel, distance = "bray")
#
#
# halpern_impact_all <- out_rescaled_2013_median_long %>%
#   mutate(ohi_variable = gsub("_2013_all_layers", "", ohi_variable)) %>%
#   mutate(ohi_variable = gsub("_2013_minus_2008", "diff", ohi_variable)) %>%
#   filter(!is.na(median)) %>%
#   separate("ohi_variable", c("variable", "buffer"), sep = "_2013_", remove = TRUE) %>%
#   mutate(buffer = gsub("s_median", "", buffer)) %>%
#   filter(buffer == "10km")
#
# test <- melt(scores(pca, display = "sites")) %>%
#   dplyr::rename(label = Var1) %>%
#   filter(Var2 == "NMDS1") %>%
#   left_join(halpern_impact_all) %>%
#   dplyr::select(label, value, variable, median) %>%
#   group_by(variable) %>%
#   do(fit = Hmisc::rcorr(.$value, .$median))
#
#
#
# test %>% broom::tidy(fit) %>%
#   filter(p.value < 0.01)
#
#

library(tidyverse)
library(phyloseq)
my_db<- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect") %>%
  collect(n = Inf)
osd2014_metadata <- read_tsv("~/ownCloud/OSD_paper/OSD-GC/data/osd2014_metadata_18-01-2017.tsv",
                             col_names = TRUE, trim_ws = TRUE,
                             col_types = list(mrgid = col_character(), biome_id = col_character(),
                                              feature_id = col_character(),material_id = col_character()))

osd2014_metadata_mg <- osd2014_metadata %>% filter(label %in% osd2014_amp_mg_intersect$label)
out_rescaled_2013_median_long <- tbl(my_db, "osd2014_halpern_scaled_median") %>%
  collect(n = Inf)

out_rescaled_2013_median_long$buffer <- factor(out_rescaled_2013_median_long$buffer, levels = c("1km", "5km", "10km", "50km", "100km"))

p <- ggplot(out_rescaled_2013_median_long %>% filter(ohi_variable == "global_cumul_impact"), aes(median)) +
  geom_density(fill = "#333333", alpha = 0.8) +
  #ggforce::facet_wrap_paginate(~ohi_variable, scales = "free", ncol = 5, nrow = 1, page = i) +
  facet_grid(buffer~ohi_variable, scales = "free_y") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 8))
print(p)

# Very low impact: <1.4
# Low impact: 1.4-4.95
# Medium impact: 4.95-8.47
# Medium high impact: 8.47-12
# High impact: 12-15.52
# Very high impact: >15.52

# Classify samples in impacted and non-impacted
# halpern_impact <- out_rescaled_2013_median_long %>%
#   filter(buffer == "1km", ohi_variable == "global_cumul_impact") %>%
#   mutate(class = ifelse(median <= 1.4, "non-impacted",
#                          ifelse( (median > 1.4 & median <= 4.8), "low_impact",
#                                  "impacted")))


halpern_impact <- out_rescaled_2013_median_long %>%
  filter(buffer == "1km", ohi_variable == "global_cumul_impact") %>%
  mutate(class = ifelse(median <= 0, "non_impacted",
                        ifelse((median > 0.01 & median <= 4.8), "low_impacted",
                               "impacted"))) %>%
  filter(median != 0, class != "non_impacted")  %>% filter(class == "low_impacted" | (class == "impacted" & median > 6)) %>%
  left_join(osd2014_metadata_mg) %>%
  filter(dist_coast_m <= 2000)


# halpern_impact <- out_rescaled_2013_median_long %>%
#   filter(buffer == "1km", ohi_variable == "global_cumul_impact") %>%
#   mutate(class = ifelse(median <= 1.4, "non_impacted",
#                         ifelse( (median > 1.4 & median <= 4.95), "low_impacted",
#                                 "impacted")))
#
#
# halpern_impact <- out_rescaled_2013_median_long %>%
#   filter(buffer == "1km", ohi_variable == "global_cumul_impact") %>%
#   mutate(class = ifelse(median <= 4.95, "non_impacted", "impacted"))

halpern_impact$class %>% table


halpern_impact_class_min <- halpern_impact %>%
  group_by(class) %>%
  count %>%
  ungroup() %>%
  dplyr::slice(which.min(n)) %>%
  .$n

sites_rand_non <- halpern_impact %>%
  filter(class == "non_impacted") %>%
  .$label %>%
  sample(.,halpern_impact_class_min)

sites_rand_low <- halpern_impact %>%
  filter(class == "low_impacted") %>%
  arrange(median) %>%
  head(n = halpern_impact_class_min) %>%
  .$label

sites_rand_imp <- halpern_impact %>%
  filter(class == "impacted") %>%
  arrange(desc(median)) %>%
  head(n = halpern_impact_class_min) %>%
  .$label


sites <- unique(c( sites_rand_low, sites_rand_imp))


# check stressors

sites_rand_imp <- out_rescaled_2013_median_long %>%
  filter(buffer == "1km", label %in% sites_rand_imp, ohi_variable != "sst", ohi_variable != "uv") %>%
  left_join(osd2014_metadata)

sites_rand_imp$label <- factor(sites_rand_imp$label, levels = sites_rand_imp %>% filter(ohi_variable == "global_cumul_impact") %>% arrange(median) %>% .$label %>% unique)

ggplot(sites_rand_imp, aes(label, median)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ohi_variable, scales = "free") +
  theme_light() +
  theme()

sites_rand_low<- out_rescaled_2013_median_long %>%
  filter(buffer == "1km", label %in% sites_rand_low, ohi_variable != "sst", ohi_variable != "uv") %>%
  left_join(osd2014_metadata)

sites_rand_low$label <- factor(sites_rand_low$label, levels = sites_rand_low %>% filter(ohi_variable == "global_cumul_impact") %>% arrange(median) %>% .$label %>% unique)

ggplot(sites_rand_low, aes(label, median)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ohi_variable, scales = "free")

# plot sites in a map


load(file = "~/ownCloud/OSD2014_data_4_Pier/phyloseq_exports/osd2014_16S_export.Rda", verbose = TRUE)


#sites <- halpern_impact$label
#sites <- out_rescaled_2013_median_long %>% filter(!(label %in% sites)) %>% .$label %>% unique

osd2014_16s_otuXsample_physeq_filt_prev_beta_mg <- prune_samples(sites, osd2014_16s_otuXsample_physeq_filt_prev_beta_mg)

simpletrim = function(physeq, minobs) {
  ## remove the rare otus

  Ji = nsamples(physeq)
  # Force orientation to be sample-by-OTU
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  # `prevalence` is the fraction of total samples in which an OTU is observed
  # at least `minobs` times.
  prevalence = apply(as(otu_table(physeq), "matrix"), 2, function(x, minobs) {
    return(sum(x > minobs))
  }, minobs)/(Ji)
  # Will only keep OTUs that appear in more than X% of samples and have total
  # reads greater than half the number of samples.
  ## keepOTUs = prevalence > 0.05 & taxa_sums(physeq) > (0.5 * Ji)
  keepOTUs = prevalence > 0.1
  return(prune_taxa(keepOTUs, physeq))
}
# Simple prune
osd2014_16s_otuXsample_physeq_filt_prev_beta_mg <- simpletrim(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg, 2)

osd2014_16s_otuXsample_physeq_filt_prev_beta_mg


library(metagenomeSeq)
cssTrans<-function(f.physeq.p = f.physeq.p, norm = norm, log = log){

  if (taxa_are_rows(f.physeq.p)) {
    f.physeq.p <- (f.physeq.p)
  }else{
    f.physeq.p <- t(f.physeq.p)
  }

  OTU <- as((otu_table(f.physeq.p)), "matrix")
  MGS <- newMRexperiment(
    counts = (OTU)
  )
  MGS <- cumNorm(MGS, p = cumNormStat(MGS))
  f.norm.p <- f.physeq.p
  otu_table(f.norm.p) <- otu_table((as.matrix(MRcounts(
    MGS,
    norm = norm,
    log = log,
    sl = median(unlist(normFactors(MGS)))
  ))), taxa_are_rows = T)
  return(f.norm.p)
}


#osd2014_16s_otuXsample_physeq_filt_prev_beta_mg <- cssTrans(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg, norm = T, log = T)

osd2014_16s_otuXsample_physeq_filt_prev_beta_mg <- transform_sample_counts(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg, function(x) x/sum(x))

otutable <- (as(otu_table(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg), "matrix"))
#otutable <- as(compositions::ilr(otutable), "matrix")
otutable <- vegan::decostand(otutable, "hellinger")

otutable <- base::as.data.frame(otutable)
names_otus <- names(otutable)
names_otus_short <- paste("otu", seq(1:length(names_otus)), sep = "_")
names(otutable) <- names_otus_short
otutable$label <- base::row.names(otutable)
#sites <- unique(c(sites_rand_non, sites_rand_low, sites_rand_imp))

otutable <- otutable[sites,]



otutable <- otutable  %>% left_join(halpern_impact %>% dplyr::select(label, class))
class_impact <- otutable %>% dplyr::select(label, class)
base::row.names(otutable) <- otutable$label
otutable$label <- NULL



library(randomForest)

#samp <- sample(nrow(otutable), 2/3 * nrow(otutable))
otutable$class <- as.factor(otutable$class)
#train <- as(otutable[samp, ], "data.frame")
#train$class <- as.factor(train$class)

#train <- train[, colSums(train != 0) > 0]

#test <-  as(otutable[-samp,], "data.frame")
#test$class <- as.factor(test$class)

#train_class <- class_impact %>% filter(label %in% rownames(train))
#train_class <- train_class[match(rownames(train), train_class$label),]

#test_class <- class_impact %>% filter(label %in% rownames(test))
#test_class <- test_class[match(rownames(test), test_class$label),]

library(doParallel)
library(caret)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
# Random Search
seed <- 1234
set.seed(seed)

#length is = (n_repeats*nresampling)+1
seeds <- vector(mode = "list", length = 50)

#(3 is the number of tuning parameter, mtry for rf, here equal to ncol(iris)-2)
for(i in 1:50) seeds[[i]]<- sample.int(n=1000, 15)

#for the last model
seeds[[51]]<-sample.int(1000, 1)

control <- trainControl(method="CV", number=10, repeats=5, search="random", seeds = seeds)
rf_random <- train(y = otutable$class, x=otutable[1:dim(otutable)[2]-1], method="rf", metric="Accuracy", tuneLength=15,
                   trControl=control, verboseIter = TRUE, savePredictions = TRUE,
                   prox=TRUE, allowParallel=TRUE, num.trees=200)
stopCluster(cl)
print(rf_random)
plot(rf_random)

#mtry = 329

n<-0.25
nt <- 100
repeat{
  halpern_classify <- randomForest::randomForest(y = otutable$class,
                                                 x = otutable[1:dim(otutable)[2]-1], ntree=nt, importance=TRUE, do.trace=nt, keep.inbag = TRUE, mtry = rf_random$bestTune$mtry)
  err <- halpern_classify$err.rate[,1][nt]
  diff <- abs(halpern_classify$err.rate[,2][nt] - halpern_classify$err.rate[,3][nt])
  if(err < n && diff < 0.06){
    break
  }
}

#save(names_otus, rf_random, names_otus_short, halpern_classify, halpern_impact, otutable, osd2014_16s_otuXsample_physeq_filt_prev_beta_mg, file = "~/ownCloud/OSD_paper/HALPERN/halpern_classify_RF_model.Rda")
plot(halpern_classify)

print(halpern_classify)

imp <- randomForest::importance(halpern_classify)
imp <- data.frame(predictors = rownames(imp), imp)

imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 10 predictors
imp.sort <- imp.sort[1:5, ]


# ggplot
ggplot(imp.sort, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying into low or mid/hight impacted OSD samples") +
  theme_light() +
  ylab("Total decrease in node impurities (Mean decrease Gini)")

library(ggstance)
gg <- ggplot(imp.sort, aes(MeanDecreaseGini, predictors, label = round(MeanDecreaseGini, 1))) +
  geom_segment(aes(x = 0, y = predictors, xend = MeanDecreaseGini, yend = predictors), color = "black", alpha = 1) +
  geom_point(size = 1, color = "black") +
  geom_rect(aes(ymin = as.numeric(predictors) - 0.5 , ymax = as.numeric(predictors) + 0.5, xmin=-0.01,
                xmax=1.85, color = NULL), alpha=.5) +
  scale_x_continuous(position = "top", expand = c(0, 0)) +
  xlab("Total decrease in node impurities (Mean decrease Gini)")  +
  ylab("OTU") +
  guides(color = FALSE, alpha = FALSE) +
  scale_fill_manual(values = c("#ECD078", "#D95B43", "#C02942", "#542437", "#53777A"), name = "BGC class") +
  scale_color_manual(values = c("#ECD078", "#D95B43", "#C02942", "#542437", "#53777A")) +
  theme_light() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))


predictors <- gsub("otu_", "", imp.sort$predictors) %>% as.numeric()
predictors_name <- names_otus[predictors]
tax_table(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg)[predictors_name,]



load(file = "~/ownCloud/OSD2014_data_4_Pier/phyloseq_exports/osd2014_16S_export.Rda", verbose = TRUE)

osd2014_16s_otuXsample_physeq_filt_prev_beta_mg <- transform_sample_counts(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg, function(x) x/sum(x))

osd2014_16s_otuXsample_physeq_filt_prev_beta_mg_filt <- prune_taxa( predictors_name, osd2014_16s_otuXsample_physeq_filt_prev_beta_mg)

plot_bar(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg_filt, x = "Sample", y = "Abundance", fill = "order")


tax_table(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg_filt) %>%
  as.data.frame()

f<-otu_table(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg_filt) %>%
  reshape2::melt() %>%
  tbl_df %>%
  dplyr::rename(label = Var2, otu_id = Var1) %>%
  left_join(out_rescaled_2013_median_long) %>%
  filter(ohi_variable == "global_cumul_impact", buffer == "1km")  %>%
  filter(label %in% sites, value > 0) %>%
  left_join(halpern_impact)

f$otu_id <- factor(f$otu_id, levels =  predictors_name)
l_ord <- halpern_impact %>% arrange(median) %>% .$label

f$class <- factor(f$class, levels = c("low_impacted", "impacted"))
g <- ggplot(f, aes(label, value, fill = class, group = 1)) +
  #geom_density(data = f %>% filter(class == "low_impacted"), aes(label, value, group = 1), stat = "identity", fill = "grey", alpha = 0.5, size = 0) +
  #geom_density(data = f %>% filter(class == "impacted"), aes(label, value, group = 1), stat = "identity", fill = "grey", alpha = 0.5, size = 0) +
  geom_bar(width = 1, stat = "identity", alpha = 0.8, color = "black", size = 0.1) +
  #geom_line(data = f %>% filter(class == "low_impacted"), aes(label, value, group = 1))+
  #geom_line(data = f %>% filter(class == "impacted"), aes(label, value, group = 1))+
  #geom_point(data = f, shape = 21, alpha = 0.7, aes(fill = class)) +
  facet_wrap( ~otu_id, scales = "free_y", ncol = 5) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(limits = l_ord) +
  scale_fill_manual(values = c("#4A4A4A","#C84359")) +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 9),
        legend.key.size = unit(0.25, "cm"))+
  ylab("OTU relative abundance") +
  xlab("")

ggplot(f, aes(median, value, fill = class, group = 1)) +
  #geom_density(data = f %>% filter(class == "low_impacted"), aes(label, value, group = 1), stat = "identity", fill = "grey", alpha = 0.5, size = 0) +
  #geom_density(data = f %>% filter(class == "impacted"), aes(label, value, group = 1), stat = "identity", fill = "grey", alpha = 0.5, size = 0) +
  #geom_bar(width = 1, stat = "identity") +
  #geom_line(data = f %>% filter(class == "low_impacted"), aes(label, value, group = 1))+
  #geom_line(data = f %>% filter(class == "impacted"), aes(label, value, group = 1))+
  geom_point(data = f, shape = 21, alpha = 0.8, aes(fill = class), size = 2) +
  facet_wrap( ~otu_id, scales = "free_y", ncol = 5) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#4A4A4A","#C84359")) +
  #scale_x_discrete(limits = l_ord) +
  theme_light() +
  theme(legend.position = "top",
        text = element_text(size = 9)) +
  ylab("OTU relative abundance") +
  xlab("Median global cumulative impact")


load(file = "~/ownCloud/OSD2014_data_4_Pier/phyloseq_exports/osd2014_16S_export.Rda", verbose = TRUE)


#sites <- halpern_impact$label
#sites <- out_rescaled_2013_median_long %>% filter(!(label %in% sites)) %>% .$label %>% unique

sites <- halpern_pred$label

osd2014_16s_otuXsample_physeq_filt_pred <- prune_samples(sites, osd2014_16s_otuXsample_physeq_filt_prev_beta_mg)

simpletrim = function(physeq, minobs) {
  ## remove the rare otus

  Ji = nsamples(physeq)
  # Force orientation to be sample-by-OTU
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  # `prevalence` is the fraction of total samples in which an OTU is observed
  # at least `minobs` times.
  prevalence = apply(as(otu_table(physeq), "matrix"), 2, function(x, minobs) {
    return(sum(x > minobs))
  }, minobs)/(Ji)
  # Will only keep OTUs that appear in more than X% of samples and have total
  # reads greater than half the number of samples.
  ## keepOTUs = prevalence > 0.05 & taxa_sums(physeq) > (0.5 * Ji)
  keepOTUs = prevalence > 0.1
  return(prune_taxa(keepOTUs, physeq))
}
# Simple prune
osd2014_16s_otuXsample_physeq_filt_pred <- simpletrim(osd2014_16s_otuXsample_physeq_filt_pred, 2)

osd2014_16s_otuXsample_physeq_filt_pred


library(metagenomeSeq)
cssTrans<-function(f.physeq.p = f.physeq.p, norm = norm, log = log){

  if (taxa_are_rows(f.physeq.p)) {
    f.physeq.p <- (f.physeq.p)
  }else{
    f.physeq.p <- t(f.physeq.p)
  }

  OTU <- as(t(otu_table(f.physeq.p)), "matrix")
  MGS <- newMRexperiment(
    counts = (OTU)
  )
  MGS <- cumNorm(MGS, p = cumNormStat(MGS))
  f.norm.p <- f.physeq.p
  otu_table(f.norm.p) <- otu_table((as.matrix(MRcounts(
    MGS,
    norm = norm,
    log = log,
    sl = median(unlist(normFactors(MGS)))
  ))), taxa_are_rows = T)
  return(f.norm.p)
}



#osd2014_16s_otuXsample_physeq_filt_prev_beta_mg <- cssTrans(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg, norm = T, log = T)

osd2014_16s_otuXsample_physeq_filt_prev_beta_mg <- transform_sample_counts(osd2014_16s_otuXsample_physeq_filt_prev_beta_mg, function(x) x/sum(x))

otutable_pred <- (as(otu_table(osd2014_16s_otuXsample_physeq_filt_pred), "matrix"))
names_otus_df <- data.frame(name = names_otus, short = names_otus_short)
pred_names_short <- names_otus_df %>%
  filter(name %in% intersect(colnames(otutable_pred),
                             names_otus_df %>% filter(short %in% colnames(otutable)) %>% .$name %>% as.character()
  ))

otutable_pred <- otutable_pred[,pred_names_short$name %>% as.character()]

colnames(otutable_pred) <- pred_names_short$short %>% as.character()

#otutable <- as(compositions::ilr(otutable), "matrix")
otutable_pred <- vegan::decostand(otutable_pred, "hellinger")
#sites <- unique(c(sites_rand_non, sites_rand_low, sites_rand_imp))







pred <- predict(halpern_classify, newdata = otutable_pred)

pred_impacted_metadata <- data.frame(class = pred) %>% mutate(label = names(pred)) %>% filter(class == "impacted") %>% left_join(osd2014_metadata)
pred_low_impacted_metadata <- data.frame(class = pred) %>% mutate(label = names(pred)) %>% filter(class != "impacted") %>% left_join(osd2014_metadata)

library(ggmap)

### Set a range

### Get a map

tmp <- pred_impacted_metadata %>% dplyr::select(label, start_lat,start_lon, dist_coast_iso3_code) %>%
  mutate(filename = paste("~/Downloads/", label, "_", dist_coast_iso3_code, "_impacted.pdf", sep = ""))


for (i in 1:dim(tmp)[1]) {

  map <- get_map(location = c(lon = tmp[i,]$start_lon, lat = tmp[i,]$start_lat), zoom = 13,
                 maptype = "satellite", source = "google")

  ### When you draw a figure, you limit lon and lat.
  foo <- ggmap(map) +
    geom_point(data = tmp[i,], aes(x = start_lon, y = start_lat), size = 3, shape = 21, fill = "red", alpha = 0.8)
  pdf(file = tmp[i,]$filename)
  print(foo)
  dev.off()
}


tmp <- pred_low_impacted_metadata %>% dplyr::select(label, start_lat,start_lon, dist_coast_iso3_code) %>%
  mutate(filename = paste("~/Downloads/", label, "_", dist_coast_iso3_code, "_low_impacted.pdf", sep = ""))


for (i in 1:dim(tmp)[1]) {

  map <- get_map(location = c(lon = tmp[i,]$start_lon, lat = tmp[i,]$start_lat), zoom = 13,
                 maptype = "satellite", source = "google")

  ### When you draw a figure, you limit lon and lat.
  foo <- ggmap(map) +
    geom_point(data = tmp[i,], aes(x = start_lon, y = start_lat), size = 3, shape = 21, fill = "red", alpha = 0.8)
  pdf(file = tmp[i,]$filename)
  print(foo)
  dev.off()
}


library(grid)
grid.newpage()
print(g, vp=viewport(width = 0.99, height = 0.12, x=0.5, y=0.75))
print(gg, vp=viewport(width = 0.4, height = 0.15, x=0.7, y=0.25))
