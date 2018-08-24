library(tidyverse)
library(vegan)

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
osd2014_ancillary_data <- tbl(my_db, "osd2014_woa13_data") %>%
  collect(n = Inf)
st_100_order_terrestrial <- tbl(my_db, "osd2014_st_order_terrestrial") %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label)

st_100 <- tbl(my_db, "osd2014_st_100") %>%
  collect(n = Inf)

st_100_long <- st_100 %>%
  gather(variable, value, -label) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  mutate(class = ifelse(variable == "saline_water_gt_25km", "marine", "other")) %>%
  dplyr::select(label, class, value) %>%
  group_by(label, class) %>%
  summarise(prop=sum(value)) %>%
  filter(class == "marine") %>%
  select(-class) %>%
  rename(prop_marine = prop) %>%
  arrange(label)


# Get ancillary data ------------------------------------------------------

osd2014_mld_adata <- tbl(my_db, "osd2014_mld_data") %>%
  collect(n = Inf) %>%
  dplyr::select(osd_id, mld) %>%
  rename(label = osd_id)

osd2014_mld_adata %>% skimr::skim()

osd2014_phenology_adata <- tbl(my_db, "osd2014_phenology_data") %>%
  collect(n = Inf) %>%
  rename(label = osd_id)

osd2014_phenology_adata %>% skimr::skim()

osd2014_satellite_adata <- tbl(my_db, "osd2014_satellite_data") %>%
  collect(n = Inf) %>%
  rename(label = osd_id) %>%
  mutate(chlor_a_monthly = as.numeric(chlor_a_monthly),
         par_monthly = as.numeric(chlor_a_monthly),
         poc_monthly = as.numeric(poc_monthly),
         KD490_monthly = as.numeric(KD490_monthly),
         FLH_monthly = as.numeric(FLH_monthly),
         chlor_a_8d = as.numeric(chlor_a_8d),
         par_8d = as.numeric(par_8d),
         poc_8d = as.numeric(poc_8d),
         KD490_8d = as.numeric(KD490_8d))

osd2014_satellite_adata %>% skimr::skim()

osd2014_woa13_adata <- tbl(my_db, "osd2014_woa13_data") %>%
  collect(n = Inf) %>%
  dplyr::select(-long_woa13, -lat_woa13) %>%
  rename(label = osd_id)

osd2014_woa13_adata %>% skimr::skim()


osd2014_iron_adata <- tbl(my_db, "osd2014_iron_data") %>%
  collect(n = Inf) %>%
  select(osd_id, iron) %>%
  rename(label = osd_id)

osd2014_iron_adata %>% skimr::skim()

osd2014_halpern_adata <- tbl(my_db, "osd2014_halpern_scaled_median") %>%
  collect(n = Inf) %>%
  filter(buffer == "1km") %>%
  dplyr::select(label, ohi_variable, median) %>%
  tidyr::spread(key = ohi_variable, value = median, fill = NA)

osd2014_halpern_adata %>% skimr::skim()

l_pl_variable <- c(
  "start_lat",
  "start_lon",

  "temp", #WOD13
  "sal", #WOD13
  "silicate", #WOD13
  "oxy", #WOD13
  "phos", #WOD13

  "chlor_a_8d", #Satellite MODIS AQUA
  "KD490_8d", #Satellite MODIS AQUA
  "par_8d", #Satellite MODIS AQUA
  "poc_8d", #Satellite MODIS AQUA

  "pp_8d_0.5", #phenology SeaWiFs Sea-Viewing Wide Field-of-View Sensor

  "inorganic", #OHI
  "ocean_acidification", #OHI
  "ocean_pollution", #OHI
  "plumes_fert", #OHI
  "shipping", #OHI
  "prop_marine"
)


osd2014_adata <- osd2014_woa13_adata %>%
  left_join(osd2014_phenology_adata) %>%
  left_join(osd2014_satellite_adata) %>%
  left_join(osd2014_woa13_adata) %>%
  left_join(osd2014_iron_adata) %>%
  left_join(osd2014_halpern_adata %>%
              select(-global_cumul_impact_diff_2008, -sst, -uv)) %>%
  left_join(merge(data.frame(sample_data(osd2014_dada2_phyloseq_beta_vst)) %>% select(label, start_lat, start_lon), st_100_long)) %>%
  select(label, l_pl_variable)

# # Calculate PCA and get collinear -----------------------------------------
osd2014_asv_pca <- ordinate(osd2014_dada2_phyloseq_beta_vst, method = "RDA", distance = "bray", trymax = 1000, noshare = 0.1)

osd2014_asv_pca <- rda(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst))
#
# osd2014_asv_pca <- rda(varespec)
#
#
osd2014_adata_df <- as.data.frame(osd2014_adata %>% dplyr::select(-prop_marine)) %>% column_to_rownames("label")
#
osd2014_adata_df <- decostand(osd2014_adata_df, "standardize", na.rm = TRUE)
envf <- envfit(ord = scores(osd2014_asv_pca), env = osd2014_adata_df, permutations = 999, na.rm = TRUE, )
#
scores(envf, "vectors")
plot(osd2014_asv_pca)
plot(envf)
plot(envf, p.max = 0.05, col = "red")
#
#
l <- usdm::vifstep(osd2014_adata_df, th = 5)
#
non_coliner_variables <- l@results$Variables
#

osd2014_adata_noncol <- osd2014_adata %>% select(label, as.character(non_coliner_variables)) %>% filter(label %in% selected_samples)
# bioenv(osd2014_asv_bc, osd2014_adata_df[,as.vector(l@results$Variables)])# %>% dplyr::select(-temp ,-ocean_pollution ,-plumes_fert ,-inorganic, -oxy))
#
# arrowmat = vegan::scores(osd2014_asv_pca, display = "bp")
# arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# biplot(osd2014_asv_pca)
