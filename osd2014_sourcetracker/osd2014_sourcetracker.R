library(tidyverse)
library(RPostgreSQL)

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)

osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)

osd2014_selected_meow_provinces <- tbl(my_db, "osd2014_selected_meow_provinces") %>%
  collect(n = Inf)

osd2014_file2label <- tbl(my_db, "osd2014_file2label") %>%
  collect(n = Inf)

#osd2014_rgc_prop <- read_tsv("data/OSD2014-exclusive_prop.txt", col_names = FALSE, trim_ws = TRUE) %>% dplyr::rename(label = X1, num_orf = X2, num_orf_excl = X3, prop = X4) %>%
osd2014_rgc_prop <- tbl(my_db, "osd95_eggnog_summary") %>%
  collect(n = Inf) %>%
  mutate(cog_description = ifelse(is.na(cog_description), "Not annotated", cog_description)) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  group_by(label, class) %>%
  summarise(total_abun = sum(N), n_orfs = sum(n_orfs), prop = sum(prop)) %>%
  ungroup()

st_100 <- tbl(my_db, "osd2014_st_100") %>%
  collect(n = Inf)

st_100_map <- tbl(my_db, "osd2014_st_100_map") %>%
  collect(n = Inf)

st_100_map <- st_100_map %>% mutate(Env = ifelse(Env == "river_water"| Env == "lake_water", "freshwater", Env))

st_100 <- st_100 %>%
  filter(label %in% osd2014_rgc_prop$label)

st_100_map$Env <- factor(st_100_map$Env, levels = c("Unknown", "agricultural_soil", "natural_soil", "forest_soil", "sediment",  "human_feces","freshwater",  "saline_water_le_25km", "saline_water_gt_25km"))

st_100_map_agg <- st_100_map %>% filter(SourceSink == "source") %>%
  group_by(Env) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(Env = fct_reorder(Env, (n)))

st_100_map_agg$Env <- factor(st_100_map_agg$Env,
                             labels =c("Forest soil", "Opean ocean", "Human feces", "Natural soil", "Sediment", "Coastal", "Agricultural soil", "Freshwater"))

agricultural_soil <- "#A28876"
forest_soil <- "#C6DFAC"
human_feces <- "#A7414C"
natural_soil <- "#DFC856"
saline_water_gt_25km <- "#4C6687"
saline_water_le_25km <- "#77AEB3"
sediment <- "#443649"
Unknown <- "#839496"
freshwater <- "#638877"

colors1 <- c(forest_soil, saline_water_gt_25km, human_feces, natural_soil, sediment, saline_water_le_25km, agricultural_soil, freshwater)

emp_sources <- ggplot(st_100_map_agg, aes(Env, n, fill = Env)) +
  geom_col(alpha = 0.9, size = 0.1, color = "black", width = 0.8) +
  ylab("Number of samples") +
  xlab("EMP biome") +
  ggpubr::rotate() +
  theme_bw() +
  scale_fill_manual(values = colors1) +
  theme(legend.position = "none")
ggsave("osd2014_sourcetracker/figures/osd2014_EMP_source_samples.pdf", width = 7.9, height = 4)


st_100 <- st_100 %>% mutate(freshwater = river_water + lake_water) %>%
  dplyr::select(-river_water, -lake_water)
st_100_mat <- as.data.frame(st_100)
rownames(st_100_mat) <- st_100_mat$label

st_100_long <- st_100 %>% gather(variable, value, -label)


meow_provinces <- c("Southern New Zealand", "Northern New Zealand", "Warm Temperate Northeast Pacific", "Warm Temperate Northwest Pacific",
"West and South Indian Shelf", "Sunda Shelf","Southeast Polynesia","Hawaii",
"Benguela", "Warm Temperate Southwestern Atlantic", "Tropical Northwestern Atlantic", "Warm Temperate Northwest Atlantic",
"Cold Temperate Northwest Atlantic","Lusitanian", "Mediterranean Sea", "Red Sea and Gulf of Aden", "Black Sea", "Northern European Seas", "Arctic")

st_100_order_terrestrial <- st_100_long %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  mutate(class = ifelse(variable == "saline_water_gt_25km", "marine", "terrestrial")) %>%
  dplyr::select(label, class, value) %>%
  group_by(label, class) %>%
  dplyr::summarise(N=sum(value)) %>%
  filter(class == "marine") %>%
  left_join(osd2014_cdata) %>%
  arrange(match(meow_province, meow_provinces), start_lat, N) %>% select(label, start_lat, N, meow_province) %>%
  .$label

st_100_order_terrestrial <-  tbl(my_db, "osd2014_st_order_coastal") %>%
  collect(n = Inf)

#st_100_long$label <- factor(st_100_long$label, levels = st_100_order_terrestrial)

st_100_long$variable <- factor(st_100_long$variable, levels = c("Unknown", "agricultural_soil", "natural_soil", "forest_soil", "sediment",  "human_feces","freshwater",  "saline_water_le_25km", "saline_water_gt_25km"))

colors <- c(Unknown, agricultural_soil, natural_soil, forest_soil, sediment, human_feces, freshwater, saline_water_le_25km, saline_water_gt_25km)
#colors <- c(base0, red, green, orange, violet, magenta, cyan,  yellow, blue)
#colors <- c("#839496", "#b58900", "#268bd2", "#dc322f", "#2aa198", "#859900", "#cb4b16", "#6c71c4", "#F5A10B")

st_100_barplot <- ggplot(st_100_long , aes(label, value, fill = factor(variable), group = variable)) +
  geom_bar(position = "stack", stat = "identity", width = 1, alpha = 0.9, size = 0.1, color = "black") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
        axis.ticks.x = element_blank(),
        #axis.title= element_text(size=8, face="bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        #legend.text = element_text(size = 8),
        #legend.key.size = unit(0.5, "cm"),
        legend.position = "top") +
  scale_x_discrete(expand = c(0, 0), limits = st_100_order_terrestrial$label) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_fill_manual(values = colors, labels = c("Unknown", "Agricultural soil", "Natural soil", "Forest soil", "Sediment", "Human feces", "Freshwater", "Coastal", "Opean ocean")) +
  ylab("Mixing proportions") +
  xlab("Samples") +
  labs(fill="")


ggpubr::ggarrange(emp_sources, st_100_barplot, ncol = 2, legend = "none", align = "h", widths = c(0.35, 0.65))

blue <- "#268bd2"
cyan <- "#2aa198"
green <- "#859900"

st_100_long_terrestrial_plot <- st_100_long %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  mutate(class = ifelse(variable == "saline_water_gt_25km", "marine", "terrestrial")) %>%
  dplyr::select(label, class, value) %>%
  group_by(label, class) %>%
  summarise(N=sum(value)) %>%
  filter(class == "marine") %>%
  ggplot(aes(label, 1, fill = (N))) +
  geom_tile(color = "white", alpha = 0.7) +
  scale_x_discrete(limits = st_100_order_terrestrial$label) +
  #scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_fill_gradient2(low = green, mid = cyan, high = blue, midpoint = 0.4) +
  ylab("") +
  xlab("") +
  labs(fill="") +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top")


osd95_eggnog_summary <- tbl(my_db, "osd95_eggnog_summary1") %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  group_by(label, class) %>%
  summarise(tot = sum(tot)) %>%
  ungroup() %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label), class = fct_relevel(class, c("SHAR", "EXCL")))


st_100_long_rgc_prop_plot_class <- ggplot(osd95_eggnog_summary, aes(label, tot, fill = class)) +
  geom_col(width = 1, alpha = 0.9, size = 0.1, color = "black") +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(expand = c(0, 0), limits = st_100_order_terrestrial$label) +
  scale_fill_manual(values = ggthemes::canva_pal(palette = "Simple and fresh")(4)[3:4]) +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
        axis.ticks.x = element_blank(),
        #axis.title= element_text(size=8, face="bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        #legend.text = element_text(size = 8),
        #legend.key.size = unit(0.5, "cm"),
        legend.position = "top") +
  ylab("Proportion") +
  xlab("Samples")


st_100_long_rgc_prop <- st_100_long %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  mutate(class = ifelse(variable == "saline_water_gt_25km", "marine", "terrestrial")) %>%
  dplyr::select(label, class, value) %>%
  group_by(label, class) %>%
  summarise(N=sum(value)) %>%
  filter(class == "marine") %>%
  arrange(label) %>%
  left_join(osd95_eggnog_summary %>% rename(cl_class = class) %>% filter(cl_class == "EXCL")) %>%
  ungroup() %>% left_join(osd2014_cdata)

fit1 <- lm(tot ~ N, data = st_100_long_rgc_prop)
summary(fit1)

colors <- c("#B99ED0", "#D2D0B7", "#C059D1", "#DB8763", "#ABE172", "#7AD4CB")

selected_mp <- c("Cold Temperate Northwest Atlantic", "Warm Temperate Northwest Atlantic", "Tropical Northwestern Atlantic",
                 "Northern European Seas" ,"Lusitanian", "Mediterranean Sea")

st_100_long_rgc_prop_plot <- ggplot(st_100_long_rgc_prop %>% group_by(meow_province) %>% mutate(tot_rank = rank(tot), N_rank = rank(N)) %>%
                                      filter(meow_province %in% osd2014_selected_meow_provinces$meow_province) %>% ungroup() %>%
                                      mutate(meow_province = fct_relevel(meow_province, selected_mp)), aes(N_rank, tot_rank, fill = meow_province), group = 1) +
  geom_smooth(method = "lm", color = "black", size = 0.5, fill = "grey20", alpha = 0.1) +
  geom_point(shape = 21, color = "black", alpha = 1, size = 1.8) +
  expand_limits(x = 0, y = 0) +
  #scale_x_continuous(labels = scales::percent) +#, limits = c(0, 0.8)) +
  #scale_y_continuous(labels = scales::percent) + #, limits = c(0, 0.6)) +
  theme_bw() +
  facet_wrap(~meow_province, scales = "free") +
  scale_fill_manual(values = colors) +
  xlab("Pelagic water influence") +
  ylab("OSD2014-specific ORFs proportion") + ggpubr::stat_cor(method = "spearman") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        axis.text = element_text(size = 6))



ggpubr::ggarrange(emp_sources, st_100_barplot, ncol = 2, nrow = 1, legend = "none", align = "hv", widths = c(0.35, 0.65))
ggsave(plot = last_plot(), filename = "osd2014_sourcetracker/figures/osd2014_st_barplot.pdf", width = 8.33, height = 3.5)

ggpubr::ggarrange(st_100_long_rgc_prop_plot_class, st_100_long_rgc_prop_plot, ncol = 2, nrow = 1, legend = "none", widths = c(0.35, 0.4))
ggsave(plot = last_plot(), filename = "osd2014_sourcetracker/figures/osd2014_st_orf_excl.pdf",  width = 6.33, height = 2.5)


st_100_long_latitude_plot <- st_100_long %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  mutate(class = ifelse(variable == "saline_water_gt_25km", "marine", "terrestrial")) %>%
  dplyr::select(label, class, value) %>%
  group_by(label, class) %>%
  summarise(N=sum(value)) %>%
  ungroup() %>%
  filter(class == "marine") %>%
  left_join(osd2014_cdata) %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label)) %>%
  ggplot(aes(label, 1, fill = (N))) +
  geom_tile(color = "white", alpha = 0.7) +
  geom_text(aes(label = meow_province), angle = 90, size = 2) +
  scale_x_discrete(limits = st_100_order_terrestrial$label) +
  #scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_fill_gradient2(low = green, mid = cyan, high = blue, midpoint = 0.4) +
  ylab("") +
  xlab("") +
  labs(fill="") +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top")


st_100_long_latitude_plot <- st_100_long %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  mutate(class = ifelse(variable == "saline_water_gt_25km", "marine", "terrestrial")) %>%
  dplyr::select(label, class, value) %>%
  group_by(label, class) %>%
  summarise(N=sum(value)) %>%
  ungroup() %>%
  filter(class == "marine") %>%
  left_join(osd2014_cdata) %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label)) %>%
  ggplot(aes(label, 1, fill = (N))) +
  geom_tile(color = "white", alpha = 0.7) +
  geom_text(aes(label = start_lat), angle = 90, size = 2) +
  scale_x_discrete(limits = st_100_order_terrestrial$label) +
  #scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_fill_gradient2(low = green, mid = cyan, high = blue, midpoint = 0.4) +
  ylab("") +
  xlab("") +
  labs(fill="") +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top")



colors <- c("#B99ED0", "#D2D0B7", "#C059D1", "#DB8763", "#ABE172", "#7AD4CB", "grey20")

selected_mp <- c("Cold Temperate Northwest Atlantic", "Warm Temperate Northwest Atlantic", "Tropical Northwestern Atlantic",
                 "Northern European Seas" ,"Lusitanian", "Mediterranean Sea", "Other")




st_100_long_mp_plot <- st_100_long %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  mutate(class = ifelse(variable == "saline_water_gt_25km", "marine", "terrestrial")) %>%
  dplyr::select(label, class, value) %>%
  group_by(label, class) %>%
  summarise(N=sum(value)) %>%
  ungroup() %>%
  filter(class == "marine") %>%
  left_join(osd2014_cdata) %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label),
         meow_province = ifelse(meow_province %in% selected_mp, meow_province, "Other"),
         meow_province = fct_relevel(meow_province, selected_mp)) %>%
  ggplot(aes(label, 1, fill = meow_province)) +
  geom_tile(color = "white", aes(alpha = log10(N))) +
  geom_text(aes(label = meow_province), angle = 90, size = 2) +
  scale_x_discrete(limits = st_100_order_terrestrial$label) +
  #scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_fill_manual(values = colors, labels = selected_mp) +
  ylab("") +
  xlab("") +
  labs(fill="") +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "none")



ggsave(st_100_barplot, filename = "osd2014_sourcetracker/figures/osd2014_ST_barplot.pdf", width = 11.69, height = 8.27)
ggsave(st_100_long_rgc_prop_plot, filename = "osd2014_sourcetracker/figures/osd2014_marine_excl_relation_barplot.pdf", width = 11.69, height = 8.27)
ggsave(st_100_long_terrestrial_plot, filename = "osd2014_sourcetracker/figures/osd2014_ST_tileplot.pdf", width = 11.69, height = 8.27)
ggsave(st_100_long_latitude_plot, filename = "osd2014_sourcetracker/figures/osd2014_ST_tileplot_lat.pdf", width = 11.69, height = 8.27)
ggsave(st_100_long_mp_plot, filename = "osd2014_sourcetracker/figures/osd2014_ST_tileplot_mp.pdf", width = 11.69, height = 8.27)

osd2014_dist_coast <- osd2014_cdata %>%
  #filter(label %in% osd2014_amp_mg_intersect$label) %>%
  dplyr::select(label, dist_coast_m) %>%
  mutate(dist_shore_km = dist_coast_m/1000) %>%
  arrange(desc(dist_shore_km))


osd2014_dist_coast$label <- factor(osd2014_dist_coast$label, levels = st_100_order_terrestrial$label)

dist_coast_plot <- ggplot(osd2014_dist_coast, aes(label, dist_shore_km, group = 1)) +
  geom_point(size = 0.25) +
  geom_line(size = 0.25) +
  theme_bw() +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top") +
  xlab("") +
  ylab("Distance from shore (km)")

ggsave(dist_coast_plot, filename = "osd2014_sourcetracker/figures/osd2014_distance_shore.pdf", width = 11.69, height = 8.27)


# Plot samples ------------------------------------------------------------
library(tidyverse)
library(ggrepel)


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

wmap <- map_data("world")
osd2014_map_all <- ggplot(wmap, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), fill = "#E5E6E6") +
  geom_path(aes(group = group), colour = "#E5E6E6") +
  geom_point(data = bind_rows( osd2014_cdata_coord_all, osd2014_cdata_coord_paper), aes(x = start_lon, y = start_lat, color = class),
             alpha = 0.8) +
  xlab("") +
  ylab("") +
  coord_equal(ratio = 1)  +
  scale_color_manual(values = c("#C02942", "#333333"), labels = c("NPL022 paper", "NPL022 all")) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")

ggsave(osd2014_map_all, filename = "osd2014_sourcetracker/figures/osd2014_map_all.pdf", width = 11.69, height = 8.27)



osd2014_map_paper <- ggplot(wmap, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), fill = "#E5E6E6") +
  geom_path(aes(group = group), colour = "#E5E6E6") +
  geom_point(data = bind_rows(osd2014_cdata_coord_paper), aes(x = start_lon, y = start_lat), alpha = 0.8, color = "#C02942") +
  xlab("") +
  ylab("") +
  coord_equal(ratio = 1)  +
  # scale_color_manual(values = c("#C02942", "#333333"), labels = c("NPL022 paper", "NPL022 all")) +
  theme(panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top")
ggsave(osd2014_map_paper, filename = "osd2014_sourcetracker/figures/osd2014_map_paper.pdf", width = 11.69, height = 8.27)

save.image(file = "osd2014_sourcetracker/data/osd2014_ST_analysis.Rdata", compress = TRUE)

