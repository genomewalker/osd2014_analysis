library(ALDEx2)
library(tidyverse)
library(phyloseq)
library(ggparl)


load("osd2014_shotgun/data/osd2014_emapper_rand_physeq_filt_objects.Rdata")
load("osd2014_shotgun/data/tara_emapper_rand_physeq_filt_objects.Rdata")


osd2014_emapper_rand_phyloseq_alpha_prop <- transform_sample_counts(osd2014_emapper_rand_phyloseq_alpha, function(x) x/sum(x))
tara_emapper_rand_phyloseq_alpha_prop <- transform_sample_counts(tara_emapper_rand_phyloseq_alpha, function(x) x/sum(x))

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label, water_depth <= 5, dist_coast_m <= 1000) %>%
  filter(grepl("Mediterranean Sea|Warm Temperate Northwest Atlantic|Tropical Northwestern Atlantic|Cold Temperate Northwest Atlantic|Lusitanian", meow_province))

tara_cdata <- tbl(my_db, "tara_metadata_prok") %>%
  collect(n = Inf) %>%
  filter(label %in% sample_names(tara_emapper_assm_phyloseq_alpha)) %>%
  filter(grepl("North Atlantic Ocean|Mediterranean", iho_label))

eggnog4_funcat <- tbl(my_db, "eggnog4_funcat") %>%
  collect(n = Inf)

# We only keep those results that fall in one category
osd2014_emapper_rand_phyloseq_long <- as(otu_table(osd2014_emapper_rand_phyloseq_alpha_prop), "matrix") %>%
  as_tibble(rownames = "group_nam") %>%
  gather(label, abundance, -group_nam) %>%
  group_by(label, group_nam) %>%
  dplyr::summarise(N = sum(abundance)) %>%
  filter(!is.na(group_nam)) %>%
  filter(label %in% osd2014_cdata$label) %>%
  inner_join(as(tax_table(osd2014_emapper_rand_phyloseq_alpha_prop), "matrix") %>% as_tibble())


osd2014_emapper_rand_phyloseq_long %>%
  group_by(label) %>%
  dplyr::summarise(p = sum(N)) %>%
  skimr::skim()


tara_emapper_rand_phyloseq_long <- as(otu_table(tara_emapper_rand_phyloseq_alpha_prop), "matrix") %>%
  as_tibble(rownames = "group_nam") %>%
  gather(label, abundance, -group_nam) %>%
  group_by(label, group_nam) %>%
  dplyr::summarise(N = sum(abundance)) %>%
  filter(!is.na(group_nam)) %>%
  filter(label %in% tara_cdata$label) %>%
  inner_join(as(tax_table(osd2014_emapper_rand_phyloseq_alpha_prop), "matrix") %>% as_tibble())

tara_emapper_rand_phyloseq_long %>%
  group_by(label) %>%
  dplyr::summarise(p = sum(N)) %>%
  skimr::skim()

osd2014_emapper_rand_phyloseq_beta_wide <- as(otu_table(osd2014_emapper_rand_phyloseq_alpha), "matrix") %>%
  as_tibble(rownames = "group_nam") %>%
  gather(label, abundance, -group_nam) %>%
  left_join( as(tax_table(tara_emapper_rand_phyloseq_alpha_prop), "matrix") %>%
               as_tibble()) %>%
  #filter(!(grepl('\\|',cog_description))) %>%
  filter(label %in% osd2014_cdata$label) %>%
  group_by(label, group_nam) %>%
  dplyr::summarise(N = sum(abundance)) %>%
  filter(!is.na(group_nam))

tara_emapper_rand_phyloseq_beta_wide <- as(otu_table(tara_emapper_rand_phyloseq_alpha), "matrix") %>%
  as_tibble(rownames = "group_nam") %>%
  gather(label, abundance, -group_nam) %>%
  left_join(as(tax_table(tara_emapper_rand_phyloseq_alpha_prop), "matrix") %>%
              as_tibble()) %>%
  #filter(!(grepl('\\|',cog_description)), group_nam %in% tara_osd_intersect_gnam) %>%
  filter(label %in% tara_cdata$label) %>%
  group_by(label, group_nam) %>%
  dplyr::summarise(N = sum(abundance)) %>%
  filter(!is.na(group_nam))


tara_osd_emapper_rand_phyloseq_beta_wide <- bind_rows(tara_emapper_rand_phyloseq_beta_wide,
                                                      osd2014_emapper_rand_phyloseq_beta_wide) %>%
  spread(label, N, fill = 0) %>%
  filter(!is.na(group_nam)) %>%
  as.data.frame() %>%
  column_to_rownames(var = "group_nam")



# DeSEQ2 analysis ---------------------------------------------------------
library(DESeq2)
physeq <- phyloseq(
  otu_table(tara_osd_emapper_rand_phyloseq_beta_wide, taxa_are_rows = TRUE),
  sample_data(bind_rows(osd2014_cdata %>% select(label) %>% mutate(study = "OSD"),
                        tara_cdata %>% select(label) %>% mutate(study = "TARA")) %>%
                mutate(study = as.factor(study)) %>%
                as.data.frame() %>% column_to_rownames("label"))
)
diagdds = phyloseq_to_deseq2(physeq, ~ study)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = as(sigtab, "data.frame")
head(sigtab)

sigtab <- sigtab %>%
  as.data.frame() %>%
  rownames_to_column("group_nam") %>%
  as_tibble()

x.sig_study_osd <- osd2014_emapper_rand_phyloseq_long %>% select(group_nam, label, N, cog_description) %>%
  filter(group_nam %in% sigtab$group_nam) %>%
  mutate(class = "OSD")

x.sig_study_tara <- tara_emapper_rand_phyloseq_long %>% select(group_nam, label, N, cog_description) %>%
  filter(group_nam %in% sigtab$group_nam) %>%
  mutate(class = "TARA")


kegg_sign <- bind_rows(x.sig_study_osd, x.sig_study_tara) %>%
  group_by(cog_description, class, label) %>%
  filter(!(is.na(cog_description))) %>%
  dplyr::summarise(N = sum(N)) %>%
  ungroup()

mu <- plyr::ddply(kegg_sign, c("class", "cog_description"), summarise, grp.mean=mean(N))

mu_tara <- mu %>% group_by(cog_description) %>% arrange(grp.mean) %>% top_n(1) %>% filter(class == "TARA")
mu_osd <- mu %>% group_by(cog_description) %>% arrange(grp.mean) %>% top_n(1) %>% filter(class == "OSD")


compare_means(data = kegg_sign %>% filter(cog_description %in% mu_tara$cog_description) %>% select(-label), formula = N ~ class, group.by = "cog_description")
compare_means(data = kegg_sign %>% filter(cog_description %in% mu_osd$cog_description) %>% select(-label), formula = N ~ class, group.by = "cog_description")


kegg_sign <- kegg_sign %>%
  tidyr::complete(cog_description, class) %>%
  ungroup() %>%
  #mutate(N = ifelse(N == 0, NA, N)) %>%
  mutate(functional_module = fct_relevel(cog_description,unique(kegg_sign$cog_description)))

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N1    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean(xx[[col]], na.rm=na.rm),
                     sd   = sd(xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N1)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N1-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}

tgc <- summarySE(kegg_sign, measurevar="N", groupvars=c("class","functional_module"))

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

p1 <- ggplot(kegg_sign %>% filter(N >0) %>% mutate(cog_description = fct_rev(cog_description)) ,aes(x = N, y = cog_description, fill = class)) +
  ggridges::geom_density_ridges(scale = 1, alpha = .5, rel_min_height = 0.0 ,panel_scaling = T, size = 0.2, quantile_lines = TRUE, quantiles = 2) +
  # geom_violin() +
  scale_x_log10(labels = scales::percent, breaks = base_breaks()) +
  scale_fill_manual(values = ggthemes::canva_pal(palette = "Stormy hues")(4)[3:4]) +
  scale_color_manual(values = ggthemes::canva_pal(palette = "Stormy hues")(4)[3:4]) +
  xlab("Proportion") +
  ylab("Functional processes") +
  #theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank())
p1

ggsave(plot = last_plot(), filename = "osd2014_shotgun/figures/osd2014_tara_kegg_comp_all.pdf", width = 11.69, height = 6)

# ggplot(kegg_sign, aes(N, fill = class,color = class)) +
#   geom_histogram(aes(y=..density..), alpha=0.5,
#                  position="identity", color = "black") +
#   geom_density(alpha=.2) +
#   facet_wrap(~cog_category, scale = "free", ncol = 6) +
#   scale_x_continuous(labels = scales::percent) +
#   theme_light() +
#   scale_fill_manual(values = ggthemes::canva_pal(palette = "Stormy hues")(4)[3:4]) +
#   scale_color_manual(values = ggthemes::canva_pal(palette = "Stormy hues")(4)[3:4])


# Draw map with samples ---------------------------------------------------

global <- map_data("world")

map_data <- bind_rows(osd2014_cdata %>% select(start_lat, start_lon) %>% mutate(class = "OSD"),
                      tara_cdata %>% select(latitude, longitude) %>% dplyr::rename(start_lat = latitude, start_lon = longitude) %>% mutate(class = "TARA"))

p2 <- ggplot() +
  geom_polygon(data = global, aes(x=long, y = lat, group = group)) +
  geom_point(data = map_data, aes(x = start_lon, y = start_lat, fill = class),
             size = 2.2, shape = 21, color = "black", alpha = 1) +
  xlab("") +
  ylab("") +
  coord_equal(xlim = c(-100, 50), ylim = c(0, 70)) +
  scale_fill_manual(values = ggthemes::canva_pal(palette = "Stormy hues")(4)[3:4]) +
  #scale_fill_manual(values = c("#C02942", "#2A67A0"), labels = c("NPL022 paper", "NPL022 all")) +
  theme_map(legend.position = "bottom")
#scale_fill_manual(values = colors) +
#scale_shape_manual(values = c(21,22),  labels = c("NPL022 paper", "NPL022 all"))


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
ggarrange(p1, p2, widths = c(0.6, 0.4))

ggsave(plot = last_plot(), filename = "osd2014_shotgun/figures/osd2014_tara_kegg_comp_all.pdf", width = 11, height = 8)

















# ALDeX2 analysis ---------------------------------------------------------



names_order <- c(grep("TARA", colnames(tara_osd_emapper_rand_phyloseq_beta_wide),value = TRUE), grep("OSD", colnames(tara_osd_emapper_rand_phyloseq_beta_wide),value = TRUE))

tara_osd_emapper_rand_phyloseq_beta_wide <- tara_osd_emapper_rand_phyloseq_beta_wide[,names_order]


conds <- c(rep("1_TARA", length(grep("TARA", colnames(tara_osd_emapper_rand_phyloseq_beta_wide)))), rep("2_OSD", length(grep("OSD", colnames(tara_osd_emapper_rand_phyloseq_beta_wide)))))

#conds <- rev(conds)

#y <- aldex.clr((tara_osd_eggnog4_prop), conds, mc.samples = 128, denom = "all", verbose = TRUE, useMC = TRUE)
x.all<-aldex(reads = tara_osd_emapper_rand_phyloseq_beta_wide, conditions = conds, mc.samples = 128, test = "t", effect = TRUE, include.sample.summary = TRUE,verbose = TRUE, denom = "iqlr")

# positive numbers mean more abundant in OSD, negative numbers more abundant in TARA
aldex.plot(x.all, test="wilcox")
x.sig <- x.all[x.all$wi.eBH < 0.01,] %>% dplyr::select((contains("sample")))

x.all$group_nam <- row.names(x.all)
x.sig <- x.all[x.all$wi.eBH < 0.001,] %>% dplyr::select(-(contains("sample")))
x.sig$group_nam <- row.names(x.sig)

x.sig <- x.all %>%
  filter(wi.eBH < 0.001) %>%
  dplyr::select(-(contains("sample")))

sig.tara <- rownames(x.all)[which((x.all$wi.eBH < 0.001 & x.all$effect < 2))]
sig.osd <- rownames(x.all)[which((x.all$wi.eBH < 0.001 & x.all$effect > 2))]


x.sig$group_nam <- factor(x.sig$group_nam, levels = rev(sort(x.sig$group_nam %>% unique())))

x.sig_study <- x.sig %>%
  as_tibble() %>%
  mutate(class = ifelse(diff.btw > 0, "OSD", "TARA"))


x.sig_study_osd <- osd2014_emapper_rand_phyloseq_long %>% select(group_nam, label, N, cog_description) %>%
  filter(group_nam %in% x.sig_study$group_nam) %>%
  mutate(class = "OSD")

x.sig_study_tara <- tara_emapper_rand_phyloseq_long %>% select(group_nam, label, N, cog_description) %>%
  filter(group_nam %in% x.sig_study$group_nam) %>%
  mutate(class = "TARA")


kegg_sign <- bind_rows(x.sig_study_osd, x.sig_study_tara) %>%
  group_by(cog_description, class, label) %>%
  filter(!(is.na(cog_description))) %>%
  dplyr::summarise(N = sum(N)) %>%
  ungroup()

mu <- plyr::ddply(kegg_sign, c("class", "cog_description"), summarise, grp.mean=mean(N))

mu_tara <- mu %>% group_by(cog_description) %>% arrange(grp.mean) %>% top_n(1) %>% filter(class == "TARA")
mu_osd <- mu %>% group_by(cog_description) %>% arrange(grp.mean) %>% top_n(1) %>% filter(class == "OSD")


compare_means(data = kegg_sign %>% filter(cog_description %in% mu_tara$cog_description) %>% select(-label), formula = N ~ class, group.by = "cog_description")
compare_means(data = kegg_sign %>% filter(cog_description %in% mu_osd$cog_description) %>% select(-label), formula = N ~ class, group.by = "cog_description")


kegg_sign <- kegg_sign %>%
  tidyr::complete(cog_description, class) %>%
  ungroup() %>%
  #mutate(N = ifelse(N == 0, NA, N)) %>%
  mutate(cog_description = fct_relevel(cog_description,unique(kegg_sign$cog_description)))

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N1    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N1)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N1-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}

tgc <- summarySE(kegg_sign, measurevar="N", groupvars=c("class","cog_description"))

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

ggplot(kegg_sign %>% filter(N >0) %>% mutate(cog_description = fct_rev(cog_description)) ,aes(x = N, y = cog_description, fill = class)) +
  ggridges::geom_density_ridges(scale = 1, alpha = .5, rel_min_height = 0.0 ,panel_scaling = T, size = 0.2, quantile_lines = TRUE, quantiles = 2) +
  # geom_violin() +
  scale_x_log10(labels = scales::percent, breaks = base_breaks()) +
  scale_fill_manual(values = ggthemes::canva_pal(palette = "Stormy hues")(4)[3:4]) +
  scale_color_manual(values = ggthemes::canva_pal(palette = "Stormy hues")(4)[3:4]) +
  xlab("Proportion") +
  ylab("Functional processes") +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
  theme(legend.position = "top",
        legend.title = element_blank())


