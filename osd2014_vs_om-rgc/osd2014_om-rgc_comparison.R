library(tidyverse)
library(RPostgreSQL)


# OSD2014 and OMR-RGC incremental clustering with CD-HIT ------------------------------


# Draw Venn diagram with the interection between OSD and OM-RGC

#6647689 -- 2173514
#        |_ 4474175 -- 3542383 (-- 3427068)

OMRGC <- 40154822
OSD <- 3403948 + 865298
SHARED <- 865298

library(Vennerable)
VenTest <- Venn(SetNames= c("OM-RGC","OSD"), Weight=c(OMRGC+OSD,OMRGC-SHARED, OSD-SHARED, SHARED ))
plot(VenTest)
pdf(file="osd2014_vs_om-rgc/figures/osd2014_incremental.pdf")
plot(VenTest)
dev.off()

# OSD2014 and OM-RGC de-novo clustering with MMSEQS2 ----------------------

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

osd95 <- tbl(my_db, "osd2014_assm_cdhit_95")%>%
  collect(n=Inf)

osd2014_eggnog <- tbl(my_db, "osd2014_assm_emapper_results")%>%
  collect(n=Inf)

eggnog_funcat <- tbl(my_db, "eggnog4_funcat") %>%
  collect(n=Inf)

eggnog_annotation <- tbl(my_db, "eggnog4_annotation") %>%
  collect(n=Inf)

eggnog_funcat <- eggnog_annotation %>%
  left_join(eggnog_funcat)

osd2014_eggnog <- osd2014_eggnog %>%
  select(-cog_functional_category) %>%
  left_join(eggnog_funcat)

osd2014_abun <- tbl(my_db, "osd2014_orf_abundance") %>%
  collect(n=Inf)

osd95_abun <- osd95 %>%
  left_join(osd2014_abun %>% rename(id = gene_id)) %>%
  tbl_df()

rm(osd95)
rm(osd2014_abun)
gc()

osd2014_omrgc_mmseqs <- tbl(my_db, "osd2014_omrgc_mmseqs") %>%
  collect(n=Inf)

osd2014_omrgc_mmseqs_prop <- osd2014_omrgc_mmseqs %>%
  mutate(V1 = ifelse(rep == members, TRUE, FALSE),
         class = ifelse(grepl("OSD", members), "OSD", "OM")) %>%
  group_by(cluster, class) %>%
  count() %>%
  ungroup() %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup()


# 10280221 mmseqs clusters
# 1349249 OSD mmseqs clusters
# 8143533 OM mmseqs clusters
# 787439 SHAR mmseqs clusters


OMRGC <- 8143533
OSD <- 1349249 + 787439
SHARED <- 787439

library(Vennerable)
VenTest <- Venn(SetNames= c("OM-RGC","OSD"), Weight=c(OMRGC+OSD,OMRGC-SHARED, OSD-SHARED, SHARED ))
plot(VenTest)
pdf(file="osd2014_vs_om-rgc/figures/osd2014_mmseqs_clusters.pdf")
plot(VenTest)
dev.off()


mmseqs_prop_osd <- osd2014_omrgc_mmseqs_prop %>%
  filter(class == "OSD", prop == 1.0)

mmseqs_prop_om <- osd2014_omrgc_mmseqs_prop %>%
  filter(class == "OM", prop == 1.0)

mmseqs_prop_shared <- osd2014_omrgc_mmseqs_prop %>%
  filter(prop != 1.0)

# All ORFs = 6,647,689
# All ORFs at 95% = 4,413,151
# RGC = 4,269,246
# Removed from being to short = 144,984
# 2,833,809 + 1,435,437 = 4,269,246
# 1,435,437
mmseqs_prop_osd_orfs <- osd2014_omrgc_mmseqs %>%
  filter(cluster %in% mmseqs_prop_osd$cluster)
# 2,833,809
mmseqs_prop_shared_orfs <- osd2014_omrgc_mmseqs %>%
  filter(cluster %in% mmseqs_prop_shared$cluster) %>%
  filter(grepl("OSD", members))


rm(osd2014_omrgc_mmseqs)
gc()

# Now we need to get the ORFs that were clustered before at 95%
# The ones in the RGC are representatives

# 1,582,385
osd95_osd_clstr <- osd95_abun %>%
  filter(id %in% mmseqs_prop_osd_orfs$members) %>%
  left_join(mmseqs_prop_osd_orfs %>% rename(id = members))
osd95_osd <- osd95_abun %>%
  filter(clstr %in% osd95_osd_clstr$clstr) %>%
  left_join(osd95_osd_clstr %>% select(clstr, cluster))

# 4,920,320
osd95_shared_clstr <- osd95_abun %>%
  filter(id %in% mmseqs_prop_shared_orfs$members) %>%
  left_join(mmseqs_prop_shared_orfs %>% rename(id = members))
osd95_shared <- osd95_abun %>%
  filter(clstr %in% osd95_shared_clstr$clstr) %>%
  left_join(osd95_shared_clstr %>% select(clstr, cluster))


osd_orfs_cl <- bind_rows(osd95_osd %>%
                           select(cluster, label, abun) %>%
                           group_by(cluster, label) %>%
                           summarise(N = sum(abun)) %>%
                           ungroup %>%
                           mutate(class = "EXCL"),
                         osd95_shared %>%
                           select(cluster, label, abun) %>%
                           group_by(cluster, label) %>%
                           summarise(N = sum(abun)) %>%
                           ungroup %>%
                           mutate(class = "SHAR"))

osd_orfs_cl$cluster %>% unique %>% length()

osd_orfs_cl_prop <- osd_orfs_cl %>%
  group_by(label) %>%
  mutate(T = sum(N), prop = N/T) %>%
  ungroup() %>%
  group_by(cluster) %>%
  mutate(occ = n())

all_cl_stats <- osd_orfs_cl_prop %>%
  select(cluster, class, prop, occ) %>%
  group_by(cluster, class, occ) %>%
  summarise(N = mean(prop)) %>%
  ungroup()

all_cl_stats_excl <- all_cl_stats %>%
  ungroup() %>%
  filter(class == "EXCL") %>%
  arrange(desc(N)) %>%
  mutate(rank = row_number())

all_cl_stats_shar <- all_cl_stats %>%
  ungroup() %>%
  filter(class == "SHAR") %>%
  arrange(desc(N)) %>%
  mutate(rank = row_number())

p1 <- ggplot(all_cl_stats %>% filter(N > 0) %>% unique %>% mutate(class = fct_relevel(class, c("EXCL", "SHAR")), aes(occ, N, color = class, fill = class)) +
  geom_hex(aes(alpha=log(..count..)), bins = 60) +
  scale_y_log10(breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1), labels = function(x) paste0(x*100, "%")) +
  scale_x_log10() +
  theme_bw() +
  xlab("Occurrence") +
  ylab("Proportion") +
  guides(alpha=guide_legend(title="# clusters", label=T)) +
  scale_fill_manual(values = (c("#29202A", "#C44D33"))) +
  scale_color_manual(values = (c("#29202A", "#C44D33")))


p2 <- ggplot(bind_rows(all_cl_stats_excl, all_cl_stats_shar) %>% filter(N > 0) %>% mutate(class = fct_relevel(class, c("EXCL", "SHAR"))),aes(x=rank,y=N)) +
  geom_line(aes(color=class),) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1), labels = function(x) paste0(x*100, "%")) +
  ylab("Proportion") +
  xlab("Rank") +
  scale_fill_manual(values = (c("#29202A", "#C44D33"))) +
  scale_color_manual(values = (c("#29202A", "#C44D33"))) +
  theme_bw() +
  theme(legend.position = "none")

p_rabun <- ggpubr::ggarrange(p2, p1, labels = c("A", "B"), legend = "none")

ggsave("osd2014_vs_om-rgc/figures/osd2014_assm_orf_rabun_occ.pdf", width = 8.33, height = 3.5)
leg <- ggpubr::get_legend(p1)
ggsave("osd2014_vs_om-rgc/figures/osd2014_assm_orf_rabun_occ_legend.pdf", plot = ggpubr::as_ggplot(leg), width = 8.33, height = 3.5)



# Explore the eggNOG categories from exclusive and shared -----------------

osd95_eggnog <- bind_rows(osd95_osd  %>%
                            left_join(osd2014_eggnog %>% rename(id = query_name)) %>%
                            select(id, cluster, label, abun, cog_functional_category, cog_description, cog_category) %>%
                            unique() %>%
                            mutate(class = "EXCL"),
                          osd95_shared %>%
                            left_join(osd2014_eggnog %>% rename(id = query_name)) %>%
                            select(id, cluster, label, abun, cog_functional_category, cog_description, cog_category) %>%
                            unique() %>%
                            mutate(class = "SHAR"))


counts_sample <- osd95_eggnog %>%
  select(label, abun) %>%
  group_by(label) %>%
  summarise(Total = sum(abun))

orfs_sample <- osd95_eggnog %>%
  select(label) %>%
  group_by(label) %>%
  summarise(n_orfs = n())

st_100_order_terrestrial <- tbl(my_db, "osd2014_st_order_terrestrial") %>%
  collect(n = Inf)

osd95_eggnog_summary <- osd95_eggnog %>%
  select(label, class, cog_description, abun) %>%
  mutate(cog_description = ifelse(is.na(cog_description), "Not annotated", cog_description)) %>%
  group_by(label, class, cog_description) %>%
  summarise(N = sum(abun), n_orfs = n()) %>%
  arrange(label) %>%
  ungroup() %>%
  group_by(label) %>%
  mutate(Total = sum(N), prop = N/Total) %>%
  ungroup() %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label),
         cog_description = fct_relevel(cog_description, c("Translation, ribosomal structure and biogenesis",
                                                          "RNA processing and modification",
                                                          "Transcription",
                                                          "Replication, recombination and repair",
                                                          "Chromatin structure and dynamics",
                                                          "Cell cycle control, cell division, chromosome partitioning",
                                                          "Nuclear structure",
                                                          "Defense mechanisms",
                                                          "Signal transduction mechanisms",
                                                          "Cell wall/membrane/envelope biogenesis",
                                                          "Cell motility",
                                                          "Cytoskeleton",
                                                          "Extracellular structures",
                                                          "Intracellular trafficking, secretion, and vesicular transport",
                                                          "Posttranslational modification, protein turnover, chaperones",
                                                          "Energy production and conversion",
                                                          "Carbohydrate transport and metabolism",
                                                          "Amino acid transport and metabolism",
                                                          "Nucleotide transport and metabolism",
                                                          "Coenzyme transport and metabolism",
                                                          "Lipid transport and metabolism",
                                                          "Inorganic ion transport and metabolism",
                                                          "Secondary metabolites biosynthesis, transport and catabolism",
                                                          #"General function prediction only",
                                                          "Function unknown",
                                                          "Not annotated")))




osd95_eggnog_summary$cog_description <- plyr::mapvalues(osd95_eggnog_summary$cog_description,
                                              c("Translation, ribosomal structure and biogenesis",
                                                "RNA processing and modification",
                                                "Transcription",
                                                "Replication, recombination and repair",
                                                "Chromatin structure and dynamics",
                                                "Cell cycle control, cell division, chromosome partitioning",
                                                "Nuclear structure",
                                                "Defense mechanisms",
                                                "Signal transduction mechanisms",
                                                "Cell wall/membrane/envelope biogenesis",
                                                "Cell motility",
                                                "Cytoskeleton",
                                                "Extracellular structures",
                                                "Intracellular trafficking, secretion, and vesicular transport",
                                                "Posttranslational modification, protein turnover, chaperones",
                                                "Energy production and conversion",
                                                "Carbohydrate transport and metabolism",
                                                "Amino acid transport and metabolism",
                                                "Nucleotide transport and metabolism",
                                                "Coenzyme transport and metabolism",
                                                "Lipid transport and metabolism",
                                                "Inorganic ion transport and metabolism",
                                                "Secondary metabolites biosynthesis, transport and catabolism",
                                                #"General function prediction only",
                                                "Function unknown",
                                                "Not annotated"),
                                              c("Translation, ribosomal struct. and biogenesis",
                                                "RNA processing and modification",
                                                "Transcription",
                                                "Replication, recombination and repair",
                                                "Chromatin structure and dynamics",
                                                "Cell cycle control, cell division, chrom. part.",
                                                "Nuclear structure",
                                                "Defense mechanisms",
                                                "Signal transduction mechanisms",
                                                "Cell wall/membrane/envelope biogenesis",
                                                "Cell motility",
                                                "Cytoskeleton",
                                                "Extracellular structures",
                                                "Intracellular trafficking, secr., and vesic. transp.",
                                                "Posttrans. mod., protein turnover, chaperones",
                                                "Energy production and conversion",
                                                "Carbohydrate transport and metabolism",
                                                "Amino acid transport and metabolism",
                                                "Nucleotide transport and metabolism",
                                                "Coenzyme transport and metabolism",
                                                "Lipid transport and metabolism",
                                                "Inorganic ion transport and metabolism",
                                                "Secondary metab. biosyn., transp. and catab.",
                                                #"General function prediction only",
                                                "Function unknown",
                                                "Not annotated")
)


osd95_eggnog_summary1 <- osd95_eggnog_summary %>%
  group_by(cog_description) %>%
  mutate(all = sum(prop)) %>%
  ungroup() %>%
  mutate(cog_description1 = ifelse(cog_description == "Not annotated" | cog_description == "Function unknown", "Unknown function", "Known Function"),
         class = fct_relevel(class, c("SHAR", "EXCL"))) %>%
  group_by(label, class, cog_description1) %>%
  summarise(tot = sum(prop)) %>%
  ungroup() %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label))


ggplot(osd95_eggnog_summary1, aes(label, tot, fill = cog_description1)) +
  geom_col(width = 1, color = "white", size = 0.1) +
  theme_light() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  facet_wrap(~class, nrow = 1, ncol = 2) +
  scale_fill_manual(values = ggthemes::canva_pal(palette = "Light and natural")(2)) +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) +
  ylab("Proportion") +
  xlab("Samples")

ggsave(plot = last_plot(), filename = "osd2014_vs_om-rgc/figures/osd2014_mmseqs_clusters_eggnog.pdf", width = 8.33, height = 3.5)

# Let’s see which percentage of the ORFs might be spurious ----------------
osd2014_shadow_orfs <- tbl(my_db, "osd2014_assm_shadow_orfs") %>%
  collect(n=Inf)

# ORFs with same length: 386
shadow_orfs_same <- osd2014_shadow_orfs %>%
  filter(shadow == "SAME")
shadow_orfs_same <- c(shadow_orfs_same$orf.x, shadow_orfs_same$orf.y) %>% unique()
shadow_orfs_same %>% length()

c(shadow_orfs_same,orfs_in_shadow) %>% unique() %>% length()

# Shadow ORFs no NA: 67336
shadow_orfs_no_same <- osd2014_shadow_orfs %>%
  filter(shadow != "SAME")
shadow_orfs_no_same <- c(shadow_orfs_no_same$orf.x, shadow_orfs_no_same$orf.y) %>% unique()
shadow_orfs_no_same %>% length()

# 386
shadow_in_clstr_same <- osd95_eggnog %>%
  filter(id %in% shadow_orfs_same) %>%
  left_join(bind_rows(osd95_osd, osd95_shared)) %>%
  select(id, cluster, label, cog_description, clstr_size, class)

# 67336
shadow_in_clstr <- osd95_eggnog %>%
  filter(id %in% shadow_orfs_no_same) %>%
  left_join(bind_rows(osd95_osd, osd95_shared)) %>%
  select(id, cluster, label, cog_description, clstr_size, class)

knitr::kable(shadow_in_clstr_same %>% group_by(class) %>% count())

knitr::kable(shadow_in_clstr %>% group_by(class) %>% count())


save.image(file = "osd2014_vs_om-rgc/data/osd2014_om-rgc_comparison.Rdata", compress = TRUE)
save(all_cl_stats, all_cl_stats_excl, all_cl_stats_shar, file = "osd2014_vs_om-rgc/data/osd2014_om-rgc_comparison_2plot.Rdata", compress = TRUE)

