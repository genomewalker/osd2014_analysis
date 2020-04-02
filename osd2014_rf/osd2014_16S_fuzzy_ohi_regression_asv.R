library(DESeq2)
library(tidyverse)
library(phyloseq)
library(fuzzyforest)
library(igraph)
library(tidygraph)
library(caret)
library(ggpubr)
source("osd2014_16S_asv/lib/fuzzyforest_lib.R")

load("osd2014_16S_asv/data/osd2014_16S_asv_networks_results.Rdata", verbose = TRUE)
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata")
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

# Get OHI data and define impacted/non-impacted ---------------------------
osd2014_order_terrestrial <- tbl(my_db, "osd2014_st_order_terrestrial") %>%
  collect(n = Inf)

osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)

osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)

osd2014_silva_dada2_names <- tbl(my_db, "osd2014_silva_dada2") %>%
  collect(n = Inf) %>%
  select(asv, asv_name)

# Get Halpern data
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


# Select those samples that are impacted and at least in 2km from the shore
halpern_impact <- out_rescaled_2013_median_long %>%
  filter(buffer == "5km", ohi_variable == "global_cumul_impact") %>%
  mutate(class = ifelse(median <= 0, "non_impacted",
                        ifelse((median > 0.01 & median <= 4.8), "low_impacted",
                               "impacted"))) %>%
  filter(median != 0, class != "non_impacted")  %>%
  filter(class == "low_impacted" | (class == "impacted")) %>% # & median > 6)) %>%
  left_join(osd2014_cdata) %>%
  #filter(dist_coast_m <= 5000) %>%
  filter(label %in% osd2014_amp_mg_intersect$label)


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

halpern_impact_class_min <- halpern_impact %>%
  group_by(class) %>%
  dplyr::count() %>%
  ungroup() %>%
  dplyr::slice(which.min(n)) %>%
  .$n

# sites_rand_non <- halpern_impact %>%
#   filter(class == "non_impacted") %>%
#   .$label %>%
#   sample(.,halpern_impact_class_min)

sites_rand_low <- halpern_impact %>%
  filter(class == "low_impacted") %>%
  sample_n(halpern_impact_class_min) %>%
  .$label

sites_rand_imp <- halpern_impact %>%
  filter(class == "impacted") %>%
  arrange(desc(median)) %>%
  sample_n(halpern_impact_class_min) %>%
  .$label


sites <- unique(c(sites_rand_low, sites_rand_imp))


# # check stressors
# sites_rand_imp <- out_rescaled_2013_median_long %>%
#   filter(buffer == "1km", label %in% sites_rand_imp, ohi_variable != "sst", ohi_variable != "uv") %>%
#   left_join(osd2014_metadata)
#
# sites_rand_imp$label <- factor(sites_rand_imp$label, levels = sites_rand_imp %>% filter(ohi_variable == "global_cumul_impact") %>% arrange(median) %>% .$label %>% unique)
#
# ggplot(sites_rand_imp, aes(label, median)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~ohi_variable, scales = "free") +
#   theme_light() +
#   theme()
#
# sites_rand_low<- out_rescaled_2013_median_long %>%
#   filter(buffer == "1km", label %in% sites_rand_low, ohi_variable != "sst", ohi_variable != "uv") %>%
#   left_join(osd2014_metadata)
#
# sites_rand_low$label <- factor(sites_rand_low$label, levels = sites_rand_low %>% filter(ohi_variable == "global_cumul_impact") %>% arrange(median) %>% .$label %>% unique)
#
# ggplot(sites_rand_low, aes(label, median)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~ohi_variable, scales = "free")



osd2014_dada2_phyloseq_alpha_filt <- subset_samples(osd2014_dada2_phyloseq_alpha, label %in% halpern_impact$label)
osd2014_dada2_phyloseq_alpha_filt <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_alpha_filt) > 0, osd2014_dada2_phyloseq_alpha_filt)
osd2014_dada2_phyloseq_alpha_filt_h <- microbiome::transform(osd2014_dada2_phyloseq_alpha_filt, transform = "hellinger")

osd2014_dada2_phyloseq_alpha_prop <- transform_sample_counts(osd2014_dada2_phyloseq_alpha_filt, function(x) x/sum(x))
tax.mean <- taxa_sums(osd2014_dada2_phyloseq_alpha_prop)/nsamples(osd2014_dada2_phyloseq_alpha_prop)
osd2014_dada2_phyloseq_beta_filt <- prune_taxa(tax.mean > 2e-5, osd2014_dada2_phyloseq_alpha_filt_h)

osd2014_dada2_phyloseq_beta_filt <- prune_taxa(taxa_names(osd2014_dada2_phyloseq_beta), osd2014_dada2_phyloseq_beta_filt)
osd2014_dada2_phyloseq_beta_filt <- prune_taxa(df_nodes$asv, osd2014_dada2_phyloseq_beta_filt)
#osd2014_dada2_phyloseq_beta_filt <- microbiome::transform(osd2014_dada2_phyloseq_beta_filt, transform = "hellinger")






osd2014_dada2_phyloseq_alpha_filt <- subset_samples(osd2014_dada2_phyloseq_alpha, label %in% halpern_impact$label)
osd2014_dada2_phyloseq_alpha_filt <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_alpha_filt) > 0, osd2014_dada2_phyloseq_alpha_filt)
rare_taxa <- microbiome::rare_members(osd2014_dada2_phyloseq_alpha_filt, prevalence = 20/100)
osd2014_dada2_phyloseq_alpha_norare <- microbiome::remove_taxa(rare_taxa, osd2014_dada2_phyloseq_alpha_filt)
#osd2014_dada2_phyloseq_alpha_norare <- proportion(osd2014_dada2_phyloseq_alpha_norare)

osd2014_dada2_phyloseq_beta_filt <- prune_taxa(df_nodes$asv, osd2014_dada2_phyloseq_alpha_norare)
osd2014_dada2_phyloseq_beta_filt <- microbiome::transform(osd2014_dada2_phyloseq_beta_filt, transform = "compositional")
osd2014_dada2_phyloseq_beta_filt <- microbiome::transform(osd2014_dada2_phyloseq_beta_filt, transform = "log10p")




# START -------------------------------------------------------------------
#
# asv_otutable <- phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_filt)
# asv_otutable_pa <- vegan::decostand(data.frame(impacted = colSums(asv_otutable[sites_rand_imp,]), low = colSums(asv_otutable[sites_rand_low,])), method = "pa")
# dim(asv_otutable_pa)
#
# asv_otutable <- phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_filt)
#
# UpSetR::upset(asv_otutable_pa, order.by = c("degree"), decreasing = c(TRUE))
#
# asv_shared_impact_low <- dplyr::as_data_frame(asv_otutable_pa) %>%
#   rownames_to_column(var = "asv") %>%
#   filter(impacted == 1, low == 1) %>%
#   tbl_df()
#
# asv_impact_low <- dplyr::as_data_frame(asv_otutable_pa) %>%
#   rownames_to_column(var = "asv") %>%
#   filter(impacted == 0, low == 1) %>%
#   tbl_df()
#
# asv_impact_imp <- dplyr::as_data_frame(asv_otutable_pa) %>%
#   rownames_to_column(var = "asv") %>%
#   filter(impacted == 1, low == 0) %>%
#   tbl_df()
#
# asv_shared_impact_low_occ <- dplyr::as_data_frame(asv_otutable) %>%
#   rownames_to_column(var = "label") %>%
#   reshape2::melt() %>%
#   tbl_df %>%
#   rename(asv = variable) %>%
#   filter(asv %in% asv_shared_impact_low$asv) %>%
#   filter(value > 0) %>%
#   group_by(asv) %>%
#   count() %>%
#   mutate(occ = n/nsamples(osd2014_dada2_phyloseq_beta_filt)) %>%
#   filter(occ >= 0.25)
#
#
# #otutable_filt <- asv_otutable[,c(as.character(asv_shared_impact_low_occ$asv), asv_impact_low$asv, asv_impact_imp$asv)]
# #asv_otutable <- decostand(asv_otutable, "hellinger")
# asv_otutable_pa <- as.data.frame(asv_otutable)
#
#
# #names(asv_otutable_pa) <- asv_name_conv$new
# asv_otutable_pa$label <- rownames(asv_otutable_pa)
# asv_otutable_pa <- asv_otutable_pa  %>% left_join(halpern_impact %>% dplyr::select(label, class))
# class_impact <- asv_otutable_pa %>% dplyr::select(label, class)
# base::row.names(asv_otutable_pa) <- asv_otutable_pa$label
# asv_otutable_pa$label <- NULL
# asv_otutable_filt_pa <- asv_otutable_pa[sites,c(as.character(asv_shared_impact_low_occ$asv), asv_impact_low$asv, asv_impact_imp$asv,"class")]
#
# #samp <- sample(nrow(asv_otutable), 2/3 * nrow(asv_otutable))
# asv_otutable_filt_pa$class <- as.factor(asv_otutable_filt_pa$class)
# dim(asv_otutable_filt_pa)
# END ---------------------------------------------------------------------




#osd2014_dada2_phyloseq_beta_filt <- subset_samples(osd2014_dada2_phyloseq_beta_vst, label %in% osd2014_cdata$label)
#osd2014_dada2_phyloseq_beta_filt <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_beta_filt) > 0, osd2014_dada2_phyloseq_beta_filt)
osd2014_dada2_phyloseq_beta_df <- (as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix")) %>% as_tibble(rownames = "label") %>%
  inner_join(halpern_impact %>%
               select(label, median)) %>%
  as.data.frame() %>%
  column_to_rownames("label")
dim(osd2014_dada2_phyloseq_beta_df)


osd2014_dada2_phyloseq_beta_df$median <- as.numeric(osd2014_dada2_phyloseq_beta_df$median)
dim(osd2014_dada2_phyloseq_beta_df)
# osd2014_dada2_phyloseq_beta_df<-asv_otutable_filt_pa

# Run fuzzforest with 5 repetitions and 10 k-folds ------------------------

ff_training <- vector(mode = "list")
ff_fit_r <- vector(mode = "list")
ff_predict <- vector(mode = "list")
ff_accuracy <- vector(mode = "list")
set.seed(11)
for (i in 1:5){
  training.samples <- osd2014_dada2_phyloseq_beta_df$median %>%
    createFolds(returnTrain = TRUE)
  #createDataPartition(p = 0.8, list = FALSE)
  ff_training[[i]] <- training.samples
  for (o in 1:10){
    cat(paste0("Rep: ", i, " Fold: ", o, " # Start\n"))
    train.data <- osd2014_dada2_phyloseq_beta_df[training.samples[[o]], ]
    test.data <- osd2014_dada2_phyloseq_beta_df[-training.samples[[o]], ]

    module_membership <- df_nodes %>%
      filter(asv %in% colnames(osd2014_dada2_phyloseq_beta_df)) %>%
      dplyr::slice(match(colnames(train.data), asv)) %>% .$com


    mtry_factor <- 1; min_ntree <- 500;  drop_fraction <- .5; ntree_factor <- 1
    nodesize <- 1; final_ntree <- 500

    screen_params <- screen_control(drop_fraction = drop_fraction,
                                    keep_fraction = .25, min_ntree = min_ntree,
                                    ntree_factor = ntree_factor,
                                    mtry_factor = mtry_factor)
    select_params <- select_control(drop_fraction = drop_fraction,
                                    number_selected = 500,
                                    min_ntree = min_ntree,
                                    ntree_factor = ntree_factor,
                                    mtry_factor = mtry_factor)
    train.data$median <- as.numeric(train.data$median)
    cat(paste0("Rep: ", i, " Fold: ", o, " # FuzzyForest\n"))
    ff_fit <- ff(train.data[,1:ncol(train.data) - 1], train.data[,ncol(train.data)], module_membership = module_membership,
                 screen_params = screen_params, select_params=select_params,
                 final_ntree = 500, num_processors = 1)

    ff_fit_r[[paste0("iter_",i,"_fold_",o)]] <- ff_fit
    cat(paste0("Rep: ", i, " Fold: ", o, " # Predict\n"))
    pred_asv <- predict(ff_fit$final_rf, newdata = test.data[,1:ncol(test.data) - 1])
    ff_predict[[paste0("iter_",i,"_fold_",o)]] <- pred_asv
    test.data$rightPred <- pred_asv == test.data$median
    cat(paste0("Rep: ", i, " Fold: ", o, " # RMSE\n"))
    rmse <- RMSE(pred_asv, test.data$median)
    ff_rmse[[paste0("iter_",i,"_fold_",o)]] <- rmse
    cat(paste0("Rep: ", i, " Fold: ", o, " # Done\n"))
  }
}


get_rmse <- function(X){
  ff_rmse[[X]] %>% as_tibble() %>% mutate(run = names(ff_predict[X]))
}
bind_rows(map(1:50, get_rmse)) %>% ggplot(aes(y = value, x = "accuracy")) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 1,
                        color = "black", alpha = 1, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4, errorbar.length = 0.2) +
  theme_bw() +
  ylab("RMSE") +
  xlab("")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Get top features --------------------------------------------------------
# we use the broken stick model to select the most important features, then for each repetition and fold we select
# those that occur more than 40%
top_features <- map_df(1:50, get_ftable_bs) %>%
  left_join(df_nodes %>% dplyr::rename(feature_name = asv)) %>%
  dplyr::group_by(run) %>%
  dplyr::top_n(n = 500, wt = variable_importance) %>%
  ungroup() %>% select(name, com, feature_name, variable_importance) %>%
  group_by(name, com, feature_name) %>% dplyr::summarise(n = n(), mean_imp = mean(variable_importance), median_imp = median(variable_importance)) %>%
  mutate(prop = n/50) %>%
  filter(prop > 0.4)

osd2014_dada2_phyloseq_beta_filt <- osd2014_dada2_phyloseq_alpha
osd2014_dada2_phyloseq_beta_filt <- transform_sample_counts(osd2014_dada2_phyloseq_beta_filt, function(x)x/sum(x))
osd2014_dada2_phyloseq_beta_filt <- prune_samples(rownames(osd2014_dada2_phyloseq_beta_df), osd2014_dada2_phyloseq_beta_filt)
osd2014_dada2_phyloseq_beta_filt <- prune_taxa(colnames(osd2014_dada2_phyloseq_beta_df), osd2014_dada2_phyloseq_beta_filt)

osd2014_order_terrestrial <- osd2014_order_terrestrial %>%
  filter(label %in% halpern_impact$label)

all_comps <- bind_rows(
  get_g("com_1") %>%
    activate(nodes) %>%
    as_tibble() %>%
    inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
                 as_tibble(rownames = "label") %>%
                 gather(asv, prop, -label)) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)),
  get_g("com_2") %>%
    activate(nodes) %>%
    as_tibble() %>%
    inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
                 as_tibble(rownames = "label") %>%
                 gather(asv, prop, -label)) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)),
  get_g("com_3") %>%
    activate(nodes) %>%
    as_tibble() %>%
    inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
                 as_tibble(rownames = "label") %>%
                 gather(asv, prop, -label)) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)),
  get_g("com_4") %>%
    activate(nodes) %>%
    as_tibble() %>%
    inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
                 as_tibble(rownames = "label") %>%
                 gather(asv, prop, -label)) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)),
  get_g("com_5") %>%
    activate(nodes) %>%
    as_tibble() %>%
    inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
                 as_tibble(rownames = "label") %>%
                 gather(asv, prop, -label)) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)),
  get_g("com_7") %>%
    activate(nodes) %>%
    as_tibble() %>%
    inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
                 as_tibble(rownames = "label") %>%
                 gather(asv, prop, -label)) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)),
  get_g("com_8") %>%
    activate(nodes) %>%
    as_tibble() %>%
    inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
                 as_tibble(rownames = "label") %>%
                 gather(asv, prop, -label)) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)),
  get_g("com_9") %>%
    activate(nodes) %>%
    as_tibble() %>%
    inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
                 as_tibble(rownames = "label") %>%
                 gather(asv, prop, -label)) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial$label))
)  %>%
  group_by(Order) %>%
  mutate(agg_prop = sum(prop)) %>% ungroup() %>% mutate(order_mod = ifelse(agg_prop > 0.01, Order, "Other")) %>%
  inner_join(top_features %>% select(name, mean_imp, median_imp))


# all_comps_order <- all_comps %>% select(order_mod) %>% unique() %>%
#   mutate(colour = o_colors)

order_mod <- c("SAR11_clade", "Synechococcales", "Flavobacteriales", "Rhodobacterales", "Betaproteobacteriales", "Cellvibrionales",
               "Oceanospirillales", "Parvibaculales", "Micrococcales", "Microtrichales", "Puniceispirillales","Sphingobacteriales",  "SAR86_clade", "Other")
o_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
              "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
              "#cab2d6", "#D0D0D0", "#2E5158", "#6a3d9a","#F0D999", "#666666")
all_comps_order <- tibble(order_mod = order_mod, colour = o_colors)
# Get the graph with all components and save file for gephi ---------------
bind_graphs(
  get_g("com_1") %>%activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  get_g("com_2") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  get_g("com_3") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  get_g("com_4") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  get_g("com_5") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  #get_g("com_6")
  get_g("com_7") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  # get_g("com_8") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
  #   inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  get_g("com_9") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp))
) %>% write.graph(file = "osd2014_rf/data/osd2014_fuzzyforest_ohi.graphml", format = "graphml")



# Plot the abundance of each ASV in each LC and MEOW province -------------

# osd2014_impacted <- halpern_impact %>%
#   filter(class == "impacted")
# osd2014_impacted <- osd2014_impacted %>%
#   mutate(label = fct_relevel(label, osd2014_order_terrestrial %>% filter(label %in% osd2014_impacted$label) %>% arrange(position) %>% .$label))
#
# osd2014_low_impacted <- halpern_impact %>%
#   filter(class == "low impacted")
# osd2014_low_impacted <- osd2014_impacted %>%
#   mutate(label = fct_relevel(label, osd2014_order_terrestrial %>% filter(label %in% osd2014_low_impacted$label) %>% arrange(position) %>% .$label))

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

all_comps %>%
  filter(name %in% (top_features %>% arrange(desc(mean_imp))  %>% .$name)) %>%
  inner_join(halpern_impact) %>%
  select(label, Order, com, class, name, prop, mean_imp) %>% unique %>%
  filter(prop > 0) %>% #compare_means(formula = prop ~ class, group.by = "com", p.adjust.method = "fdr") %>% View
  ggplot(aes(class, (prop), fill = class)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.7, jitter.size = 1,
                        color = "black", alpha = 0.7, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4, errorbar.length = 0.2) +
  scale_y_continuous(trans = scales::log_trans(),  labels = scales::percent) +
  facet_wrap(~com, scales = "free",nrow = 2) +
  stat_compare_means(comparisons = list(c("impacted", "low_impacted")), aes(label=..p.adj..)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#2E3239", "#1F629A"))

ggsave(plot = last_plot(), filename = "osd2014_rf/figures/osd2014_fuzzyforests_bplots_ohi.pdf", width = 9, height = 4)

plot_asv <- function(X){
  all_comps %>%
    filter(name %in% (top_features %>% arrange(desc(mean_imp))  %>% .$name)) %>%
    inner_join(halpern_impact) %>%
    select(label, Order, com, class, name, prop, mean_imp) %>% unique %>%
    filter(prop > 0, name == X) %>% #compare_means(formula = prop ~ class, group.by = "com", p.adjust.method = "fdr") %>% View
    ggplot(aes(class, (prop), fill = class)) +
    ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.7, jitter.size = 1,
                          color = "black", alpha = 0.7, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4, errorbar.length = 0.2) +
    scale_y_continuous(trans = scales::log_trans(),  labels = scales::percent) +
    facet_wrap(~name, scales = "free",nrow = 2) +
    stat_compare_means(comparisons = list(c("impacted", "low_impacted")), aes(label=..p.adj..)) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#2E3239", "#1F629A"))
}
ggarrange(plot_asv("asv_218"),
          plot_asv("asv_23"),
          plot_asv("asv_86"),
          plot_asv("asv_22"),
          plot_asv("asv_1361"), ncol = 4, nrow = 2)

ggsave(plot = last_plot(), filename = "osd2014_rf/figures/osd2014_fuzzyforests_bplots_ohi_asv.pdf", width = 9, height = 4)


l<-all_comps %>%
  filter(name %in% (top_features %>% arrange(desc(mean_imp))  %>% .$name)) %>%
  inner_join(halpern_impact) %>%
  select(label, Order, com, class, name, prop, mean_imp) %>% unique %>%
  filter(prop > 0) %>% compare_means(data = ., prop ~ class, group.by = "com", p.adjust.method = "fdr")



all_comps %>%
  inner_join(halpern_impact) %>%
  group_by(label, Order, com, class) %>%
  dplyr::summarise(prop = sum(prop)) %>%
  ungroup() %>%
  filter(prop >0) %>%
  group_by(Order) %>%
  mutate(agg_prop = sum(prop)) %>%
  ungroup() %>%
  mutate(Order = ifelse((agg_prop > 0.01 & !is.na(Order)), Order, "Other")) %>%
  mutate(agg_prop = sum(agg_prop)) %>%
  ungroup() %>%
  #filter(!(com %in% c("com_8")))%>%#, meow_province %in% c("Mediterranean Sea", "Northern European Seas")) %>%
  #mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)) %>%
  #mutate(meow_province = fct_relevel(meow_province, meow_provinces)) %>%
  select(label, prop, Order, class, com) %>%
  #filter(meow_province == X) %>%
  droplevels() %>%
  complete(label,com, class, Order, fill = list(prop = 0)) %>%
  filter(!is.na(com)) %>%  #filter(meow_province == "Cold Temperate Northwest Atlantic") %>% View
  #inner_join(all_comps_order %>% dplyr::rename(Order = order_mod)) %>%
  ggplot(aes(label, prop, fill = Order)) +
  geom_col(color = "grey20", size = 0.2, width = 1) +
  facet_grid(class~com , scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  xlab("") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(limits = all_comps_order$order_mod, values = all_comps_order$colour)



plots <- map(meow_provinces, plot_com)
ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, ncol = 1, nrow = length(plots))
save.image(file = "osd2014_rf/data/osd2014_fuzzyforest_ohi.Rda")






#
#
# asv_table <- asv_otutable_filt_pa
# dim(asv_table)
#
# library(doParallel)
# library(caret)
# cl <- makeCluster(detectCores())
# registerDoParallel(cl)
# # Random Search
# seed <- 1234
# set.seed(seed)
#
# #length is = (n_repeats*nresampling)+1
# seeds <- vector(mode = "list", length = 50)
#
# #(3 is the number of tuning parameter, mtry for rf, here equal to ncol(iris)-2)
# for(i in 1:50) seeds[[i]] <- sample.int(n=1000, 15)
#
# #for the last model
# seeds[[51]] <- sample.int(1000, 1)
# nt <- 1000
# control <- trainControl(method="repeatedcv", number=10, repeats=5, search="random", seeds = seeds)
# asv_rf_random <- train(y = asv_table$class, x=asv_table[1:dim(asv_table)[2]-1], method="rf", metric="Accuracy", tuneLength=15,
#                        trControl=control, verboseIter = TRUE, savePredictions = TRUE,
#                        prox=TRUE, allowParallel=TRUE, num.trees=nt)
# stopCluster(cl)
# print(asv_rf_random)
# plot(asv_rf_random)
#
# #mtry = 329
#
# n <- 0.30
#
# repeat{
#   asv_halpern_classify <- randomForest::randomForest(y = asv_table$class,
#                                                      x = asv_table[1:dim(asv_table)[2]-1], ntree=nt, importance=TRUE, do.trace=nt, keep.inbag = TRUE, mtry = asv_rf_random$bestTune$mtry)
#   err <- asv_halpern_classify$err.rate[,1][nt]
#   diff <- abs(asv_halpern_classify$err.rate[,2][nt] - asv_halpern_classify$err.rate[,3][nt])
#   if(err < n && diff < 0.06){
#     break
#   }
# }
#
# #save(names_otus, rf_random, names_otus_short, halpern_classify, halpern_impact, asv_otutable, osd2014_16s_otuXsample_physeq_filt_prev_beta_mg, file = "~/ownCloud/OSD_paper/HALPERN/halpern_classify_RF_model.Rda")
# plot(asv_halpern_classify)
#
# print(asv_halpern_classify)
#
# asv_imp <- randomForest::importance(asv_halpern_classify)
# asv_imp <- data.frame(predictors = rownames(asv_imp), asv_imp)
#
# asv_imp.sort <- arrange(asv_imp, desc(MeanDecreaseGini))
# asv_imp.sort$predictors <- factor(asv_imp.sort$predictors, levels = asv_imp.sort$predictors)
#
# # Select the top 10 predictors
# asv_imp.sort <- asv_imp.sort[1:100, ]
#
# asv_imp.sort_s <- asv_imp.sort %>%
#   as_tibble() %>%
#   left_join(osd2014_silva_dada2_names %>% rename(predictors = asv)) %>%
#   mutate(asv_name = fct_reorder(asv_name, -MeanDecreaseGini))
# # ggplot
# ggplot(asv_imp.sort_s, aes(x = asv_name, y = MeanDecreaseGini)) +
#   geom_bar(stat = "identity", fill = "indianred") +
#   coord_flip() +
#   ggtitle("Most important OTUs for classifying into low or mid/hight impacted OSD samples") +
#   theme_light() +
#   ylab("Total decrease in node impurities (Mean decrease Gini)")
#
# asv_otutable_prop <- transform_sample_counts(osd2014_dada2_phyloseq_alpha, function(x) x/sum(x))
#
#
#
# f <- qpsmelt(asv_otutable_prop) %>%
#   tbl_df %>%
#   select(label, OTU, Abundance, asv_name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
#   dplyr::rename(asv = OTU) %>%
#   right_join(asv_imp.sort_s) %>%
#   filter(label %in% sites) %>%
#   unique() %>%
#   left_join(halpern_impact %>% select(label, median))
#
#
# f <- f %>%
#   mutate(label = fct_reorder(label, -median), asv_name = fct_reorder(asv_name, MeanDecreaseGini))
#
#
# base_breaks <- function(n = 10){
#   function(x) {
#     axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
#   }
# }
#
# textcol <- "grey40"
#
# #modified ggplot
# p_asv <- ggplot(f,aes(y=label,x=asv_name, fill = ((Abundance) +0.001)*100))+
#   geom_tile()+
#   #redrawing tiles to remove cross lines from legend
#   geom_tile(colour="white",size=0.25, show.legend = FALSE)+
#   #remove axis labels, add title
#   labs(x="",y="",title="")+
#   #remove extra space
#   scale_y_discrete(expand=c(0,0))+
#   #custom breaks on x-axis
#   scale_x_discrete(expand=c(0,0))+
#   #custom colours for cut levels and na values
#   #scale_fill_gradientn(colours =rev(c("#d53e4f","#f46d43","#fdae61",
#   # "#fee08b","#e6f598","#abdda4","#ddf1da")), na.value="grey90", trans = "sqrt", labels = scales::percent) +
#   viridis::scale_fill_viridis(option = "D", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
#   #scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(7,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
#   #mark year of vaccination
#   #geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
#   #equal aspect ratio x and y axis
#   coord_fixed() +
#   #set base size for all font elements
#   theme_grey(base_size=10)+
#   #theme options
#   theme(
#     legend.position = "bottom",
#     #remove legend title
#     legend.title=element_blank(),
#     #remove legend margin
#     legend.spacing = grid::unit(0,"cm"),
#     #change legend text properties
#     legend.text=element_text(colour=textcol,size=7,face="bold"),
#     #change legend key height
#     legend.key.height=grid::unit(0.2,"cm"),
#     #set a slim legend
#     legend.key.width=grid::unit(0.8,"cm"),
#     #set x axis text size and colour
#     axis.text.x=element_text(hjust = 1, vjust = 0.5, colour=textcol, angle = 90, size = 6),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     #set y axis text colour and adjust vertical justification
#     #axis.text.y=element_text(vjust = 0.2,colour=textcol, size = 6),
#     #change axis ticks thickness
#     axis.ticks=element_line(size=0.4),
#     #change title font, size, colour and justification
#     plot.title=element_blank(),
#     #remove plot background
#     plot.background=element_blank(),
#     #remove plot border
#     panel.border=element_blank())
#
# f1 <- f %>% select(label, median) %>% unique
#
# f1$label <- factor(f1$label, levels = f1 %>% arrange(median) %>% .$label)
#
# p1 <- ggplot( f1, aes(x = label, y = median)) +
#   geom_bar(stat = "identity", size = 0.25) +
#   coord_flip() +
#   #remove axis labels, add title
#   labs(x="",y="",title="") +
#   #remove extra space
#   #custom breaks on x-axis
#   scale_x_discrete(expand=c(0,0)) +
#   #custom colours for cut levels and na values
#   #scale_fill_gradientn(colours =rev(c("#d53e4f","#f46d43","#fdae61",
#   # "#fee08b","#e6f598","#abdda4","#ddf1da")), na.value="grey90", trans = "sqrt", labels = scales::percent) +
#   #scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(7,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
#   #mark year of vaccination
#   #geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
#   #equal aspect ratio x and y axis
#   #coord_fixed() +
#   #set base size for all font elements
#   theme_grey(base_size=10)+
#   #theme options
#   theme(
#     legend.position = "bottom",
#     #remove legend title
#     legend.title=element_blank(),
#     #remove legend margin
#     legend.spacing = grid::unit(0,"cm"),
#     #change legend text properties
#     legend.text=element_text(colour=textcol,size=7,face="bold"),
#     #change legend key height
#     legend.key.height=grid::unit(0.2,"cm"),
#     #set a slim legend
#     legend.key.width=grid::unit(0.8,"cm"),
#     #set x axis text size and colour
#     axis.text.y=element_blank(),
#     axis.ticks.y = element_blank(),
#     #set y axis text colour and adjust vertical justification
#     axis.text.x=element_text(colour=textcol),
#     #change axis ticks thickness
#     axis.ticks=element_line(size=0.4),
#     #change title font, size, colour and justification
#     plot.title=element_blank(),
#     #remove plot background
#     plot.background=element_blank(),
#     #remove plot border
#     panel.border=element_blank(),
#     panel.background = element_blank())
#
# f2 <- f %>% select(asv_name, MeanDecreaseGini) %>% unique
#
# f2$asv_name <- factor(f2$asv_name, levels = f2 %>% arrange(MeanDecreaseGini) %>% .$asv_name)
#
# p2 <- ggplot(f2, aes(x = asv_name, y = MeanDecreaseGini)) +
#   geom_bar(stat = "identity", size = 0.25) +
#   #remove axis labels, add title
#   labs(x="",y="",title="")+
#   #remove extra space
#   #custom breaks on x-axis
#   scale_x_discrete(expand=c(0,0))+
#   #custom colours for cut levels and na values
#   #scale_fill_gradientn(colours =rev(c("#d53e4f","#f46d43","#fdae61",
#   # "#fee08b","#e6f598","#abdda4","#ddf1da")), na.value="grey90", trans = "sqrt", labels = scales::percent) +
#   #scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(7,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
#   #mark year of vaccination
#   #geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
#   #equal aspect ratio x and y axis
#   #set base size for all font elements
#   theme_grey(base_size=10)+
#   #theme options
#   theme(
#     legend.position = "bottom",
#     #remove legend title
#     legend.title=element_blank(),
#     #remove legend margin
#     legend.spacing = grid::unit(0,"cm"),
#     #change legend text properties
#     legend.text=element_text(colour=textcol,size=7,face="bold"),
#     #change legend key height
#     legend.key.height=grid::unit(0.2,"cm"),
#     #set a slim legend
#     legend.key.width=grid::unit(0.8,"cm"),
#     #set x axis text size and colour
#     axis.text.x=element_blank(),
#     axis.ticks.x = element_blank(),
#     #set y axis text colour and adjust vertical justification
#     axis.text.y=element_text(vjust = 0.2,colour=textcol),
#     #change axis ticks thickness
#     axis.ticks=element_line(size=0.4),
#     #change title font, size, colour and justification
#     plot.title=element_blank(),
#     #remove plot background
#     plot.background=element_blank(),
#     #remove plot border
#     panel.border=element_blank(),
#     panel.background = element_blank())
#
# ggpubr::ggarrange(p_asv, p1, p2, ncol = 1, nrow = 3, heights = c(0.5, .25, .25))
#
#
#
#
#
# k <- f %>%
#   filter(asv %in% asv_imp.sort$predictors) %>%
#   left_join(halpern_impact) %>%
#   rename(type = class) %>%
#   select(type, Abundance, asv_name, MeanDecreaseAccuracy, asv)
#
# asv_wt <- ggpubr::compare_means(Abundance ~ type, data = k, group.by = c("asv_name"), method = "wilcox.test", p.adjust.method = "BH") %>%
#   filter(p.signif != "ns") %>% arrange(p.adj) %>% left_join(k)
#
# asv_wt %>%
#   arrange(desc(MeanDecreaseAccuracy)) %>%
#   filter(p.signif != "ns") %>%
#   select(asv_name, MeanDecreaseAccuracy) %>% unique()
#
# asv_imp.sort_s[1:25,]
#
# asv_wt$type <- factor(asv_wt$type, levels = c("low_impacted", "impacted"))
#
# asv_wt$asv <- factor(asv_wt$asv, levels = asv_imp.sort_s$predictors)
#
#
# p <- ggplot(asv_wt %>% filter(asv_name == "asv_33"), aes(type, Abundance, fill = type, color = type)) +
#   geom_boxplot(position = "dodge", width = 0.5, alpha = 0.6) +
#   geom_jitter(alpha = 0.6, color = "#4A4A4A", width = 0.25, ) +
#   ggpubr::stat_compare_means() +
#   scale_y_continuous(labels = scales::percent) +
#   theme_bw() +
#   xlab("") +
#   scale_fill_manual(values = c("#2E9598", "#A8216B")) +
#   scale_color_manual(values = c("#2E9598", "#A8216B")) +
#   theme(legend.position = "none")
#
# p1 <- ggplot(asv_wt %>% filter(asv_name == "ASV33"), aes(type, value, fill = type, color = type)) +
#   geom_boxplot(position = "dodge", width = 0.5, alpha = 0.6) +
#   geom_jitter(alpha = 0.6, color = "#4A4A4A", width = 0.25, ) +
#   ggpubr::stat_compare_means() +
#   scale_y_continuous(labels = scales::percent) +
#   theme_bw() +
#   xlab("") +
#   scale_fill_manual(values = c("#2E9598", "#A8216B")) +
#   scale_color_manual(values = c("#2E9598", "#A8216B")) +
#   theme(legend.position = "none")
#
#
# p2 <- ggplot(asv_wt %>% filter(asv_s == "ASV1543"), aes(type, value, fill = type, color = type)) +
#   geom_boxplot(position = "dodge", width = 0.5, alpha = 0.6) +
#   geom_jitter(alpha = 0.6, color = "#4A4A4A", width = 0.25, ) +
#   ggpubr::stat_compare_means() +
#   scale_y_continuous(labels = c("0.0%", "1.0%", "2.0%", "3.0%", "4.0%"), breaks = c(0, 0.01, 0.02, 0.03, 0.04)) +
#   #scale_y_continuous(labels = scales::percent) +
#   theme_bw() +
#   xlab("") +
#   scale_fill_manual(values = c("#2E9598", "#A8216B")) +
#   scale_color_manual(values = c("#2E9598", "#A8216B")) +
#   theme(legend.position = "none")
#
# p3 <- ggplot(asv_wt %>% filter(asv_s == "ASV1463"), aes(type, value, fill = type, color = type)) +
#   geom_boxplot(position = "dodge", width = 0.5, alpha = 0.6) +
#   geom_jitter(alpha = 0.6, color = "#4A4A4A", width = 0.25, ) +
#   ggpubr::stat_compare_means(label = "p.format") +
#   scale_y_continuous(labels = scales::percent) +
#   theme_bw() +
#   xlab("") +
#   scale_fill_manual(values = c("#2E9598", "#A8216B")) +
#   scale_color_manual(values = c("#2E9598", "#A8216B")) +
#   theme(legend.position = "none")
#
# ggpubr::ggarrange(p3,p2,p1,p, ncol = 4, nrow = 1)
# ggsave("~/Downloads/asv_bp.pdf", width = 29.39, height = 6.92, units = "cm")
#
#
# tax_table(osd2014_16s_otuXsample_physeq_filt_prev_prop) %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "old") %>%
#   tbl_df() %>%
#   left_join(asv_name_conv) %>%
#   filter(new %in% c(as.character(asv_imp.sort_s$predictors[1:3]), as.character(asv_imp.sort_s$predictors[9]))) %>%
#   select(-old, -new, -genus)
#
# # Let's predict
#
# # First get samples not used for training
# pred_sites <- out_rescaled_2013_median_long %>% filter(!(label %in% sites)) %>% .$label %>% unique
#
#
# # Trim the original matrix to contain the non-training and the modules used for
# asv_otutable_pred <- as.data.frame(as(otu_table(osd2014_16s_otuXsample_physeq_filt_prev_prop), "matrix"))
# asv_otutable_pred <- base::as.data.frame(asv_otutable_pred)
# asv_otutable_pred <- asv_otutable_pred[pred_sites, c(as.character(asv_shared_impact_low_occ$asv), asv_impact_low$asv, asv_impact_imp$asv)]
#
# pred_asv <- predict(halpern_classify, newdata = asv_otutable_pred)
#
# pred_impacted_metadata_asv <- data.frame(class = pred_asv) %>% mutate(label = names(pred_asv)) %>% filter(class == "impacted") %>% left_join(osd2014_metadata)
# pred_low_impacted_metadata_asv <- data.frame(class = pred_asv) %>% mutate(label = names(pred_asv)) %>% filter(class != "impacted") %>% left_join(osd2014_metadata)
#
# library(ggmap)
#
# ### Set a range
#
# ### Get a map
#
# tmp <- pred_impacted_metadata %>% dplyr::select(label, start_lat,start_lon, dist_coast_iso3_code) %>%
#   mutate(filename = paste("~/Downloads/", label, "_", dist_coast_iso3_code, "_asv_impacted.pdf", sep = ""))
#
#
# for (i in 1:dim(tmp)[1]) {
#
#   map <- get_map(location = c(lon = tmp[i,]$start_lon, lat = tmp[i,]$start_lat), zoom = 13,
#                  maptype = "satellite", source = "google")
#
#   ### When you draw a figure, you limit lon and lat.
#   foo <- ggmap(map) +
#     geom_point(data = tmp[i,], aes(x = start_lon, y = start_lat), size = 3, shape = 21, fill = "red", alpha = 0.8)
#   pdf(file = tmp[i,]$filename)
#   print(foo)
#   dev.off()
# }
#
#
# tmp <- pred_low_impacted_metadata %>% dplyr::select(label, start_lat,start_lon, dist_coast_iso3_code) %>%
#   mutate(filename = paste("~/Downloads/", label, "_", dist_coast_iso3_code, "_asv_low_impacted.pdf", sep = ""))
#
#
# for (i in 1:dim(tmp)[1]) {
#
#   map <- get_map(location = c(lon = tmp[i,]$start_lon, lat = tmp[i,]$start_lat), zoom = 13,
#                  maptype = "satellite", source = "google")
#
#   ### When you draw a figure, you limit lon and lat.
#   foo <- ggmap(map) +
#     geom_point(data = tmp[i,], aes(x = start_lon, y = start_lat), size = 3, shape = 21, fill = "red", alpha = 0.8)
#   pdf(file = tmp[i,]$filename)
#   print(foo)
#   dev.off()
# }
#
#
# pred_impacted_metadata_kegg_mg <- pred_impacted_metadata_kegg %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label
# pred_impacted_metadata_asv_mg <- pred_impacted_metadata_asv %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label
#
# pred_impacted_metadata_kegg_mg_l <- length(pred_impacted_metadata_kegg_mg)
# pred_impacted_metadata_asv_mg_l <- length(pred_impacted_metadata_asv_mg)
#
# length(intersect(pred_impacted_metadata_kegg_mg,pred_impacted_metadata_asv_mg))/(pred_impacted_metadata_kegg_mg_l + pred_impacted_metadata_asv_mg_l)
#
# pred_low_impacted_metadata_kegg_mg <- pred_low_impacted_metadata_kegg %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label
# pred_low_impacted_metadata_asv_mg <- pred_low_impacted_metadata_asv %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label
#
# pred_low_impacted_metadata_kegg_mg_l <- length(pred_low_impacted_metadata_kegg_mg)
# pred_low_impacted_metadata_asv_mg_l <- length(pred_low_impacted_metadata_asv_mg)
#
# length(intersect(pred_low_impacted_metadata_kegg_mg,pred_low_impacted_metadata_asv_mg))/(pred_low_impacted_metadata_kegg_mg_l + pred_low_impacted_metadata_asv_mg_l)
#
#
# length(as.character(osd2014_amp_mg_intersect$label)[as.character(pred_impacted_metadata_kegg$label)])
# length(pred_impacted_metadata_asv$label)
# length(pred_low_impacted_metadata_kegg$label)
# length(pred_low_impacted_metadata_asv$label)
#
# intersect(pred_impacted_metadata_kegg$label, pred_impacted_metadata_asv$label)
# intersect(pred_low_impacted_metadata_kegg$label, pred_low_impacted_metadata_asv$label)
# length(intersect(pred_impacted_metadata_kegg$label, pred_impacted_metadata_asv$label))
# length(intersect(pred_low_impacted_metadata_kegg$label, pred_low_impacted_metadata_asv$label))
#


