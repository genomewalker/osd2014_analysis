library(fuzzyforest)
library(tidyverse)
library(randomForest)
library(caret)
library(phyloseq)
library(igraph)
library(tidygraph)
library(ggraph)
source("osd2014_16S_asv/lib/fuzzyforest_lib.R")

load("osd2014_16S_asv/data/osd2014_16S_asv_networks_results.Rdata", verbose = TRUE)
load("~/Downloads/ff_fit.Rda")
#set seed so that results are reproducible
set.seed(1)

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)
osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)
osd2014_cdata <- osd2014_cdata %>%
  filter(label %in% osd2014_amp_mg_intersect$label, meow_province %in% osd2014_meow_regions$meow_province)

osd2014_asv_connectedness <- tbl(my_db, "osd2014_asv_connectedness") %>%
  collect(n = Inf)

osd2014_order_terrestrial <- tbl(my_db, "osd2014_st_order_terrestrial") %>%
  collect(n = Inf)
osd2014_dada2_phyloseq_alpha_prop <- transform_sample_counts(osd2014_dada2_phyloseq_alpha, function(x) x/sum(x))
osd2014_order_terrestrial <- osd2014_order_terrestrial %>%
  filter(label %in% osd2014_cdata$label)

osd2014_dada2_phyloseq_alpha_filt <- subset_samples(osd2014_dada2_phyloseq_alpha_prop, label %in% osd2014_cdata$label)
osd2014_dada2_phyloseq_alpha_filt <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_alpha_filt) > 0, osd2014_dada2_phyloseq_alpha_filt)
osd2014_dada2_phyloseq_beta_filt <- prune_taxa(taxa_names(osd2014_dada2_phyloseq_beta), osd2014_dada2_phyloseq_alpha_filt)


#osd2014_dada2_phyloseq_beta_filt <- subset_samples(osd2014_dada2_phyloseq_beta_vst, label %in% osd2014_cdata$label)
#osd2014_dada2_phyloseq_beta_filt <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_beta_filt) > 0, osd2014_dada2_phyloseq_beta_filt)
osd2014_dada2_phyloseq_beta_df <- (as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix")) %>% as_tibble(rownames = "label") %>%
  inner_join(osd2014_cdata %>%
               select(label, meow_province)) %>%
  as.data.frame() %>%
  column_to_rownames("label")

osd2014_dada2_phyloseq_beta_df$meow_province <- as.factor(osd2014_dada2_phyloseq_beta_df$meow_province)


# Run fuzzforest with 5 repetitions and 10 k-folds ------------------------

ff_training <- vector(mode = "list")
ff_fit_r <- vector(mode = "list")
ff_predict <- vector(mode = "list")
ff_accuracy <- vector(mode = "list")

for (i in 1:5){
  training.samples <- osd2014_dada2_phyloseq_beta_df$meow_province %>%
    createFolds(returnTrain = TRUE)
  #createDataPartition(p = 0.8, list = FALSE)
  ff_training[[i]] <- training.samples
  for (o in 1:10){
    cat(paste("Rep:", i, "Fold:", o, "\n"))
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
    train.data$meow_province <- as.factor(train.data$meow_province)
    ff_fit <- ff(train.data[,1:ncol(train.data) - 1], train.data[,ncol(train.data)], module_membership = module_membership,
                 screen_params = screen_params, select_params=select_params,
                 final_ntree = 500, num_processors = 4)

    ff_fit_r[[paste0("iter_",i,"_fold_",o)]] <- ff_fit
    pred_asv <- predict(ff_fit$final_rf, newdata = test.data[,1:ncol(test.data) - 1])
    ff_predict[[paste0("iter_",i,"_fold_",o)]] <- pred_asv

    test.data$rightPred <- pred_asv == test.data$meow_province
    accuracy <- sum(test.data$rightPred)/nrow(test.data)
    ff_accuracy[[paste0("iter_",i,"_fold_",o)]] <- accuracy
  }
}



# Check how accurate are the models ---------------------------------------

bind_rows(map(1:50, get_accuracy)) %>% ggplot(aes(y = value, x = "accuracy")) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 1,
                        color = "black", alpha = 1, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4, errorbar.length = 0.2) +
  geom_dotplot(binaxis = "y", dotsize = 0.2, stackdir = "down", binwidth = 0.1,
               position = position_nudge(-0.025)) +
   theme_bw() +
  ylab("Accuracy") +
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



# Define some colors for phylum -------------------------------------------
# p_colors <- tibble(Phylum = c("Proteobacteria", "Bacteroidetes", "Cyanobacteria", "Actinobacteria",
#                               "Verrucomicrobia", "Planctomycetes", "Euryarchaeota", "Marinimicrobia_(SAR406_clade)",
#                               "Firmicutes", "Lentisphaerae",  "Tenericutes", "Chlamydiae", "Epsilonbacteraeota",
#                               "Patescibacteria", "Other"),
#                    colour = c("#cf4149", "#8f62ca", "#66b14a", "#6b8bcd", "#c2ad4b", "#ce7531",
#                               "#4baf90", "#c85d9d", "#542437", "#c26f65", "#A5C990", "#767834",
#                               "#559279", "#D3BBC3", "#7f8c8d"))
#
#



# Get all components ------------------------------------------------------

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
  # get_g("com_3") %>%
  #   activate(nodes) %>%
  #   as_tibble() %>%
  #   inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
  #                as_tibble(rownames = "label") %>%
  #                gather(asv, prop, -label)) %>%
  #   mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)),
  get_g("com_4") %>%
    activate(nodes) %>%
    as_tibble() %>%
    inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
                 as_tibble(rownames = "label") %>%
                 gather(asv, prop, -label)) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)),
  # get_g("com_5") %>%
  #   activate(nodes) %>%
  #   as_tibble() %>%
  #   inner_join(as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>%
  #                as_tibble(rownames = "label") %>%
  #                gather(asv, prop, -label)) %>%
  #   mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)),
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
  mutate(agg_prop = sum(prop)) %>% ungroup() %>% mutate(order_mod = ifelse(agg_prop > 0.01, Order, "Other"))





# all_comps_order <- all_comps %>% select(order_mod) %>% unique() %>%
#   mutate(colour = o_colors)




order_mod <- c("SAR11_clade", "Synechococcales", "Flavobacteriales", "Rhodobacterales", "Verrucomicrobiales", "Betaproteobacteriales",
               "Oceanospirillales", "Cellvibrionales", "Micrococcales", "Rhodospirillales", "KI89A_clade", "Cytophagales ",
               "Arenicellales", "Chitinophagales", "Tenderiales", "Other")
o_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
              "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
              "#cab2d6","#6a3d9a","#D0D0D0", "#2E5158", "#6A6D51",
              "#F0D999", "#43233A","#666666")
all_comps_order <- tibble(order_mod = order_mod, colour = o_colors)
# Get the graph with all components and save file for gephi ---------------
bind_graphs(
  get_g("com_1") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  get_g("com_2") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  #get_g("com_3"),
  get_g("com_4") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  #get_g("com_5"),
  #get_g("com_6")
  get_g("com_7") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp)),
  #get_g("com_8") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>% inner_join(all_comps_order),
  get_g("com_9") %>% activate(nodes) %>% mutate(order_mod = ifelse(Order %in% all_comps_order$order_mod, Order, "Other")) %>%
    inner_join(all_comps_order) %>% inner_join(top_features %>% select(name, mean_imp, median_imp))
) %>% write.graph(file = "osd2014_16S_asv/data/osd2014_fuzzyforest.graphml", format = "graphml")






# Plot the abundance of each ASV in each LC and MEOW province -------------

meow_provinces <- c("Tropical Northwestern Atlantic", "Warm Temperate Northwest Atlantic", "Cold Temperate Northwest Atlantic",
                    "Lusitanian", "Mediterranean Sea", "Northern European Seas")

plots <- map(meow_provinces, plot_com)
ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, ncol = 1, nrow = length(plots))
ggsave(plot = last_plot(), filename = "osd2014_16S_asv/figures/osd2014_fuzzyforests_bplots.pdf", width = 11.69, height = 8.27)


# Test with only the selected ASVs ----------------------------------------

top_features$feature_name



osd2014_dada2_phyloseq_alpha_prop <- transform_sample_counts(osd2014_dada2_phyloseq_alpha, function(x) x/sum(x))
osd2014_order_terrestrial <- osd2014_order_terrestrial %>%
  filter(label %in% osd2014_cdata$label)

osd2014_dada2_phyloseq_alpha_filt <- subset_samples(osd2014_dada2_phyloseq_alpha_prop, label %in% osd2014_cdata$label)
osd2014_dada2_phyloseq_beta_filt <- prune_taxa(top_features$feature_name, osd2014_dada2_phyloseq_alpha_filt)
osd2014_dada2_phyloseq_beta_filt <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_beta_filt) > 0, osd2014_dada2_phyloseq_beta_filt)


#osd2014_dada2_phyloseq_beta_filt <- subset_samples(osd2014_dada2_phyloseq_beta_vst, label %in% osd2014_cdata$label)
osd2014_dada2_phyloseq_beta_df <- (as(otu_table(osd2014_dada2_phyloseq_beta_filt), "matrix")) %>% as_tibble(rownames = "label") %>%
  inner_join(osd2014_cdata %>%
               select(label, meow_province)) %>%
  as.data.frame() %>%
  column_to_rownames("label")

osd2014_dada2_phyloseq_beta_df$meow_province <- as.factor(osd2014_dada2_phyloseq_beta_df$meow_province)


ff_training_test <- vector(mode = "list")
ff_fit_r_test <- vector(mode = "list")
ff_predict_test <- vector(mode = "list")
ff_accuracy_test <- vector(mode = "list")


for (i in 1:5){
  training.samples <- osd2014_dada2_phyloseq_beta_df$meow_province %>%
    createFolds(returnTrain = TRUE)
  #createDataPartition(p = 0.8, list = FALSE)
  ff_training_test[[i]] <- training.samples
  for (o in 1:10){
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
    train.data$meow_province <- as.factor(train.data$meow_province)
    ff_fit <- ff(train.data[,1:ncol(train.data) - 1], train.data[,ncol(train.data)], module_membership = module_membership,
                 screen_params = screen_params, select_params=select_params,
                 final_ntree = 500, num_processors = 4)

    ff_fit_r_test[[paste0("iter_",i,"_fold_",o)]] <- ff_fit

    #
    # asv_imp <- randomForest::importance(ff_fit$final_rf, type = 1, scale = F)
    # asv_imp <- data.frame(predictors = rownames(asv_imp), asv_imp) %>% as_tibble()
    #
    # asv_imp.sort <- arrange(asv_imp, desc(MeanDecreaseAccuracy))
    # asv_imp.sort$predictors <- factor(asv_imp.sort$predictors, levels = asv_imp.sort$predictors)
    #
    # # Select the top 10 predictors
    # asv_imp.sort <- asv_imp.sort[1:100, ]
    #
    # asv_imp.sort_s <- asv_imp.sort %>%
    #   as_tibble() %>%
    #   left_join(df_nodes %>% dplyr::rename(predictors = asv)) %>%
    #   mutate(name = fct_reorder(name, -MeanDecreaseAccuracy))
    # # ggplot
    # ggplot(asv_imp.sort_s, aes(x = name, y = MeanDecreaseAccuracy)) +
    #   geom_bar(stat = "identity", fill = "indianred") +
    #   coord_flip() +
    #   ggtitle("Most important OTUs for classifying into low or mid/hight impacted OSD samples") +
    #   theme_light() +
    #   ylab("Total decrease in node impurities (Mean decrease Gini)")

    pred_asv <- predict(ff_fit$final_rf, newdata = test.data[,1:ncol(test.data) - 1])
    ff_predict_test[[paste0("iter_",i,"_fold_",o)]] <- pred_asv

    test.data$rightPred <- pred_asv == test.data$meow_province
    accuracy <- sum(test.data$rightPred)/nrow(test.data)
    ff_accuracy_test[[paste0("iter_",i,"_fold_",o)]] <- accuracy
  }
}

get_accuracy_test <- function(X){
  ff_accuracy_test[[X]] %>% as_tibble() %>% mutate(run = names(ff_accuracy_test[X]))
}
bind_rows(map(1:50, get_accuracy)) %>% ggplot(aes(y = value, x = "accuracy")) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 1,
                        color = "black", alpha = 1, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4, errorbar.length = 0.2) +
  theme_bw() +
  ylab("Accuracy") +
  xlab("")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

save.image("osd2014_16S_asv/data/osd2014_fuzzyforest.Rda")
