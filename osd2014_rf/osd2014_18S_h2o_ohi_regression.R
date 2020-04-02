library(h2o)
library(DESeq2)
library(tidyverse)
library(phyloseq)
library(fuzzyforest)
library(igraph)
library(tidygraph)
library(caret)
library(ggpubr)
library(propr)
library(ape)
source("osd2014_16S_asv/lib/fuzzyforest_lib.R")
source("osd2014_rf/lib/h2o_model.R")
source("osd2014_rf/lib/philr_var.R")

load("osd2014_18S_asv/data/osd2014_18S_asv_physeq_filt_objects_with_phylo.Rdata")
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

skimr::skim(out_rescaled_2013_median_long %>% group_by(ohi_variable, buffer) %>% filter(ohi_variable == "global_cumul_impact", median > 0))

p <- ggplot(out_rescaled_2013_median_long %>% filter(ohi_variable == "global_cumul_impact", median > 0), aes(median)) +
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

# Select those samples that are impacted and at least in 2km from the shore
halpern_impact <- out_rescaled_2013_median_long %>%
  filter(buffer == "1km", ohi_variable == "global_cumul_impact", median > 0) %>%
  filter(median > 0, median <= 8.47) %>%
  left_join(osd2014_cdata) %>%
  filter(label %in% osd2014_amp_mg_intersect$label)


asv_names <- as(tax_table(osd2014_dada2_phyloseq_alpha), "matrix") %>% as_tibble(rownames = "asv") %>% select(asv, asv_name)

# Filter samples that have an OHI value
osd2014_dada2_phyloseq_alpha_filt <- subset_samples(osd2014_dada2_phyloseq_alpha, label %in% halpern_impact$label)
# remove empty ASVs
osd2014_dada2_phyloseq_alpha_filt <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_alpha_filt) > 0, osd2014_dada2_phyloseq_alpha_filt)

# Filter ASVs that doesn't have more than 3 counts in 10% of the samples
osd2014_dada2_phyloseq_alpha_filt <- filter_taxa(osd2014_dada2_phyloseq_alpha_filt, function(x) sum(x > 3) > (0.1*length(x)), TRUE)

# Explore the coeficient of variation
osd2014_dada2_phyloseq_alpha_prop <- transform_sample_counts(osd2014_dada2_phyloseq_alpha_filt, function(x) x/sum(x))
tax.mean <- taxa_sums(osd2014_dada2_phyloseq_alpha_prop)/nsamples(osd2014_dada2_phyloseq_alpha_prop)

tax.cv <- base::apply(as(otu_table(osd2014_dada2_phyloseq_alpha_prop), "matrix"),2, function(x){sd(x)/mean(x)})
tax.cv.c <- base::apply(as(otu_table(osd2014_dada2_phyloseq_alpha_prop), "matrix"),2, function(X)DescTools::CoefVar(X, unbiased = TRUE))

tibble(asv = names(tax.mean), mean = tax.mean) %>%
  left_join(tibble(asv = names(tax.mean), cv = tax.cv)) %>%
  ggplot(aes(mean, cv)) +
  geom_point()

tibble(asv = names(tax.cv.c), cv.c = tax.cv.c) %>%
  left_join(tibble(asv = names(tax.mean), cv = tax.cv)) %>%
  ggplot(aes(cv, cv.c)) +
  geom_point()

tibble(asv = names(tax.mean), mean = tax.mean) %>%
  left_join(tibble(asv = names(tax.cv.c), cv.c = tax.cv.c)) %>%
  ggplot(aes(mean, cv.c)) +
  geom_point()


# the sequencing depth of each sample was standardized to the abundance of the median sampling depth,
total = median(sample_sums(osd2014_dada2_phyloseq_alpha_filt))
standf = function(x, t=total) round(t * (x / sum(x)))
osd2014_dada2_phyloseq_alpha_filt = transform_sample_counts(osd2014_dada2_phyloseq_alpha_filt, standf)

# Filter ASVs with a CV <= 3.0
#osd2014_dada2_phyloseq_alpha_norare = filter_taxa(osd2014_dada2_phyloseq_alpha_filt, function(x) sd(x)/mean(x) > 3.0, TRUE)
osd2014_dada2_phyloseq_alpha_norare <- filter_taxa(osd2014_dada2_phyloseq_alpha_filt, function(X)DescTools::CoefVar(X, unbiased = TRUE) > 3.0, TRUE)

# phi <- propr::propr((as(otu_table(osd2014_dada2_phyloseq_alpha_norare), "matrix")), metric = "rho")
# hc <- caret::findCorrelation(phi@matrix, cutoff = 0.9)
# hc <- sort(hc)
# reduced_Data <- colnames(phi@matrix)[c(hc)]
#
#
# osd2014_dada2_phyloseq_alpha_norare <- microbiome::remove_taxa(reduced_Data, osd2014_dada2_phyloseq_alpha_norare)
#osd2014_dada2_phyloseq_alpha_norare <- microbiome::transform(osd2014_dada2_phyloseq_alpha_norare, transform = "compositional")
#osd2014_dada2_phyloseq_alpha_norare <- microbiome::transform(osd2014_dada2_phyloseq_alpha_norare, transform = "log10p")

# pseudocount of 1 was added prior to PhILR transformation to avoid taking log-ratios with  zero counts.
osd2014_dada2_phyloseq_alpha_norare <- transform_sample_counts(osd2014_dada2_phyloseq_alpha_norare, function(x) x + 1)

(otus.outgroup <- row.names(subset(tax_table(osd2014_dada2_phyloseq_alpha_norare),tax_table(osd2014_dada2_phyloseq_alpha_norare)[,'Phylum']=="Euryarchaeota")))
phy_tree(osd2014_dada2_phyloseq_alpha_norare) <- root(phy_tree(osd2014_dada2_phyloseq_alpha_norare), otus.outgroup,resolve.root = TRUE)
tax_table(osd2014_dada2_phyloseq_alpha_norare) <- tax_table(osd2014_dada2_phyloseq_alpha_norare)[,-ncol(tax_table(osd2014_dada2_phyloseq_alpha_norare))]

is.rooted(phy_tree(osd2014_dada2_phyloseq_alpha_norare))
is.binary.tree(phy_tree(osd2014_dada2_phyloseq_alpha_norare))
phy_tree(osd2014_dada2_phyloseq_alpha_norare) <- makeNodeLabel(phy_tree(osd2014_dada2_phyloseq_alpha_norare), method="number", prefix='n')
name.balance(phy_tree(osd2014_dada2_phyloseq_alpha_norare), tax_table(osd2014_dada2_phyloseq_alpha_norare), 'n1')
otu.table <- (as(otu_table(osd2014_dada2_phyloseq_alpha_norare), "matrix"))
tree <- phy_tree(osd2014_dada2_phyloseq_alpha_norare)
metadata <- sample_data(osd2014_dada2_phyloseq_alpha_norare)
tax <- tax_table(osd2014_dada2_phyloseq_alpha_norare)

gp.philr <- philr(otu.table, tree,
                  part.weights='enorm.x.gm.counts',
                  ilr.weights='blw.sqrt')


# philr -------------------------------------------------------------------

osd2014_dada2_phyloseq_alpha_norare


osd2014_dada2_phyloseq_beta_df <- gp.philr %>% as_tibble(rownames = "label") %>%
  inner_join(halpern_impact %>%
               select(label, median)) %>%
  mutate(class = ohi_median2class_short(median)) %>%
  as.data.frame() %>%
  column_to_rownames("label")

remote_h2o <- h2o.init(nthreads = 64, ip = "localhost", port = 55577)
#remote_h2o <- h2o.init(nthreads = 4, ip = "localhost", port = 55599)#, port = 55577)

#df <- as.h2o(osd2014_dada2_phyloseq_beta_df)
df <- osd2014_dada2_phyloseq_beta_df


#seed <- round(as.numeric(as.POSIXct(Sys.time())))
#seed <- seeds[[o]]
#seed_c <- round(seeds[[o]])
#cat("Seed:", seed)

create_df <- function(X, physeq){
  # physeq <- osd2014_dada2_phyloseq_alpha_norare
  # X <- df_folds[[1]]
  train_fold <- sample_names(physeq)[X]
  test_fold <-  sample_names(physeq)[-X]

  physeq_train <- prune_samples(samples = train_fold, physeq)
  physeq_test <- prune_samples(samples = test_fold, physeq)

  otu_table_train <- (as(otu_table(physeq_train), "matrix"))
  otu_table_test <- (as(otu_table(physeq_test), "matrix"))

  philr_train <- philr(otu_table_train, phy_tree(physeq_train),
                       part.weights='enorm.x.gm.counts',
                       ilr.weights='blw.sqrt')

  philr_train_df <- philr_train %>% as_tibble(rownames = "label") %>%
    inner_join(halpern_impact %>%
                 select(label, median)) %>%
    #mutate(class = ohi_median2class_short(median)) %>%
    as.data.frame() %>%
    column_to_rownames("label")


  philr_test <- philr(otu_table_test, phy_tree(physeq_test),
                      part.weights='enorm.x.gm.counts',
                      ilr.weights='blw.sqrt')

  philr_test_df <- philr_test %>% as_tibble(rownames = "label") %>%
    inner_join(halpern_impact %>%
                 select(label, median)) %>%
    #mutate(class = ohi_median2class_short(median)) %>%
    as.data.frame() %>%
    column_to_rownames("label")

  list(train = philr_train_df, test = philr_test_df)
}



#split_h2o <- h2o.splitFrame(df, c(0.7, 0.15), seed = seed)
# set.seed(seed)
#df_train <- caret::createDataPartition(as.factor(as.data.frame(df)$class), p = 0.67, list = FALSE, times = 5)
#seed <- 12345
seed <- 222674
set.seed(seed)
k <- 5
h2o_automl_results <- vector(mode = "list")
df_folds <- caret::createFolds(halpern_impact %>% mutate(class = ohi_median2class_short(median)) %>% .$class %>% as.factor, k = k, returnTrain = TRUE)

#df <- df[,1:(ncol(df) - 1)]
for (o in 1:k) {
  cat(paste0("Iter: ", o, "\n"))
  h2o_automl_results[[o]] <- o

  df_list <- create_df(df_folds[[o]], osd2014_dada2_phyloseq_alpha_norare)

  train_conv_h2o <- as.h2o(df_list$train)
  test_conv_h2o  <- as.h2o(df_list$test)

  #train_conv_h2o <- as.h2o(df[ df_folds[[o]],])
  #test_conv_h2o  <- as.h2o(df[-df_folds[[o]],])

  #train_conv_h2o <- h2o.assign(split_h2o[[1]], paste0("train_", seed)) # 70%
  #valid_conv_h2o <- h2o.assign(split_h2o[[2]], paste0("valid_", seed)) # 15%
  #test_conv_h2o  <- h2o.assign(h2o.rbind(split_h2o[[2]], split_h2o[[3]]), paste0("test_", seed))  # 15%
  #test_conv_h2o  <- h2o.assign(split_h2o[[3]], paste0("test_", seed))  # 15%

  #

  lambda  <- forecast::BoxCox.lambda(train_conv_h2o %>% as.tibble() %>% .$median)

  train_medians <- tibble(median = train_conv_h2o %>% as_tibble() %>% .$median, class = "1_train") %>%
    mutate(median_bc = forecast::BoxCox(median, lambda), median_log = log(median), median_sqrt = sqrt(median))
  #valid_medians <- tibble(median = valid_conv_h2o %>% as_tibble() %>% .$median, class = "2_valid") %>%
  #   mutate(median_bc = forecast::BoxCox(median, lambda), median_log = log(median), median_sqrt = sqrt(median))
  test_medians <- tibble(median = test_conv_h2o %>% as_tibble() %>% .$median, class = "3_test") %>%
    mutate(median_bc = forecast::BoxCox(median, lambda), median_log = log(median), median_sqrt = sqrt(median))

  train_medians$median %>% skimr::skim()
  #valid_medians$median %>% skimr::skim()
  test_medians$median %>% skimr::skim()

  split_plot <- ggpubr::ggarrange(
    bind_rows(train_medians, test_medians) %>% ggplot(aes(median, color = class)) +
      geom_density() + ggtitle("RAW") +
      theme_bw(),
    bind_rows(train_medians, test_medians) %>% ggplot(aes(median_bc, color = class)) +
      geom_density() + ggtitle("BOX_COX") +
      theme_bw(),
    bind_rows(train_medians, test_medians) %>% ggplot(aes(median_log, color = class)) +
      geom_density() + ggtitle("LOG") +
      theme_bw(),
    bind_rows(train_medians, test_medians) %>% ggplot(aes(median_sqrt, color = class)) +
      geom_density() + ggtitle("SQRT") +
      theme_bw(),
    nrow = 1, ncol = 4, common.legend = TRUE
  )

  y <- "median"
  x <- setdiff(names(df_list$train), y)

  # For binary classification, response should be a factor
  if (y == "median"){
    train_conv_h2o[,y] <- h2o::as.numeric(train_conv_h2o[,y])
    #valid_conv_h2o[,y] <- h2o::as.numeric(valid_conv_h2o[,y])
    test_conv_h2o[,y] <- h2o::as.numeric(test_conv_h2o[,y])
  }else{
    train_conv_h2o[,y] <- h2o::as.factor(train_conv_h2o[,y])
    #valid_conv_h2o[,y] <- h2o::as.factor(valid_conv_h2o[,y])
    test_conv_h2o[,y] <- h2o::as.factor(test_conv_h2o[,y])
  }

  aml <- vector(mode = "list")
  #for (i in 1:10) {
  aml <- h2o.automl(y = y, x = x,
                    training_frame    = train_conv_h2o,
                    leaderboard_frame = test_conv_h2o,
                    #validation_frame = valid_conv_h2o,
                    #leaderboard_frame = test,
                    nfolds = 5,
                    #stopping_metric = "AUTO",
                    max_models = 100,
                    keep_cross_validation_predictions = TRUE,
                    #project_name = as.character(Sys.time()),
                    project_name = as.character(seed),
                    #project_name = "test_log2",
                    seed = seed, #sample(1:seed, 1)#,
                    exclude_algos = c("DeepLearning", "GLM")
  )
  #}


  h2o_automl_results[[o]]$aml <- aml
  h2o_automl_results[[o]]$seed <- seed
  #h2o_automl_results[[o]]$split <- split_h2o
  h2o_automl_results[[o]]$train <- as_tibble(train_conv_h2o)
  #h2o_automl_results[[o]]$valid <- as_tibble(valid_conv_h2o)
  h2o_automl_results[[o]]$test <- as_tibble(test_conv_h2o)
  h2o_automl_results[[o]]$split_plot <- split_plot
}

names(h2o_automl_results) <- paste0("rep_", 1:k)
automl_leader  <- lapply(h2o_automl_results, function(X) X$aml@leader)

# Try Ensemble models
#ensembles <- lapply(h2o_automl_results, get_ensemble)

# BEGIN: Save the models --------------------------------------------------
save_model_data(type = "18S", local_path = "osd2014_rf/data/h2o_models", remote_path = "/home/afernand/h2o_models")
# END: Save the models ----------------------------------------------------


# BEGIN: Reload the models ------------------------------------------------
# h2o.loadModel(path = "/home/afernand/h2o_models/18S/2018-11-28/rep_1/StackedEnsemble_model_R_1542698203489_190813")
#
# h2o.loadModel(path = "/home/afernand/h2o_models/18S/2018-11-28/rep_2/StackedEnsemble_model_R_1542698203489_190814")
#
# h2o.loadModel(path = "/home/afernand/h2o_models/18S/2018-11-28/rep_3/StackedEnsemble_model_R_1542698203489_190815")
#
# h2o.loadModel(path = "/home/afernand/h2o_models/18S/2018-11-28/rep_4/StackedEnsemble_model_R_1542698203489_190816")
#
# h2o.loadModel(path = "/home/afernand/h2o_models/18S/2018-11-28/rep_5/StackedEnsemble_model_R_1542698203489_190817")
#

h2o.loadModel(path = "/Users/ufo/Desktop/osd2014_analyses_repo/osd2014_analysis/osd2014_rf/data/h2o_models/18S/2018-12-14/rep_5/DeepLearning_grid_1_AutoML_20181213_153326_model_5")
h2o.loadModel(path = "/Users/ufo/Desktop/osd2014_analyses_repo/osd2014_analysis/osd2014_rf/data/h2o_models/18S/2018-12-14/rep_4/StackedEnsemble_AllModels_AutoML_20181213_145215")
h2o.loadModel(path = "/Users/ufo/Desktop/osd2014_analyses_repo/osd2014_analysis/osd2014_rf/data/h2o_models/18S/2018-12-14/rep_3/GBM_grid_1_AutoML_20181213_142127_model_2")
h2o.loadModel(path = "/Users/ufo/Desktop/osd2014_analyses_repo/osd2014_analysis/osd2014_rf/data/h2o_models/18S/2018-12-14/rep_2/DeepLearning_grid_1_AutoML_20181213_140041_model_3")
h2o.loadModel(path = "/Users/ufo/Desktop/osd2014_analyses_repo/osd2014_analysis/osd2014_rf/data/h2o_models/18S/2018-12-14/rep_1/DeepLearning_grid_1_AutoML_20181213_131919_model_3")

load(file = "osd2014_rf/data/h2o_automl_results_20181214.Rda")

# ensembles <- list(
#   h2o.getModel("StackedEnsemble_model_R_1542698203489_190813"),
#   h2o.getModel("StackedEnsemble_model_R_1542698203489_190814"),
#   h2o.getModel("StackedEnsemble_model_R_1542698203489_190815"),
#   h2o.getModel("StackedEnsemble_model_R_1542698203489_190816"),
#   h2o.getModel("StackedEnsemble_model_R_1542698203489_190817")
# )

# END: Reload the models --------------------------------------------------



ohi_thrs <- c("Very low impact", "Low impact", "Medium impact", "Medium high impact", "High impact", "Very high impact")
#
# for (i in paste0("rep_", 1:5)){
#   cat(i)
#   seed <- h2o_automl_results[[i]]$seed
#   split_h2o <- h2o.splitFrame(df, c(0.7, 0.15), seed = seed)
#   h2o_automl_results[[i]]$train <- h2o.assign(split_h2o[[1]], paste0("train_", seed)) # 70%
#   h2o_automl_results[[i]]$valid <- h2o.assign(split_h2o[[2]], paste0("valid_", seed)) # 15%
#   h2o_automl_results[[i]]$test  <- h2o.assign(split_h2o[[3]], paste0("test_", seed))  # 15%
# }

# Predictions and performance for individual/ensemble models --------------
automl_preds  <- lapply(h2o_automl_results, function(X) predict(object = X$aml@leader, newdata = as.h2o(X$test)))
automl_perf <- lapply(h2o_automl_results, function(X) h2o.performance(X$aml@leader, newdata = as.h2o(X$test)))

# ensemble_preds  <- lapply(ensembles, function(X) h2o.predict(object = X$ensemble, newdata = as.h2o(X$test)))
# ensemble_perf <- lapply(ensembles, function(X) h2o.performance( X$ensemble, newdata = as.h2o(X$test)))


# Estimate errors for the different models --------------------------------
automl_errors <- lapply(names(h2o_automl_results), get_error, data = h2o_automl_results)
names(automl_errors) <- paste0("rep_", 1:k)

# automl_errors_ensemble <- lapply(ensembles, get_error_ensemble)

automl_errors_summary <- map_df(automl_errors, summarise_errors, .id = "rep")
automl_errors_summary %>% arrange(MAE)
# automl_errors_summary_ensemble <- lapply(automl_errors_ensemble, summarise_errors)

automl_plots <- lapply(automl_errors, plot_error_point)
ggarrange(plotlist = automl_plots)
# automl_plots_ensemble <- lapply(automl_errors_ensemble, plot_error_point)
# ggarrange(plotlist = automl_plots_ensemble)


automl_plots_residuals <- lapply(automl_errors, plot_residuals_point)
ggarrange(plotlist = automl_plots_residuals)

# automl_plots_residuals_ensemble <- lapply(automl_errors_ensemble, plot_residuals_point)

automl_confusion <- lapply(automl_errors, function(X) caret::confusionMatrix(X$class_pred, X$class_actual))
# automl_confusion_ensemble <- lapply(automl_errors_ensemble, function(X) caret::confusionMatrix(X$class_pred, X$class_actual))

map_df(automl_confusion,broom::tidy, .id = "repetition") %>%
  dplyr::filter(term %in% c("accuracy", "precision", "recall", "sensitivity", "specificity")) %>%
  ggplot(aes(repetition, estimate)) +
  geom_col() +
  ggpubr::rotate() +
  facet_wrap(~term, nrow = 1)

# # Stats for predicted class -----------------------------------------------
#
# test_h2o_df <- as_tibble(test_conv_h2o)
#
# test_h2o_2 = test_h2o_df %>%
#   as.data.frame() %>%
#   mutate(sample_id = rownames(test_h2o_df ))
#
# test_correct_rnumbers <- test_performance %>%
#   mutate(sample_id = rownames(test_performance)) %>%
#   filter(correct == 'correct')
#
# test_correct <- test_h2o_2 %>%
#   as_tibble() %>%
#   filter(row_number() %in% test_correct_rnumbers$sample_id)
#
# test_wrong_rnumbers <- test_performance %>%
#   mutate(sample_id = rownames(test_performance)) %>%
#   filter(correct != 'correct')
#
# test_wrong <- test_h2o_2 %>%
#   as_tibble() %>%
#   filter(row_number() %in% test_wrong_rnumbers$sample_id)
#
# explainer <- lime::lime(
#   as.data.frame(train_conv_h2o) %>% dplyr::select(-class),
#   model          = automl_leader,
#   bin_continuous = FALSE)
#
# explanation_corr <- lime::explain(
#   x = as.data.frame(test_correct) %>% dplyr::select(-class, -sample_id),
#   explainer      = explainer,
#   n_labels = 1,
#   n_features = 10,
#   kernel_width = 0.5)
#
# plot_features(explanation_corr, ncol = 4)
#
# explanation_wrong <- lime::explain(
#   x = as.data.frame(test_wrong) %>% select(-class, -sample_id),
#   explainer      = explainer,
#   n_labels = 1,
#   n_features = 10,
#   kernel_width = 0.5)
#
# plot_features(explanation_wrong, ncol = 4)
#
# importance_h2o <- h2o.varimp(automl_leader)
#
#
#
# explainer <- lime::lime(
#   as.data.frame(train_conv_h2o) %>% select(-median),
#   model          = automl_leader,
#   bin_continuous = TRUE)
#
# explanation <- lime::explain(
#   x = as.data.frame(test_conv_h2o) %>% select(-median),
#   explainer      = explainer,
#   n_features = 10,
#   kernel_width = 0.5)
#
# plot_features(explanation_corr, ncol = 4)

# Get the most importance balances ----------------------------------------

# We will get the balances, and their position in the tree and the level

# Get importants and merge the

# Get importance from the individual models

importance_h2o <- h2o.varimp(h2o_automl_results$rep_7$aml@leader)

importance_h2o <- map_df(h2o_automl_results, function(X){
  f <-  h2o.varimp(X$aml@leader)
  f1 <- tibble(pos = brStick(f$scaled_importance)[["Use PCs:"]]) %>%
    mutate(pos_difference = pos - lag(pos), pos_difference = ifelse(is.na(pos_difference), 1, pos_difference)) %>%
    filter(pos_difference > 1)
  if (nrow(f1) == 0) {
    sel <- f[brStick(f$scaled_importance)[["Use PCs:"]], ]
  }else{
    sel <- head(f, f1[1,]$pos - 1)
  }

  h2o.varimp(X$aml@leader) %>% mutate(rank = row_number()) %>% as_tibble() %>% filter(variable %in% sel$variable)
}, .id = "fold")

importance_h2o_filt <- importance_h2o %>% select(variable, rank) %>% group_by(variable) %>% summarise(N = n(), median_rank = median(rank)) %>% arrange(median_rank, desc(N))


importance_h2o %>%
  as_tibble() %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(rank, relative_importance, group = 1)) +
  geom_line()




importance_h2o_filt <- importance_h2o %>%
  as_tibble() %>%
  mutate(rank = row_number()) %>% filter(rank <= 29)

importance_h2o_filt <- importance_h2o %>%
  as_tibble() %>%
  filter(relative_importance > 0)
importance_h2o_filt

balance_names <- map_df(importance_h2o_filt$variable,
                        function(x) {
                          tibble(node = x,
                                 balance = name.balance(tree, tax[,-ncol(tax)], x)) %>%
                            separate(col = "balance", into = c("numerator", "denominator"), sep = "/", remove = FALSE)})

explore_balances <- function(X, tree = tree){

  sub_tree <- tree %>% tree_subset_mod(node = X, levels_back = 0)

  if (!is.binary(sub_tree)){
    sub_tree <- multi2di(sub_tree)
  }

  balance_node_id <- name.to.nn(sub_tree, X)

  balance_name <- name.balance(sub_tree, tax[,-1], X)
  balance_node_order <- philr:::get.ud.nodes(sub_tree, X)

  up_balance <- plyr::mlply(.data=name.to.nn(sub_tree, balance_node_order$up), .fun=getChildren, tree=sub_tree, .parallel=FALSE)[[1]] %>%
    as_tibble() %>%
    mutate(balance = paste0(X, "+")) %>%
    dplyr::rename(label = value)

  down_balance <- plyr::mlply(.data=name.to.nn(sub_tree, balance_node_order$down), .fun=getChildren, tree=sub_tree, .parallel=FALSE)[[1]] %>%
    as_tibble() %>%
    mutate(balance = paste0(X, "-")) %>%
    dplyr::rename(label = value)

  sub_tree_df <- bind_rows(down_balance, up_balance) %>%
    right_join(sub_tree %>% tidytree::as_tibble()) %>%
    left_join(as(tax, "matrix") %>% as_tibble(rownames = "label")) %>%
    left_join(asv_names %>% dplyr::rename(label = asv))

  sub_tree_tips <- sub_tree_df %>%
    filter(!is.na(asv_name))

  sub_tree$tip.label <- plyr::mapvalues(sub_tree$tip.label, from = asv_names$asv, to = asv_names$asv_name, warn_missing = FALSE)

  r_color <- randomcoloR::randomColor(count = 1)

  p <- ggtree(sub_tree) +
    geom_balance(node = balance_node_id, fill = r_color, alpha=0.6) +
    geom_nodelab() +
    geom_tiplab() +
    ggtitle(balance_name)

  p <- annotate_balance(sub_tree, X, p=p, labels = c(paste0(X, "+"), paste0(X, "-")), bar=TRUE)

  list(balance_node_id = balance_node_id, sub_tree = sub_tree, sub_tree_plot = p,
       balance_name = balance_name, balance_node_order = balance_node_order,
       sub_tree_df = sub_tree_df, sub_tree_tips = sub_tree_tips)
}


library(ggtree)
library(tidytree)
library(treeio)
balances_strees <- pbmcapply::pbmclapply(balance_names$node, explore_balances, tree = tree)
names(balances_strees) <- balance_names$node
n <- 10
gp.philr.long <- osd2014_dada2_phyloseq_beta_df %>%
  rownames_to_column("label") %>%
  as_tibble() %>%
  gather(coord, value, c(-class, -median, -label)) %>%
  mutate(class_1 = ifelse(grepl("low", class, ignore.case = TRUE), "LOW", "MEDIUM"),
         #coord = fct_relevel(coord, importance_h2o$variable),
         class_2 = ohi_median2class(median))

ggplot(gp.philr.long %>% filter(coord %in% importance_h2o_filt$variable[1:n]), aes(x=class_1, y=value)) +
  geom_boxplot(fill='lightgrey') +
  # ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 1,
  #                       color = "black", alpha = 1, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4, errorbar.length = 0.2) +
  geom_point() +
  facet_wrap(.~coord, scales='free') +
  xlab('Impact') + ylab('Balance Value') +
  theme_bw()

# Calculate variance

BYSITE <- list()
for (site in c("LOW", "MEDIUM")){
  #for (site in c("Medium impact","Low impact", "Very low impact")){
  sel_samples <- gp.philr.long %>% select(label, class_1) %>% unique() %>% filter(class_1 == site) %>% .$label
  SITE <- subset_samples(osd2014_dada2_phyloseq_alpha_norare, sample_names(osd2014_dada2_phyloseq_alpha_norare) %in% sel_samples)
  SITE <- filter_taxa(SITE, function(x) sum(x > 1) > (0.2*length(x)) , TRUE)
  SITE <- prune_samples(rowSums(otu_table(SITE)) > 50, SITE)
  BYSITE[[site]] <- SITE
}

bysite = list()
for (site in c("LOW", "MEDIUM")){
  #for (site in c("Medium impact","Low impact", "Very low impact")){
  bysite[[site]] <- run.var.analysis(BYSITE[[site]], n=20000, n.support = 10, ncores=3)
}


# Get p-values
one.sided.pval <- bysite %>%
  map(~ecdf(.x[['nmr']][['log_mdtt']])(.x[['lm']][['coefficients']][2,1])) %>%
  as_vector()
two.sided.pval <- ifelse(one.sided.pval < 0.5, one.sided.pval*2, (1-one.sided.pval)*2)

#sink(file='results/permutation_pval.txt')
cat('BH corrected two-sided tests from Permutation:\n')
cat("\n")
(two.sided.pval <- p.adjust(two.sided.pval, method='BH'))



# Plot trees for each type ------------------------------------------------

plot_tree_var <- function(site, to.plot){
  # Simply Phyla Levels
  tax_var <- as(tax_table(BYSITE[[site]]), "matrix") %>% as.data.frame()
  tax_var <- tax_var[,'Supergroup', drop=F]
  simple.tax <- simplify.phylum(tax_var$Phylum)
  tax$Phylum <- simple.tax$tax

  # create tree and color
  tr <- phy_tree(BYSITE[[site]])
  p <- ggtree(tr, aes(color=var.tz, size=isTip)) %<+% as.data.frame(bysite[[site]]$var.tz) +
    scale_size_manual(values=c(0.5, .25)) +
    labs(color='Balance Variance')

  # Add labels to highlihgted points
  if (!is.null(to.plot)){
    data <- p$data %>% filter(label %in% to.plot)
    #p <- p + geom_label_repel(aes(label=label), size=4, data=data)
    p <- p + geom_point2(aes(label = label), size = 3, data=data)
  }

  # Don't Add Phyla bar
  p <- p +
    theme(legend.position='bottom') +
    scale_color_gradient2(low='red', mid='darkgrey', high='blue', trans='log10',
                          midpoint = 0.5) +
    guides(size='none') +
    geom_treescale()
  p_medium <- p
}

tree_var_low <- plot_tree_var('LOW', to.plot = NULL)
tree_var_medium <- plot_tree_var('MEDIUM', to.plot = NULL)


# Plot variance by distance to tips ---------------------------------------
pvd_low <- figure.pvd_18S('LOW', NULL, plot.rug = F, plot.rank.lines = T, point.size = 1.5) +
  scale_x_log10() +
  scale_y_log10(limits = c(0.08, 15)) +
  xlab("Mean Distance to Tips") +
  ylab("Balance Variance") +
  ggtitle("LOW")

pvd_high <- figure.pvd_18S('MEDIUM', NULL, plot.rug = F, plot.rank.lines = T, point.size = 1.5) +
  scale_x_log10() +
  scale_y_log10(limits = c(0.08, 15)) +
  xlab("Mean Distance to Tips") +
  ylab("Balance Variance") +
  ggtitle("MEDIUM")

# Explore variances for the most important --------------------------------

bind_rows(bysite[["LOW"]]$var.tz %>% select(coord, var.tz, mean.dist.to.tips) %>% mutate(class = "LOW"),
          bysite[["MEDIUM"]]$var.tz %>% select(coord, var.tz, mean.dist.to.tips) %>% mutate(class = "MEDIUM")) %>%
  filter(coord %in% balance_names$node) %>%
  ggplot(aes(mean.dist.to.tips, var.tz, fill = class)) +
  # ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 1,
  #                       color = "black", alpha = 1, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4, errorbar.length = 0.2) +
  geom_point(shape = 21, color = "black") +
  theme_bw()


ggarrange(tree_var_low, pvd_low,
          tree_var_medium, pvd_high,
          ncol = 2, nrow = 2, common.legend = TRUE)



# philR functions for plotting --------------------------------------------



ggpubr::ggarrange(pvd_low, pvd_high, ncol = 2, nrow = 1, hjust = T, vjust = T)

figure.gms <- function(site, low, high){
  gm <- calc.gm.coords(site, c(high, low))
  gm$coord <- factor(gm$coord, levels=c(high, low))
  gms <- gm %>%
    na.omit() %>%
    ggplot(aes(x=gm.down, y=gm.up)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(coord~., scales = 'free') +
    theme_bw() +
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(size=16),
          axis.text.y=element_text(size=16))
}


cs <- balance_names$node
calc.gm.coords("LOW", balance_names$node) %>% as_tibble()
calc.gm.coords("MEDIUM", balance_names$node) %>% as_tibble()

site <- "MEDIUM"
calc.gm.coords(site, balance_names$node) %>% t() %>% as_tibble() %>% select(coord, var) %>% unique() %>% arrange(desc(var))
calc.gm.coords(site, balance_names$node) %>%
  as_tibble() %>%
  filter(coord == "n3") %>%
  ggplot(aes(gm_up, gm_down)) +
  geom_point() + scale_x_log10() + scale_y_log10()

# Plot variances in the tree for the two classes --------------------------
bottom.10_low <- bysite[["LOW"]]$var.tz %>%
  filter(mean.dist.to.tips > 0) %>%
  top_n(-10, var.tz) %>%
  arrange(var.tz) %>%
  mutate(class = "LOW")
bottom.10_medium <- bysite[["MEDIUM"]]$var.tz %>%
  filter(mean.dist.to.tips > 0) %>%
  top_n(-10, var.tz) %>%
  arrange(var.tz) %>%
  mutate(class = "MEDIUM")

top.10_low <- bysite[["LOW"]]$var.tz %>%
  filter(mean.dist.to.tips > 0) %>%
  top_n(10, var.tz) %>%
  arrange(var.tz) %>%
  mutate(class = "LOW")
top.10_medium <- bysite[["MEDIUM"]]$var.tz %>%
  filter(mean.dist.to.tips > 0) %>%
  top_n(10, var.tz) %>%
  arrange(var.tz) %>%
  mutate(class = "MEDIUM")


# Get network data --------------------------------------------------------

load("osd2014_18S_asv/data/osd2014_18S_asv_networks_results.Rdata", verbose = TRUE)

test_balances <- balances_strees[[1]]$sub_tree_tips

test_balances %>% inner_join(df_nodes %>% dplyr::rename(label = asv)) %>% group_by(com) %>% dplyr::count()

g_comp <- g %>% as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree(), weighted_degree = centrality_degree() / local_ave_degree()) %>%
  inner_join(as(tax_table(osd2014_dada2_phyloseq_alpha), "matrix") %>% as_tibble(rownames = "asv")) %>%
  arrange(desc(weighted_degree))

get_g_balances <- function(X){

  vids <- g_comp %>% activate(nodes) %>% filter(asv %in% X$label) %>% as_tibble() %>% .$name
  list_of_edges <- E(g_comp)[from(vids) | to(vids)]
  g_subgraph <- subgraph.edges(g_comp, list_of_edges) %>% as_tbl_graph() %>% activate(nodes)

  g_subgraph %>% as_tibble() %>% filter(name %in% vids)
  g_subgraph %>% as_tibble() %>% ggplot(aes(weighted_degree)) + geom_histogram()

}


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

