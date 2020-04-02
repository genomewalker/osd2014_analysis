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
source("osd2014_16S_asv/lib/fuzzyforest_lib.R")

load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_with_phylo.Rdata")
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

halpern_impact$class %>% table()
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


osd2014_dada2_phyloseq_alpha_filt <- subset_samples(osd2014_dada2_phyloseq_alpha, label %in% halpern_impact$label)
osd2014_dada2_phyloseq_alpha_filt <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_alpha_filt) > 0, osd2014_dada2_phyloseq_alpha_filt)
rare_taxa <- microbiome::rare_members(osd2014_dada2_phyloseq_alpha_filt, prevalence = 20/100)
osd2014_dada2_phyloseq_alpha_norare <- microbiome::remove_taxa(rare_taxa, osd2014_dada2_phyloseq_alpha_filt)
#osd2014_dada2_phyloseq_alpha_norare <- proportion(osd2014_dada2_phyloseq_alpha_norare)

phi <- propr((as(otu_table(osd2014_dada2_phyloseq_alpha_norare), "matrix")), metric = "rho")
hc <- caret::findCorrelation(phi@matrix, cutoff = 0.9)
hc <- sort(hc)
reduced_Data <- colnames(phi@matrix)[c(hc)]


osd2014_dada2_phyloseq_alpha_norare <- microbiome::remove_taxa(reduced_Data, osd2014_dada2_phyloseq_alpha_norare)

#osd2014_dada2_phyloseq_alpha_norare <- microbiome::transform(osd2014_dada2_phyloseq_alpha_norare, transform = "clr")
osd2014_dada2_phyloseq_alpha_norare <- transform_sample_counts(osd2014_dada2_phyloseq_alpha_norare, function(x) x+0.0001)


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
               #select(label, median)) %>%
               select(label, class)) %>%
  as.data.frame() %>%
  column_to_rownames("label")

#osd2014_dada2_phyloseq_beta_df_tax <- (as(tax_table(osd2014_dada2_phyloseq_alpha_norare), "matrix")) %>% as_tibble(rownames = "asv")

#colnames(osd2014_dada2_phyloseq_beta_df) <- plyr::mapvalues(colnames(osd2014_dada2_phyloseq_beta_df), from = osd2014_dada2_phyloseq_beta_df_tax$asv, to = osd2014_dada2_phyloseq_beta_df_tax$asv_name)



# osd2014_dada2_phyloseq_beta_df <- t(otu_table_rare_removed_norm) %>% as_tibble(rownames = "label") %>%
#   #select(top_features$feature_name, label) %>%
#   inner_join(halpern_impact %>%
#                select(label, median)) %>%
#   as.data.frame() %>%
#   column_to_rownames("label")
#
# dim(osd2014_dada2_phyloseq_beta_df)

h2o.init(nthreads = 64, ip = "localhost", port = 55577)


df <- as.h2o(osd2014_dada2_phyloseq_beta_df)

split_h2o <- h2o.splitFrame(df, c(0.7, 0.15), seed = 1242)
train_conv_h2o <- h2o.assign(split_h2o[[1]], "train" ) # 60%
valid_conv_h2o <- h2o.assign(split_h2o[[2]], "valid" ) # 20%
test_conv_h2o  <- h2o.assign(split_h2o[[3]], "test" )  # 20%


y <- "class"
x <- setdiff(names(df), y)

# For binary classification, response should be a factor
if (y == "median"){
  train_conv_h2o[,y] <- h2o::as.numeric(train_conv_h2o[,y])
  valid_conv_h2o[,y] <- h2o::as.numeric(valid_conv_h2o[,y])
  test_conv_h2o[,y] <- h2o::as.numeric(test_conv_h2o[,y])
}else{
  train_conv_h2o[,y] <- h2o::as.factor(train_conv_h2o[,y])
  valid_conv_h2o[,y] <- h2o::as.factor(valid_conv_h2o[,y])
  test_conv_h2o[,y] <- h2o::as.factor(test_conv_h2o[,y])
}
# Number of CV folds (to generate level-one data for stacking)
nfolds <- 10

aml <- h2o.automl(y = y, x = x,
                  training_frame    = train_conv_h2o,
                  leaderboard_frame = test_conv_h2o,
                  validation_frame = valid_conv_h2o,
                  #leaderboard_frame = test,
                  max_runtime_secs = 600,
                  #stopping_metric = "AUTO",
                  #max_models = 20,
                  seed = 1221222,
                  project_name = as.character(Sys.time())
)



automl_leader <- aml@leader

f1_th <- h2o.find_threshold_by_max_metric(h2o.performance(automl_leader, newdata = test_conv_h2o), "f1")


pred_conversion <- predict(object = automl_leader, newdata = test_conv_h2o)

pred_conversion <- pred_conversion %>% as.tibble() %>% mutate(predict_f1 = ifelse(low_impacted > f1_th, "low_impacted", "impacted"))

perf <- h2o.performance(automl_leader, newdata = test_conv_h2o)

test_performance <- test_conv_h2o %>%
  tibble::as_tibble() %>%
  select(class) %>%
  tibble::add_column(prediction = as.vector(pred_conversion$predict_f1)) %>%
  mutate(correct = ifelse(class == prediction, "correct", "wrong")) %>%
  mutate_if(is.character, as.factor)

confusion_matrix <- test_performance %>% select(-correct) %>% table()
confusion_matrix

tn <- confusion_matrix[1]
tp <- confusion_matrix[4]
fp <- confusion_matrix[3]
fn <- confusion_matrix[2]

accuracy <- (tp + tn) / (tp + tn + fp + fn)
misclassification_rate <- 1 - accuracy
recall <- tp / (tp + fn)
precision <- tp / (tp + fp)
null_error_rate <- tn / (tp + tn + fp + fn)

library(purrr)

tibble(
  accuracy,
  misclassification_rate,
  recall,
  precision,
  null_error_rate
) %>%
  purrr::transpose()


test_h2o_df <- as_tibble(test_conv_h2o)

test_h2o_2 = test_h2o_df %>%
  as.data.frame() %>%
  mutate(sample_id = rownames(test_h2o_df ))

test_correct_rnumbers <- test_performance %>%
  mutate(sample_id = rownames(test_performance)) %>%
  filter(correct == 'correct')

test_correct <- test_h2o_2 %>%
  as_tibble() %>%
  filter(row_number() %in% test_correct_rnumbers$sample_id)

test_wrong_rnumbers <- test_performance %>%
  mutate(sample_id = rownames(test_performance)) %>%
  filter(correct != 'correct')

test_wrong <- test_h2o_2 %>%
  as_tibble() %>%
  filter(row_number() %in% test_wrong_rnumbers$sample_id)

explainer <- lime::lime(
  as.data.frame(train_conv_h2o) %>% select(-class),
  model          = automl_leader,
  bin_continuous = FALSE)

explanation_corr <- lime::explain(
  x = as.data.frame(test_correct) %>% select(-class, -sample_id),
  explainer      = explainer,
  n_labels = 1,
  n_features = 10,
  kernel_width = 0.5)

plot_features(explanation_corr, ncol = 4)

explanation_wrong <- lime::explain(
  x = as.data.frame(test_wrong) %>% select(-class, -sample_id),
  explainer      = explainer,
  n_labels = 1,
  n_features = 10,
  kernel_width = 0.5)

plot_features(explanation_wrong, ncol = 4)

importance_h2o <- h2o.varimp(automl_leader)

tc.names <- sapply(importance_h2o$names[1:2], function(x) name.balance(tree, tax[,-1], x))
tc.names

library(ggtree)
plot_tree_and_labels <- function(tree){
  ggtree(tree)+
    geom_tiplab()+
    geom_text2(aes(subset=!isTip, label=label))
}


tc.nn <- name.to.nn(tree, importance_h2o$names[1:2])
tc.colors <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99')
p <- ggtree(tree, layout='fan') +
  geom_balance(node=tc.nn[1], fill=tc.colors[1], alpha=0.6) +
  geom_balance(node=tc.nn[2], fill=tc.colors[2], alpha=0.6)
#geom_balance(node=tc.nn[3], fill=tc.colors[3], alpha=0.6) +
#geom_balance(node=tc.nn[4], fill=tc.colors[4], alpha=0.6) +
#geom_balance(node=tc.nn[5], fill=tc.colors[5], alpha=0.6)
p <- annotate_balance(tree, 'n288', p=p, labels = c('n288+', 'n288-'),
                      offset.text=0.15, bar=FALSE)
annotate_balance(tree, 'n342', p=p, labels = c('n342+', 'n342-'),
                 offset.text=0.15, bar=FALSE)


gp.philr.long <- convert_to_long(osd2014_dada2_phyloseq_beta_df[,-1], osd2014_dada2_phyloseq_beta_df$class) %>%
  filter(coord %in% importance_h2o$names[1:10]) %>% as_tibble() %>% mutate(value = as.numeric(value))

ggplot(gp.philr.long, aes(x=labels, y=value)) +
  geom_boxplot(fill='lightgrey') +
  facet_wrap(.~coord, scales='free') +
  xlab('Human') + ylab('Balance Value') +
  theme_bw()

