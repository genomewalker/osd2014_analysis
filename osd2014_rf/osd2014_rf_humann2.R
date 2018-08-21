library(DESeq2)
library(tidyverse)
library(phyloseq)
library(caret)
library(ggpubr)
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)

osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf) %>%
  mutate(name = label) %>%
  as.data.frame() %>%
  column_to_rownames("name")


humann2_mod2label <- tbl(my_db, "humann2_mod2label") %>%
  collect(n = Inf)

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
  filter(buffer == "1km", ohi_variable == "global_cumul_impact") %>%
  mutate(class = ifelse(median <= 0, "non_impacted",
                        ifelse((median > 0.01 & median <= 4.8), "low_impacted",
                               "impacted"))) %>%
  filter(median != 0, class != "non_impacted")  %>%
  filter(class == "low_impacted" | (class == "impacted" & median > 6)) %>%
  left_join(osd2014_cdata) %>%
  #filter(dist_coast_m <= 2000) %>%
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
# sites_rand_low <- halpern_impact %>%
#   filter(class == "low_impacted") %>%
#   arrange(median) %>%
#   head(n = halpern_impact_class_min) %>%
#   .$label

# sites_rand_low <- prune_taxa("M00223", osd2014_h2mod_phyloseq) %>% otu_table() %>% as.matrix() %>% t %>% as.data.frame() %>% as.tibble(rownames = "label") %>% inner_join(halpern_impact) %>% filter(class =="low_impacted") %>% arrange(M00223) %>% head(19) %>% .$label

sites_rand_imp <- halpern_impact %>%
  filter(class == "impacted") %>%
  arrange(desc(median)) %>%
  sample_n(halpern_impact_class_min) %>%
  .$label


sites <- unique(c(sites_rand_low, sites_rand_imp))


# Read humann2 results
humman2_mod_abun <- tbl(my_db, "osd2014_humman2_mod_abun") %>%
  collect(n = Inf) %>%
  filter(module != "UNINTEGRATED" & module != "UNMAPPED") %>%
  spread(label, abun, fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("module")

humman2_path_abun <- tbl(my_db, "osd2014_humman2_path_abun") %>%
  collect(n = Inf) %>%
  filter(pathway != "UNINTEGRATED" & pathway != "UNMAPPED") %>%
  spread(label, abundance, fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("pathway")

humman2_ko_abun <- tbl(my_db, "osd2014_read_ko20140317_abun") %>%
  collect(n = Inf) %>%
  spread(label, abun, fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("ko")


humman2_mod_desc <- tbl(my_db, "kegg_module_description") %>%
  collect(n = Inf) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("module")

humman2_path_desc <- tbl(my_db, "kegg_pathway_description") %>%
  collect(n = Inf) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("pathway")

humman2_ko_desc <- tbl(my_db, "kegg_ko_description") %>%
  collect(n = Inf) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("ko_id")


load("~/osd2014_h2_rf.Rda")

osd2014_h2mod_phyloseq <- phyloseq(otu_table(as.matrix(humman2_mod_abun), taxa_are_rows = TRUE),
                                   tax_table(as.matrix(humman2_mod_desc)),
                                   sample_data(osd2014_cdata))

# Select samples training
osd2014_h2mod_phyloseq_train <- prune_samples(halpern_impact$label, osd2014_h2mod_phyloseq)

osd2014_h2mod_phyloseq_train <- prune_taxa(taxa_sums(osd2014_h2mod_phyloseq_train) > 0,
                                           osd2014_h2mod_phyloseq_train)



osd2014_h2mod_phyloseq_train_prop <- transform_sample_counts(osd2014_h2mod_phyloseq_train, function(x) x/sum(x))
osd2014_h2mod_phyloseq_train_h <- microbiome::transform(osd2014_h2mod_phyloseq_train, "hellinger")
osd2014_h2mod_phyloseq_train_c <- microbiome::transform(osd2014_h2mod_phyloseq_train, "clr")

osd2014_h2mod_phyloseq_train_ds <- osd2014_h2mod_phyloseq_train
osd2014_h2mod_phyloseq_train_deseq <- phyloseq_to_deseq2(osd2014_h2mod_phyloseq_train_ds, ~1)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(osd2014_h2mod_phyloseq_train_deseq), 1, gm_mean)
diagdds = estimateSizeFactors(osd2014_h2mod_phyloseq_train_deseq)
diagdds = estimateDispersions(diagdds)
diagvst = varianceStabilizingTransformation(diagdds, blind = FALSE)
diagvst <- assay(diagvst)
diagvst[diagvst < 0] <- 0
otu_table(osd2014_h2mod_phyloseq_train_ds) <- otu_table(diagvst, taxa_are_rows = TRUE)

#btex_kegg_modules$module
tax.mean <- taxa_sums(osd2014_h2mod_phyloseq_train_prop)/nsamples(osd2014_h2mod_phyloseq_train_prop)

osd2014_h2mod_phyloseq_train_filt <- prune_taxa(tax.mean > 2e-5, osd2014_h2mod_phyloseq_train_h)
h2mod_otutable <- phyloseq:::veganifyOTU(osd2014_h2mod_phyloseq_train_filt)


#otutable_filt <- h2mod_otutable[,c(as.character(h2mod_shared_impact_low_occ$h2mod), h2mod_impact_low$h2mod, h2mod_impact_imp$h2mod)]
#h2mod_otutable <- decostand(h2mod_otutable, "hellinger")
h2mod_otutable <- as.data.frame(h2mod_otutable)
#names(h2mod_otutable_pa) <- h2mod_name_conv$new
h2mod_otutable$label <- rownames(h2mod_otutable)
h2mod_otutable <- h2mod_otutable  %>% left_join(halpern_impact %>% dplyr::select(label, class))
class_impact <- h2mod_otutable %>% dplyr::select(label, class)
base::row.names(h2mod_otutable) <- h2mod_otutable$label
h2mod_otutable$label <- NULL
library(randomForest)



library(doParallel)
library(caret)

# Random Search
seed <- 1232
set.seed(seed)

#length is = (n_repeats*nresampling)+1
seeds <- vector(mode = "list", length = 50)

#(3 is the number of tuning parameter, mtry for rf, here equal to ncol(iris)-2)
for(i in 1:50) seeds[[i]] <- sample.int(n=1000, 15)

#for the last model
seeds[[51]] <- sample.int(1000, 1)
nt <- 1000
h2mod_otutable$class <- as.factor(h2mod_otutable$class)

model_weights <- ifelse(h2mod_otutable$class == "impacted",
                        (1/table(h2mod_otutable$class)[1]) * 0.5,
                        (1/table(h2mod_otutable$class)[2]) * 0.5)


set.seed(911)

training.samples <- h2mod_otutable$class %>%
  createDataPartition(p = 0.9, list = FALSE)

train.data <- h2mod_otutable[training.samples, ]
test.data <- h2mod_otutable[-training.samples, ]
train.data$class <- as.factor(train.data$class)
test.data$class <- as.factor(test.data$class)


cl <- makeCluster(4)
registerDoParallel(cl)
set.seed(seed)
control <- trainControl( method="repeatedcv", number=10, repeats=5, summaryFunction=prSummary, classProbs = TRUE ,search = "random", sampling = "smote")
train.data$class <- factor(as.character(train.data$class))
rf_smote <- train(y = train.data$class, x = train.data[1:dim(train.data)[2]-1], method="rf",
                  trControl=control, verboseIter = TRUE, savePredictions = TRUE,
                  prox=TRUE, allowParallel=TRUE,  metric = "AUC", tuneLength = 100, ntree = 1000,
                  importance = TRUE)
# sampsize = c('impacted'=as.integer(table(train.data$class)[[1]]*1),
#              'low_impacted'=as.integer(table(train.data$class)[[1]]*1)), ,
# strata=as.factor(train.data$class), )
#, ntree = 1000)weights = model_weights,

control <- trainControl( method="repeatedcv", number=10, repeats=5, summaryFunction=prSummary, classProbs = TRUE ,search = "random", sampling = "up")
train.data$class <- factor(as.character(train.data$class))
rf_up <- train(y = train.data$class, x = train.data[1:dim(train.data)[2]-1], method="rf",
               trControl=control, verboseIter = TRUE, savePredictions = TRUE,
               prox=TRUE, allowParallel=TRUE,  metric = "AUC", tuneLength = 100, ntree = 1000,
               importance = TRUE)

control <- trainControl( method="repeatedcv", number=10, repeats=5, summaryFunction=prSummary, classProbs = TRUE ,search = "random", sampling = "down")
train.data$class <- factor(as.character(train.data$class))
rf_down <- train(y = train.data$class, x = train.data[1:dim(train.data)[2]-1], method="rf",
                 trControl=control, verboseIter = TRUE, savePredictions = TRUE,
                 prox=TRUE, allowParallel=TRUE,  metric = "AUC", tuneLength = 100, ntree = 1000,
                 importance = TRUE)

model_weights <- ifelse(train.data$class == "impacted",
                        (1/table(train.data$class)[1]) * 0.5,
                        (1/table(train.data$class)[2]) * 0.5)
control <- trainControl( method="repeatedcv", number=10, repeats=5, summaryFunction=prSummary, classProbs = TRUE ,search = "random")
rf_w <- train(y = train.data$class, x = train.data[1:dim(train.data)[2]-1], method="rf",
              trControl=control, verboseIter = TRUE, savePredictions = TRUE,
              prox=TRUE, allowParallel=TRUE,  metric = "AUC", tuneLength = 100, ntree = 1000,
              weights = model_weights, importance = TRUE)


rf_s <- train(y = train.data$class, x = train.data[1:dim(train.data)[2]-1], method="rf",
              trControl=control, verboseIter = TRUE, savePredictions = TRUE,
              prox=TRUE, allowParallel=TRUE,  metric = "AUC", tuneLength = 100, ntree = 1000,

              sampsize = c('impacted'=as.integer(table(train.data$class)[[1]]*1), 'low_impacted'=as.integer(table(train.data$class)[[1]]*1)),
              strata=as.factor(train.data$class), importance = TRUE)


rf <- train(y = train.data$class, x = train.data[1:dim(train.data)[2]-1], method="rf",
            trControl=control, verboseIter = TRUE, savePredictions = TRUE,
            prox=TRUE, allowParallel=TRUE,  metric = "AUC", tuneLength = 100, ntree = 1000,
            importance = TRUE)

stopCluster(cl)


model_list <- list(original = rf,
                   weighted = rf_w,
                   up = rf_up,
                   down = rf_down,
                   smote = rf_smote,
                   ssize = rf_s)


calc_auprc <- function(model, data){

  index_class2 <- data$class == "low_impacted"
  index_class1 <- data$class == "impacted"

  predictions <- predict(model, data, type = "prob")

  PRROC::pr.curve(predictions$low_impacted[index_class2],
                  predictions$low_impacted[index_class1],
                  curve = TRUE)

}

# Get results for all 5 models
model_list_pr <- model_list %>%
  map(calc_auprc, data = test.data)

model_list_pr %>%
  map(function(the_mod) the_mod$auc.integral)

results_list_pr <- list(NA)
num_mod <- 1

for(the_pr in model_list_pr){

  results_list_pr[[num_mod]] <-
    data_frame(recall = the_pr$curve[, 1],
               precision = the_pr$curve[, 2],
               model = names(model_list_pr)[num_mod])

  num_mod <- num_mod + 1

}

results_df_pr <- bind_rows(results_list_pr)

custom_col <- c("#000000", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#A62830")

ggplot(aes(x = recall, y = precision, group = model),
       data = results_df_pr) +
  geom_line(aes(color = model), size = 1) +
  scale_color_manual(values = custom_col) +
  geom_abline(intercept =
                sum(test.data$class == "low_impacted")/nrow(test.data),
              slope = 0, color = "gray", size = 1) +
  theme_bw()

# Plotting results --------------------------------------------------------

top_features <- rf_down$finalModel$importance %>%
  as.tibble(rownames = "module") %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  inner_join(humman2_mod_desc %>% as.tibble(rownames = "module")) %>%
  top_n(n = 10, wt = MeanDecreaseAccuracy) %>%
  mutate(module = fct_reorder(module, MeanDecreaseAccuracy))
ggplot(top_features, aes(x = module, y = MeanDecreaseAccuracy,)) +
  geom_point(shape = 22, size = 2, fill = "grey50") +
  ggpubr::rotate() +
  theme_bw() +
  scale_y_continuous() +
  scale_fill_manual(values = c("black", "red")) +
  ggpubr::theme_cleveland() +
  ylab("Average importance") +
  theme(legend.position = "top")
ggsave(plot = last_plot(), filename = "osd2014_rf/figures/osd2014_h2_rf_imp_ohi.pdf", width = 6, height = 4)

qpsmelt <- function(X) {
  if (taxa_are_rows(X)) {
    count_table <- as(otu_table(X), "matrix") %>%
      as_tibble(rownames = "OTU") %>%
      gather(label, Abundance, -OTU)
  }else{
    count_table <-as(otu_table(X), "matrix") %>%
      as_tibble(rownames = "label") %>%
      gather(OTU, Abundance, -label)
  }
  sample_table <- as(sample_data(X), "matrix") %>%
    as_tibble()
  taxa_table <- as(tax_table(X), "matrix") %>%
    as_tibble(rownames = "OTU")
  count_table %>%
    left_join(sample_table, by = "label") %>%
    left_join(taxa_table,  by = "OTU") %>%
    filter(Abundance > 0)
}

osd2014_rescaled_2013_median_long <- tbl(my_db, "osd2014_halpern_scaled_median") %>%
  collect(n = Inf) %>%
  filter(buffer == "1km") %>%
  select(-buffer) %>%
  spread(ohi_variable, median, fill =0) %>%
  select(-global_cumul_impact, -global_cumul_impact_diff_2008)

osd2014_rf_h2mod_results <- qpsmelt(osd2014_h2mod_phyloseq_train_prop) %>%
  select(label, OTU, Abundance) %>%
  dplyr::rename(module = OTU) %>%
  inner_join(top_features) %>%
  filter(label %in% sites) %>%
  unique() %>%
  inner_join(osd2014_rescaled_2013_median_long %>% select(label, plumes_fert, plumes_pest, inorganic)) %>%
  inner_join(halpern_impact %>% select(label, class, median))

base_breaks <- function(n = 3){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

osd2014_rf_h2mod_results %>%
  #filter(module %in% (top_features %>% arrange(desc(MeanDecreaseAccuracy)) %>% head(10)  %>% .$module)) %>%
  select(label, module, class, Abundance, MeanDecreaseAccuracy) %>% unique %>%
  filter(Abundance > 0) %>% #compare_means(formula = prop ~ class, group.by = "com", p.adjust.method = "fdr") %>% View
  mutate(module = fct_reorder(module, -MeanDecreaseAccuracy)) %>%
  ggplot(aes(class, (Abundance), fill = class)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.7, jitter.size = 1,
                        color = "black", alpha = 0.7, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4, errorbar.length = 0.2) +
  scale_y_log10(labels = scales::percent, breaks = base_breaks()) +
  facet_wrap(~module, scales = "free",nrow = 2) +
  stat_compare_means(comparisons = list(c("impacted", "low_impacted")), aes(label=..p.adj..)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#2E3239", "#1F629A"))

ggsave(plot = last_plot(), filename = "osd2014_rf/figures/osd2014_h2_rf_bplots_ohi.pdf", width = 8, height = 3)


save.image("osd2014_rf/data/osd2014_h2_ohi.Rda")

# SANDBOX -----------------------------------------------------------------


osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)
osd2014_rf_h2mod_results %>%
  #filter(module %in% (top_features %>% arrange(desc(MeanDecreaseAccuracy)) %>% head(10)  %>% .$module)) %>%
  select(label, module, class, Abundance, MeanDecreaseAccuracy) %>% unique %>%
  filter(Abundance > 0) %>% #compare_means(formula = prop ~ class, group.by = "com", p.adjust.method = "fdr") %>% View
  mutate(module = fct_reorder(module, -MeanDecreaseAccuracy)) %>%
  inner_join(osd2014_cdata) %>%
  filter(module == "M00568", meow_province %in% osd2014_meow_regions$meow_province) %>%
  ggplot(aes(meow_province, (Abundance), fill = meow_province)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.7, jitter.size = 1,
                        color = "black", alpha = 0.7, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4, errorbar.length = 0.2) +
  scale_y_log10(labels = scales::percent, breaks = base_breaks()) +
  facet_wrap(~module, scales = "free",nrow = 2) +
  stat_compare_means(comparisons = list(c("impacted", "low_impacted")), aes(label=..p.adj..)) +
  theme_bw() +
  theme(legend.position = "none") + NULL
scale_fill_manual(values = c("#2E3239", "#1F629A"))


osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)
osd2014_rf_h2mod_results %>%
  #filter(module %in% (top_features %>% arrange(desc(MeanDecreaseAccuracy)) %>% head(10)  %>% .$module)) %>%
  select(label, module, class, Abundance, MeanDecreaseAccuracy) %>% unique %>%
  filter(Abundance > 0) %>% #compare_means(formula = prop ~ class, group.by = "com", p.adjust.method = "fdr") %>% View
  mutate(module = fct_reorder(module, -MeanDecreaseAccuracy)) %>%
  inner_join(osd2014_cdata) %>%
  filter(module == "M00623", meow_province %in% osd2014_meow_regions$meow_province) %>%
  ggplot(aes(meow_province, (Abundance), fill = meow_province)) +
  ggpol::geom_boxjitter(jitter.shape = 21, jitter.color = "black", jitter.alpha = 0.7, jitter.size = 1,
                        color = "black", alpha = 0.7, errorbar.draw = TRUE, jitter.height = 0.05, jitter.width = 0.075, width = 0.4, errorbar.length = 0.2) +
  scale_y_log10(labels = scales::percent, breaks = base_breaks()) +
  facet_wrap(~module, scales = "free",nrow = 2) +
  stat_compare_means(comparisons = list(c("impacted", "low_impacted")), aes(label=..p.adj..)) +
  theme_bw() +
  theme(legend.position = "none") + NULL
scale_fill_manual(values = c("#2E3239", "#1F629A"))

# 1 M00568 Catechol ortho-cleavage, catechol => 3-oxoadipate ## K03381 K01856 K03464 (K01055,K14727)
# 2 M00008 Entner-Doudoroff pathway, glucose-6P => glyceraldehyde-3P + pyruvate              0.0102
# 3 M00545 Trans-cinnamate degradation, trans-cinnamate => acetyl-CoA                        0.00903
# 4 M00623 Phthalate degradation, phthalate => protocatechuate                               0.00784
# 5 M00328 Hemophore/metalloprotease transport system                                        0.00664
# 6 M00595 Thiosulfate oxidation by SOX complex, thiosulfate => sulfate                      0.00589
# 7 M00192 Putative thiamine transport system                                                0.00513
# 8 M00078 Heparan sulfate degradation                                                       0.00296
# 9 M00442 Putative hydroxymethylpyrimidine transport system                                 0.00242
# 10 M00737 Bacitracin resistance, VraDE transporter                                          0.00227

osd2014_assm_emapper_results <- tbl(my_db, "osd2014_assm_emapper_results") %>%
  collect(n = Inf)

osd2014_orf_abundance <- tbl(my_db, "osd2014_orf_abundance") %>%
  collect(n = Inf)

osd2014_assm_emapper_results %>%
  dplyr::select(query_name, KEGG_KOs) %>%
  separate_rows(KEGG_KOs, sep = ",") %>%
  rename(gene_id = query_name) %>%
  filter(KEGG_KOs %in% c("K03381", "K01856", "K03464", "K01055", "K14727")) %>%
  inner_join(osd2014_orf_abundance) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  select(label, gene_id, KEGG_KOs) %>%
  write_tsv("~/Downloads/osd2014_M00568_kos.txt", col_names = FALSE)

kegg_ko_descriptions <- tbl(my_db, "kegg_ko_description") %>%
  collect(n = Inf)
kegg_ko_descriptions %>%
filter(ko_id %in% c("K03381", "K01856", "K03464", "K01055", "K14727"))






# Phosphonate detailed analysis -------------------------------------------
# Entry	M00223
# Name	Phosphonate transport system
# Definition	K02044+K02042+K02041


osd2014_kegg_20140317_abundance <- tbl(my_db, "osd2014_kegg_assm_abundance") %>%
  collect(n = Inf)
osd2014_orf_abun_sample <- tbl(my_db, "osd2014_orf_abundance") %>%
  collect(n = Inf)


osd2014_assm_emapper_results <- tbl(my_db, "osd2014_assm_emapper_results") %>%
  collect(n = Inf)
osd2014_assm_emapper_results_phos <- osd2014_assm_emapper_results %>%
  dplyr::rename(gene_id = query_name) %>%
  left_join(osd2014_orf_abundance) %>% filter(label %in% sites, (KEGG_KOs == "K02044" | KEGG_KOs == "K02042" | KEGG_KOs == "K02041")) %>%
  left_join(halpern_impact %>% dplyr::select(label, class))

osd2014_phosphonate <- osd2014_kegg_20140317_abundance %>% filter(label %in% sites, (ko_id == "K02044" | ko_id == "K02042" | ko_id == "K02041")) %>%
  left_join(halpern_impact %>% dplyr::select(label, class))


osd2014_assm_emapper_results_phos %>% group_by(label, class, KEGG_KOs) %>%
  summarise(n = sum(abun)) %>%
  ungroup() %>%
  mutate( label = fct_reorder(label, n)) %>%
  ggplot(aes(label, n, fill = class)) +
  geom_col(width = 0.15, color = "black", fill = "black") +
  geom_point(shape = 21, size = 2) +
  ggpubr::rotate() +
  theme_light() +
  facet_wrap(~KEGG_KOs, scales = "free_x")


osd2014_assm_emapper_results_phos %>%
  select(label, gene_id, KEGG_KOs, class) %>%
  write_tsv(path = "~/Downloads/K01055_entropy.txt", col_names = FALSE)

osd2014_phosphonate %>%
  select(label, gene_id, ko_id) %>%
  inner_join(halpern_impact %>% dplyr::select(label, class)) %>%
  write_tsv(path = "~/Downloads/osd2014_phosphonate_gene_ids.tsv", col_names = FALSE)


entr <- read_tsv("~/Downloads/K01055_entropy.txt", col_names = FALSE) %>%
  dplyr::rename(contig = X1, gene_id = X2, entropy = X3) %>%
  inner_join(osd2014_orf_abundance) %>%
  inner_join(halpern_impact)

entr %>%
  group_by(label, class) %>%
  summarise(entropy = mean(entropy), n = n()) %>%
  ggplot(aes(class, entropy, size = n)) +
  ggpol::geom_boxjitter(aes(size = n)) +
  #facet_wrap(~ko_id) +
  scale_y_log10() +
  stat_compare_means()






  # Prepare Pi and entropy data ---------------------------------------------
osd2014_phosphonate_pi <- read_tsv("~/Downloads/osd2014_phosphonate_pi.tsv", col_names = FALSE)
names(osd2014_phosphonate_pi) <- c("label", "contig", "gene_id", "variants", "pi")

osd2014_phosphonate %>%
  inner_join(osd2014_phosphonate_pi) %>%
  #group_by(label, class, ko_id) %>%
  #summarise(mu_pi = mean(pi)) %>%
  filter(pi > 0) %>%
  ggplot(aes(class, pi)) +
  geom_boxplot() +
  facet_wrap(~ko_id) +
  scale_y_log10() +
  stat_compare_means(aes(label = ..p.adj..))


osd2014_phosphonate_entr <- read_tsv("~/Downloads/K02042/osd2014_phosphonate_entropy_K02042.tsv", col_names = FALSE)
names(osd2014_phosphonate_entr) <- c("contig", "gene_id", "entropy")


osd2014_phosphonate_entr %>%
  separate(contig, sep = "_spades_", into = "label", extra = "drop") %>%
  inner_join(halpern_impact) %>%
  #group_by(label, class, ko_id) %>%
  #summarise(mu_entr = mean(entropy)) %>%
  filter(entropy > 0) %>%
  ggplot(aes(class, entropy)) +
  ggpol::geom_boxjitter() +
  #facet_wrap(~ko_id) +
  scale_y_log10() +
  stat_compare_means(aes(label = ..p.adj..))

# VCF genomepop -----------------------------------------------------------
library(PopGenome)
setwd("")
snp <- readData("/Users/ufo/Downloads/K02041/vcf/", format="VCF", FAST = TRUE)

# Let's predict

# First get samples not used for training
pred_sites <- out_rescaled_2013_median_long %>% filter(!(label %in% sites)) %>% .$label %>% unique


# Trim the original matrix to contain the non-training and the modules used for
h2mod_otutable_pred <- as.data.frame(as(otu_table(osd2014_16s_otuXsample_physeq_filt_prev_prop), "matrix"))
h2mod_otutable_pred <- base::as.data.frame(h2mod_otutable_pred)
h2mod_otutable_pred <- h2mod_otutable_pred[pred_sites, c(as.character(h2mod_shared_impact_low_occ$h2mod), h2mod_impact_low$h2mod, h2mod_impact_imp$h2mod)]

pred_h2mod <- predict(halpern_classify, newdata = h2mod_otutable_pred)

pred_impacted_metadata_h2mod <- data.frame(class = pred_h2mod) %>% mutate(label = names(pred_h2mod)) %>% filter(class == "impacted") %>% left_join(osd2014_metadata)
pred_low_impacted_metadata_h2mod <- data.frame(class = pred_h2mod) %>% mutate(label = names(pred_h2mod)) %>% filter(class != "impacted") %>% left_join(osd2014_metadata)

library(ggmap)

### Set a range

### Get a map

tmp <- pred_impacted_metadata %>% dplyr::select(label, start_lat,start_lon, dist_coast_iso3_code) %>%
  mutate(filename = paste("~/Downloads/", label, "_", dist_coast_iso3_code, "_h2mod_impacted.pdf", sep = ""))


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
  mutate(filename = paste("~/Downloads/", label, "_", dist_coast_iso3_code, "_h2mod_low_impacted.pdf", sep = ""))


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


pred_impacted_metadata_kegg_mg <- pred_impacted_metadata_kegg %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label
pred_impacted_metadata_h2mod_mg <- pred_impacted_metadata_h2mod %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label

pred_impacted_metadata_kegg_mg_l <- length(pred_impacted_metadata_kegg_mg)
pred_impacted_metadata_h2mod_mg_l <- length(pred_impacted_metadata_h2mod_mg)

length(intersect(pred_impacted_metadata_kegg_mg,pred_impacted_metadata_h2mod_mg))/(pred_impacted_metadata_kegg_mg_l + pred_impacted_metadata_h2mod_mg_l)

pred_low_impacted_metadata_kegg_mg <- pred_low_impacted_metadata_kegg %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label
pred_low_impacted_metadata_h2mod_mg <- pred_low_impacted_metadata_h2mod %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label

pred_low_impacted_metadata_kegg_mg_l <- length(pred_low_impacted_metadata_kegg_mg)
pred_low_impacted_metadata_h2mod_mg_l <- length(pred_low_impacted_metadata_h2mod_mg)

length(intersect(pred_low_impacted_metadata_kegg_mg,pred_low_impacted_metadata_h2mod_mg))/(pred_low_impacted_metadata_kegg_mg_l + pred_low_impacted_metadata_h2mod_mg_l)


length(as.character(osd2014_amp_mg_intersect$label)[as.character(pred_impacted_metadata_kegg$label)])
length(pred_impacted_metadata_h2mod$label)
length(pred_low_impacted_metadata_kegg$label)
length(pred_low_impacted_metadata_h2mod$label)

intersect(pred_impacted_metadata_kegg$label, pred_impacted_metadata_h2mod$label)
intersect(pred_low_impacted_metadata_kegg$label, pred_low_impacted_metadata_h2mod$label)
length(intersect(pred_impacted_metadata_kegg$label, pred_impacted_metadata_h2mod$label))
length(intersect(pred_low_impacted_metadata_kegg$label, pred_low_impacted_metadata_h2mod$label))




