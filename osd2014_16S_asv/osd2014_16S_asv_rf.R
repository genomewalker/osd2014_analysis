library(DESeq2)
library(tidyverse)
library(phyloseq)


# BEGIN: WARNING!!!! -------------------------------------------------------------
# You can access to the data used in this analysis in several ways:
# 1. You have a copy of the PostgreSQL DB
# 2. You downloaded the .Rdata files from http://osd2014.metagenomics.eu/ and placed them
#    in the data folder
# 3. You can load the files remotely, it might take a while when the file is very large
# END: WARNING!!!! -------------------------------------------------------------


# BEGIN: WARNING!!: This will load all the data and results for the analysis --------
# Uncomment if you want to use it. Some of the analysis step might require long
# computational times and you might want to use a computer with many cores/CPUs

# load("osd2014_16S_asv/data/osd2014_16S_alpha_diversity.Rdata", verbose = TRUE)
# load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_alpha_diversity.Rdata"), verbose = TRUE)

# END: WARNING!! ---------------------------------------------------------------



# BEGIN: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------

# Load necessary data -----------------------------------------------------
# Use if you have the postgres DB in place

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata")

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
  filter(buffer == "1km", ohi_variable == "global_cumul_impact") %>%
  mutate(class = ifelse(median <= 0, "non_impacted",
                        ifelse((median > 0.01 & median <= 4.8), "low_impacted",
                               "impacted"))) %>%
  filter(median != 0, class != "non_impacted")  %>% filter(class == "low_impacted" | (class == "impacted" & median > 6)) %>%
  left_join(osd2014_cdata) %>%
  filter(dist_coast_m <= 2000) %>%
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


halpern_impact$class %>% table


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
  arrange(median) %>%
  head(n = halpern_impact_class_min) %>%
  .$label

sites_rand_imp <- halpern_impact %>%
  filter(class == "impacted") %>%
  arrange(desc(median)) %>%
  head(n = halpern_impact_class_min) %>%
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

# Select samples training
osd2014_dada2_phyloseq_train <- prune_samples(sites, osd2014_dada2_phyloseq_alpha)

osd2014_dada2_phyloseq_train <- prune_taxa(taxa_sums(osd2014_dada2_phyloseq_train) > 0,
                                           osd2014_dada2_phyloseq_train)

osd2014_dada2_phyloseq_train_prop <- transform_sample_counts(osd2014_dada2_phyloseq_train, function(x)x/sum(x))

#library(microbiome)
#osd2014_dada2_phyloseq_train <- microbiome::transform(osd2014_dada2_phyloseq_train, "hellinger")


osd2014_dada2_phyloseq_train_deseq <- phyloseq_to_deseq2(osd2014_dada2_phyloseq_train, ~1)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(osd2014_dada2_phyloseq_train_deseq), 1, gm_mean)
diagdds = estimateSizeFactors(osd2014_dada2_phyloseq_train_deseq, type = "poscounts")
diagdds = estimateDispersions(diagdds)
diagvst = varianceStabilizingTransformation(diagdds, blind = FALSE)
diagvst <- assay(diagvst)
diagvst[diagvst < 0] <- 0
otu_table(osd2014_dada2_phyloseq_train) <- otu_table(diagvst, taxa_are_rows = TRUE)


tax.mean <- taxa_sums(osd2014_dada2_phyloseq_train_prop)/nsamples(osd2014_dada2_phyloseq_train_prop)

osd2014_dada2_phyloseq_train <- prune_taxa(tax.mean > 1e-5, osd2014_dada2_phyloseq_train)

asv_otutable <- phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_train)
asv_otutable_pa <- vegan::decostand(data.frame(impacted = colSums(asv_otutable[sites_rand_imp,]), low = colSums(asv_otutable[sites_rand_low,])), method = "pa")
dim(asv_otutable_pa)

asv_otutable <- phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_train)

UpSetR::upset(asv_otutable_pa, order.by = c("degree"), decreasing = c(TRUE))

asv_shared_impact_low <- as_data_frame(asv_otutable_pa) %>%
  rownames_to_column(var = "asv") %>%
  filter(impacted == 1, low == 1) %>%
  tbl_df()

asv_impact_low <- as_data_frame(asv_otutable_pa) %>%
  rownames_to_column(var = "asv") %>%
  filter(impacted == 0, low == 1) %>%
  tbl_df()

asv_impact_imp <- as_data_frame(asv_otutable_pa) %>%
  rownames_to_column(var = "asv") %>%
  filter(impacted == 1, low == 0) %>%
  tbl_df()

asv_shared_impact_low_occ <- as_data_frame(asv_otutable) %>%
  rownames_to_column(var = "label") %>%
  reshape2::melt() %>%
  tbl_df %>%
  rename(asv = variable) %>%
  filter(asv %in% asv_shared_impact_low$asv) %>%
  filter(value > 0) %>%
  group_by(asv) %>%
  count() %>%
  mutate(occ = n/nsamples(osd2014_dada2_phyloseq_train)) %>%
  filter(occ >= 0.25)


#otutable_filt <- asv_otutable[,c(as.character(asv_shared_impact_low_occ$asv), asv_impact_low$asv, asv_impact_imp$asv)]
#asv_otutable <- decostand(asv_otutable, "hellinger")
asv_otutable_pa <- as.data.frame(asv_otutable)


#names(asv_otutable_pa) <- asv_name_conv$new
asv_otutable_pa$label <- rownames(asv_otutable_pa)
asv_otutable_pa <- asv_otutable_pa  %>% left_join(halpern_impact %>% dplyr::select(label, class))
class_impact <- asv_otutable_pa %>% dplyr::select(label, class)
base::row.names(asv_otutable_pa) <- asv_otutable_pa$label
asv_otutable_pa$label <- NULL
asv_otutable_filt_pa <- asv_otutable_pa[sites,c(as.character(asv_shared_impact_low_occ$asv), asv_impact_low$asv, asv_impact_imp$asv,"class")]
library(randomForest)

#samp <- sample(nrow(asv_otutable), 2/3 * nrow(asv_otutable))
asv_otutable_filt_pa$class <- as.factor(asv_otutable_filt_pa$class)
#train <- as(asv_otutable[samp, ], "data.frame")
#train$class <- as.factor(train$class)

#train <- train[, colSums(train != 0) > 0]

#test <-  as(asv_otutable[-samp,], "data.frame")
#test$class <- as.factor(test$class)

#train_class <- class_impact %>% filter(label %in% rownames(train))
#train_class <- train_class[match(rownames(train), train_class$label),]

#test_class <- class_impact %>% filter(label %in% rownames(test))
#test_class <- test_class[match(rownames(test), test_class$label),]

asv_table <- asv_otutable_filt_pa
dim(asv_table)
# asv_table <- as.data.frame(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_train))
# imp_class <- halpern_impact %>%
#   select(label, class) %>%
#   filter(label %in% sample_names(osd2014_dada2_phyloseq_train)) %>%
#   slice(match(label, rownames(asv_table)))
# asv_table$class <- as.factor(imp_class$class)

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
for(i in 1:50) seeds[[i]] <- sample.int(n=1000, 15)

#for the last model
seeds[[51]] <- sample.int(1000, 1)
nt <- 1000
control <- trainControl(method="repeatedcv", number=10, repeats=5, search="random", seeds = seeds)
asv_rf_random <- train(y = asv_table$class, x=asv_table[1:dim(asv_table)[2]-1], method="rf", metric="Accuracy", tuneLength=15,
                       trControl=control, verboseIter = TRUE, savePredictions = TRUE,
                       prox=TRUE, allowParallel=TRUE, num.trees=nt)
stopCluster(cl)
print(asv_rf_random)
plot(asv_rf_random)

#mtry = 329

n <- 0.30

repeat{
  asv_halpern_classify <- randomForest::randomForest(y = asv_table$class,
                                                     x = asv_table[1:dim(asv_table)[2]-1], ntree=nt, importance=TRUE, do.trace=nt, keep.inbag = TRUE, mtry = asv_rf_random$bestTune$mtry)
  err <- asv_halpern_classify$err.rate[,1][nt]
  diff <- abs(asv_halpern_classify$err.rate[,2][nt] - asv_halpern_classify$err.rate[,3][nt])
  if(err < n && diff < 0.06){
    break
  }
}

#save(names_otus, rf_random, names_otus_short, halpern_classify, halpern_impact, asv_otutable, osd2014_16s_otuXsample_physeq_filt_prev_beta_mg, file = "~/ownCloud/OSD_paper/HALPERN/halpern_classify_RF_model.Rda")
plot(asv_halpern_classify)

print(asv_halpern_classify)

asv_imp <- randomForest::importance(asv_halpern_classify)
asv_imp <- data.frame(predictors = rownames(asv_imp), asv_imp)

asv_imp.sort <- arrange(asv_imp, desc(MeanDecreaseGini))
asv_imp.sort$predictors <- factor(asv_imp.sort$predictors, levels = asv_imp.sort$predictors)

# Select the top 10 predictors
asv_imp.sort <- asv_imp.sort[1:100, ]

asv_imp.sort_s <- asv_imp.sort %>%
  as_tibble() %>%
  left_join(osd2014_silva_dada2_names %>% rename(predictors = asv)) %>%
  mutate(asv_name = fct_reorder(asv_name, -MeanDecreaseGini))
# ggplot
ggplot(asv_imp.sort_s, aes(x = asv_name, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying into low or mid/hight impacted OSD samples") +
  theme_light() +
  ylab("Total decrease in node impurities (Mean decrease Gini)")

asv_otutable_prop <- transform_sample_counts(osd2014_dada2_phyloseq_alpha, function(x) x/sum(x))



f <- psmelt(asv_otutable_prop) %>%
  tbl_df %>%
  select(label, OTU, Abundance, asv_name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  dplyr::rename(asv = OTU) %>%
  right_join(asv_imp.sort_s) %>%
  filter(label %in% sites) %>%
  unique() %>%
  left_join(halpern_impact %>% select(label, median))


f <- f %>%
  mutate(label = fct_reorder(label, -median), asv_name = fct_reorder(asv_name, MeanDecreaseGini))


base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

textcol <- "grey40"

#modified ggplot
p_asv <- ggplot(f,aes(y=label,x=asv_name, fill = ((Abundance) +0.001)*100))+
  geom_tile()+
  #redrawing tiles to remove cross lines from legend
  geom_tile(colour="white",size=0.25, show.legend = FALSE)+
  #remove axis labels, add title
  labs(x="",y="",title="")+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+
  #custom breaks on x-axis
  scale_x_discrete(expand=c(0,0))+
  #custom colours for cut levels and na values
  #scale_fill_gradientn(colours =rev(c("#d53e4f","#f46d43","#fdae61",
  # "#fee08b","#e6f598","#abdda4","#ddf1da")), na.value="grey90", trans = "sqrt", labels = scales::percent) +
  viridis::scale_fill_viridis(option = "D", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
  #scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(7,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
  #mark year of vaccination
  #geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
  #equal aspect ratio x and y axis
  coord_fixed() +
  #set base size for all font elements
  theme_grey(base_size=10)+
  #theme options
  theme(
    legend.position = "bottom",
    #remove legend title
    legend.title=element_blank(),
    #remove legend margin
    legend.spacing = grid::unit(0,"cm"),
    #change legend text properties
    legend.text=element_text(colour=textcol,size=7,face="bold"),
    #change legend key height
    legend.key.height=grid::unit(0.2,"cm"),
    #set a slim legend
    legend.key.width=grid::unit(0.8,"cm"),
    #set x axis text size and colour
    axis.text.x=element_text(hjust = 1, vjust = 0.5, colour=textcol, angle = 90, size = 6),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    #set y axis text colour and adjust vertical justification
    #axis.text.y=element_text(vjust = 0.2,colour=textcol, size = 6),
    #change axis ticks thickness
    axis.ticks=element_line(size=0.4),
    #change title font, size, colour and justification
    plot.title=element_blank(),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank())

f1 <- f %>% select(label, median) %>% unique

f1$label <- factor(f1$label, levels = f1 %>% arrange(median) %>% .$label)

p1 <- ggplot( f1, aes(x = label, y = median)) +
  geom_bar(stat = "identity", size = 0.25) +
  coord_flip() +
  #remove axis labels, add title
  labs(x="",y="",title="") +
  #remove extra space
  #custom breaks on x-axis
  scale_x_discrete(expand=c(0,0)) +
  #custom colours for cut levels and na values
  #scale_fill_gradientn(colours =rev(c("#d53e4f","#f46d43","#fdae61",
  # "#fee08b","#e6f598","#abdda4","#ddf1da")), na.value="grey90", trans = "sqrt", labels = scales::percent) +
  #scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(7,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
  #mark year of vaccination
  #geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
  #equal aspect ratio x and y axis
  #coord_fixed() +
  #set base size for all font elements
  theme_grey(base_size=10)+
  #theme options
  theme(
    legend.position = "bottom",
    #remove legend title
    legend.title=element_blank(),
    #remove legend margin
    legend.spacing = grid::unit(0,"cm"),
    #change legend text properties
    legend.text=element_text(colour=textcol,size=7,face="bold"),
    #change legend key height
    legend.key.height=grid::unit(0.2,"cm"),
    #set a slim legend
    legend.key.width=grid::unit(0.8,"cm"),
    #set x axis text size and colour
    axis.text.y=element_blank(),
    axis.ticks.y = element_blank(),
    #set y axis text colour and adjust vertical justification
    axis.text.x=element_text(colour=textcol),
    #change axis ticks thickness
    axis.ticks=element_line(size=0.4),
    #change title font, size, colour and justification
    plot.title=element_blank(),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    panel.background = element_blank())

f2 <- f %>% select(asv_name, MeanDecreaseGini) %>% unique

f2$asv_name <- factor(f2$asv_name, levels = f2 %>% arrange(MeanDecreaseGini) %>% .$asv_name)

p2 <- ggplot(f2, aes(x = asv_name, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", size = 0.25) +
  #remove axis labels, add title
  labs(x="",y="",title="")+
  #remove extra space
  #custom breaks on x-axis
  scale_x_discrete(expand=c(0,0))+
  #custom colours for cut levels and na values
  #scale_fill_gradientn(colours =rev(c("#d53e4f","#f46d43","#fdae61",
  # "#fee08b","#e6f598","#abdda4","#ddf1da")), na.value="grey90", trans = "sqrt", labels = scales::percent) +
  #scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(7,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = prettyNum) +
  #mark year of vaccination
  #geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
  #equal aspect ratio x and y axis
  #set base size for all font elements
  theme_grey(base_size=10)+
  #theme options
  theme(
    legend.position = "bottom",
    #remove legend title
    legend.title=element_blank(),
    #remove legend margin
    legend.spacing = grid::unit(0,"cm"),
    #change legend text properties
    legend.text=element_text(colour=textcol,size=7,face="bold"),
    #change legend key height
    legend.key.height=grid::unit(0.2,"cm"),
    #set a slim legend
    legend.key.width=grid::unit(0.8,"cm"),
    #set x axis text size and colour
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    #set y axis text colour and adjust vertical justification
    axis.text.y=element_text(vjust = 0.2,colour=textcol),
    #change axis ticks thickness
    axis.ticks=element_line(size=0.4),
    #change title font, size, colour and justification
    plot.title=element_blank(),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    panel.background = element_blank())

ggpubr::ggarrange(p_asv, p1, p2, ncol = 1, nrow = 3, heights = c(0.5, .25, .25))





k <- f %>%
  filter(asv %in% asv_imp.sort$predictors) %>%
  left_join(halpern_impact) %>%
  rename(type = class) %>%
  select(type, Abundance, asv_name, MeanDecreaseGini, asv)

asv_wt <- ggpubr::compare_means(Abundance ~ type, data = k, group.by = c("asv_name"), method = "wilcox.test", p.adjust.method = "BH") %>%
  filter(p.signif != "ns") %>% arrange(p.adj) %>% left_join(k)

asv_wt %>%
  arrange(desc(MeanDecreaseGini)) %>%
  filter(p.signif != "ns") %>%
  select(asv_name, MeanDecreaseGini) %>% unique()

asv_imp.sort_s[1:25,]

asv_wt$type <- factor(asv_wt$type, levels = c("low_impacted", "impacted"))

asv_wt$asv <- factor(asv_wt$asv, levels = asv_imp.sort_s$predictors)


p <- ggplot(asv_wt %>% filter(asv_name == "asv_39"), aes(type, Abundance, fill = type, color = type)) +
  geom_boxplot(position = "dodge", width = 0.5, alpha = 0.6) +
  geom_jitter(alpha = 0.6, color = "#4A4A4A", width = 0.25, ) +
  ggpubr::stat_compare_means() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = c("#2E9598", "#A8216B")) +
  scale_color_manual(values = c("#2E9598", "#A8216B")) +
  theme(legend.position = "none")

p1 <- ggplot(asv_wt %>% filter(asv_s == "ASV1115"), aes(type, value, fill = type, color = type)) +
  geom_boxplot(position = "dodge", width = 0.5, alpha = 0.6) +
  geom_jitter(alpha = 0.6, color = "#4A4A4A", width = 0.25, ) +
  ggpubr::stat_compare_means() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = c("#2E9598", "#A8216B")) +
  scale_color_manual(values = c("#2E9598", "#A8216B")) +
  theme(legend.position = "none")


p2 <- ggplot(asv_wt %>% filter(asv_s == "ASV1543"), aes(type, value, fill = type, color = type)) +
  geom_boxplot(position = "dodge", width = 0.5, alpha = 0.6) +
  geom_jitter(alpha = 0.6, color = "#4A4A4A", width = 0.25, ) +
  ggpubr::stat_compare_means() +
  scale_y_continuous(labels = c("0.0%", "1.0%", "2.0%", "3.0%", "4.0%"), breaks = c(0, 0.01, 0.02, 0.03, 0.04)) +
  #scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = c("#2E9598", "#A8216B")) +
  scale_color_manual(values = c("#2E9598", "#A8216B")) +
  theme(legend.position = "none")

p3 <- ggplot(asv_wt %>% filter(asv_s == "ASV1463"), aes(type, value, fill = type, color = type)) +
  geom_boxplot(position = "dodge", width = 0.5, alpha = 0.6) +
  geom_jitter(alpha = 0.6, color = "#4A4A4A", width = 0.25, ) +
  ggpubr::stat_compare_means(label = "p.format") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = c("#2E9598", "#A8216B")) +
  scale_color_manual(values = c("#2E9598", "#A8216B")) +
  theme(legend.position = "none")

ggpubr::ggarrange(p3,p2,p1,p, ncol = 4, nrow = 1)
ggsave("~/Downloads/asv_bp.pdf", width = 29.39, height = 6.92, units = "cm")


tax_table(osd2014_16s_otuXsample_physeq_filt_prev_prop) %>%
  as.data.frame() %>%
  rownames_to_column(var = "old") %>%
  tbl_df() %>%
  left_join(asv_name_conv) %>%
  filter(new %in% c(as.character(asv_imp.sort_s$predictors[1:3]), as.character(asv_imp.sort_s$predictors[9]))) %>%
  select(-old, -new, -genus)

# Let's predict

# First get samples not used for training
pred_sites <- out_rescaled_2013_median_long %>% filter(!(label %in% sites)) %>% .$label %>% unique


# Trim the original matrix to contain the non-training and the modules used for
asv_otutable_pred <- as.data.frame(as(otu_table(osd2014_16s_otuXsample_physeq_filt_prev_prop), "matrix"))
asv_otutable_pred <- base::as.data.frame(asv_otutable_pred)
asv_otutable_pred <- asv_otutable_pred[pred_sites, c(as.character(asv_shared_impact_low_occ$asv), asv_impact_low$asv, asv_impact_imp$asv)]

pred_asv <- predict(halpern_classify, newdata = asv_otutable_pred)

pred_impacted_metadata_asv <- data.frame(class = pred_asv) %>% mutate(label = names(pred_asv)) %>% filter(class == "impacted") %>% left_join(osd2014_metadata)
pred_low_impacted_metadata_asv <- data.frame(class = pred_asv) %>% mutate(label = names(pred_asv)) %>% filter(class != "impacted") %>% left_join(osd2014_metadata)

library(ggmap)

### Set a range

### Get a map

tmp <- pred_impacted_metadata %>% dplyr::select(label, start_lat,start_lon, dist_coast_iso3_code) %>%
  mutate(filename = paste("~/Downloads/", label, "_", dist_coast_iso3_code, "_asv_impacted.pdf", sep = ""))


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
  mutate(filename = paste("~/Downloads/", label, "_", dist_coast_iso3_code, "_asv_low_impacted.pdf", sep = ""))


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
pred_impacted_metadata_asv_mg <- pred_impacted_metadata_asv %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label

pred_impacted_metadata_kegg_mg_l <- length(pred_impacted_metadata_kegg_mg)
pred_impacted_metadata_asv_mg_l <- length(pred_impacted_metadata_asv_mg)

length(intersect(pred_impacted_metadata_kegg_mg,pred_impacted_metadata_asv_mg))/(pred_impacted_metadata_kegg_mg_l + pred_impacted_metadata_asv_mg_l)

pred_low_impacted_metadata_kegg_mg <- pred_low_impacted_metadata_kegg %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label
pred_low_impacted_metadata_asv_mg <- pred_low_impacted_metadata_asv %>% filter(label %in% osd2014_amp_mg_intersect$label) %>% .$label

pred_low_impacted_metadata_kegg_mg_l <- length(pred_low_impacted_metadata_kegg_mg)
pred_low_impacted_metadata_asv_mg_l <- length(pred_low_impacted_metadata_asv_mg)

length(intersect(pred_low_impacted_metadata_kegg_mg,pred_low_impacted_metadata_asv_mg))/(pred_low_impacted_metadata_kegg_mg_l + pred_low_impacted_metadata_asv_mg_l)


length(as.character(osd2014_amp_mg_intersect$label)[as.character(pred_impacted_metadata_kegg$label)])
length(pred_impacted_metadata_asv$label)
length(pred_low_impacted_metadata_kegg$label)
length(pred_low_impacted_metadata_asv$label)

intersect(pred_impacted_metadata_kegg$label, pred_impacted_metadata_asv$label)
intersect(pred_low_impacted_metadata_kegg$label, pred_low_impacted_metadata_asv$label)
length(intersect(pred_impacted_metadata_kegg$label, pred_impacted_metadata_asv$label))
length(intersect(pred_low_impacted_metadata_kegg$label, pred_low_impacted_metadata_asv$label))
