# Prepare data for PICRUST ------------------------------------------------

# Modified from https://github.com/vmaffei/dada2_to_picrust ---------------


# Dependencies: ShortRead & biom
library(ShortRead)
library(biom) # note: use Joey's biom latest dev version; library(devtools); install_github("joey711/biom");
library(phyloseq)
library(ggpol)
# 1) Make study db
# grab study seqs


# BEGIN: WARNING!!!! -------------------------------------------------------------
# You can access to the data used in this analysis in several ways:
# 1. You have a copy of the PostgreSQL DB
# 2. You downloaded the .Rdata files from http://osd2014.metagenomics.eu/ and placed them
#    in the data folder
# 3. You can load the files remotely, it might take a while when the file is very large
# END: WARNING!!!! -------------------------------------------------------------


# Load necessary data -----------------------------------------------------
# Use if you have the postgres DB in place
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_read_ko20140317_abun <- tbl(my_db, "osd2014_read_ko20140317_abun") %>%
  collect(n = Inf)


# If downloaded file at osd2014_16S_asv/data/ use:
load(file = "osd2014_16S_asv/data/osd2014_16S_seqtabnochim.Rdata", verbose = TRUE)
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq.Rdata", verbose = TRUE)
load("osd2014_shotgun/data/osd2014_read_ko20140317_abun.Rdata", verbose = TRUE)

# If remote use
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_seqtabnochim.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_asv_physeq.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_shotgun/data/osd2014_read_ko20140317_abun.Rdata"), verbose = TRUE)

# Load necessary data -----------------------------------------------------


seqs_study <- colnames(seqtab.nochim)
ids_study <- paste("study", 1:ncol(seqtab.nochim), sep = "_")
# merge db and study seqs
db_out <- data.frame(ids=ids_study,seqs=seqs_study,count=colSums(seqtab.nochim))
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
# write study fasta for filtering
writeFasta(fasta, file = "osd2014_16S_asv/data/gg_13_5_study_db.fasta.pre")
# filter sequences that diverge from gg_13_5 by 97%
# depending on how well greengenes covers your study sequences, consider reducing 97% to 70 or 50%
# You will need vsearch installed
system('vsearch --usearch_global osd2014_16S_asv/data/gg_13_5_study_db.fasta.pre --db osd2014_16S_asv/data/97_otus.fasta --matched osd2014_16S_asv/data/gg_13_5_study_db.fasta --id 0.90')
id_filtered <- as.character(id(readFasta("osd2014_16S_asv/data/gg_13_5_study_db.fasta")))
db_out_filt <- db_out[db_out$ids%in%id_filtered,]
seqtab_biom <- t(seqtab.nochim)
# 2) output seq variant count data as biom;
# subset seqtab and output sample count biom
seqtab_biom <- seqtab_biom[rownames(seqtab_biom)%in%db_out_filt$seqs,]
rownames(seqtab_biom) <- db_out_filt[db_out_filt$seqs%in%rownames(seqtab_biom),"ids"]
biom_object <- biom::make_biom(data = seqtab_biom)
biom::write_biom(biom_object, biom_file = "osd2014_16S_asv/data/sample_counts.biom")
# create final study db
system('cat osd2014_16S_asv/data/97_otus.fasta >> osd2014_16S_asv/data/gg_13_5_study_db.fasta')


# Picrust part ------------------------------------------------------------

# align w/ pynast using QIIME scripts; the included options lessen alignment restrictions to prevent alignment failure
# minimum sequence length set by -e
# alignment runtime greatly reduced by parallelization: parallel_align_seqs_pynast.py -O and # of cores
parallel_align_seqs_pynast.py -e 90 -p 0.1 -i ./genome_prediction/gg_13_5_study_db.fasta -o ./genome_prediction/gg_13_5_study_db.fasta.aligned -O 45
# filter alignment with default settings; consider lane filtering by entropy using -e and a low entropy value of ~0.01-0.02
# note: FastTree and/or PICRUSt produce weird errors (segfaults) if -e filters too many lanes
filter_alignment.py -i ./genome_prediction/gg_13_5_study_db.fasta.aligned/gg_13_5_study_db_aligned.fasta -o ./genome_prediction/gg_13_5_study_db.fasta.aligned.filtered/
# build tree with fasttree; options are taken from greengenes 13_5 readme notes
# tree building runtime greatly reduced by parallelization: use FastTreeMP w/ same options instead of FastTree
FastTreeMP -nt -gamma -fastest -no2nd -spr 4 ./genome_prediction/gg_13_5_study_db.fasta.aligned.filtered/gg_13_5_study_db_aligned_pfiltered.fasta > ./genome_prediction/study_tree.tree


# format 16S copy number data
format_tree_and_trait_table.py -t ./genome_prediction/study_tree.tree -i gg_16S_counts.tab -o ./genome_prediction/format/16S/
# format kegg IMG data
format_tree_and_trait_table.py -t ./genome_prediction/study_tree.tree -i gg_ko_counts.tab -o ./genome_prediction/format/KEGG/
# perform ancestral state reconstruction
ancestral_state_reconstruction.py -i ./genome_prediction/format/16S/trait_table.tab -t ./genome_prediction/format/16S/pruned_tree.newick -o ./genome_prediction/asr/16S_asr_counts.tab -c ./genome_prediction/asr/asr_ci_16S.tab
ancestral_state_reconstruction.py -i ./genome_prediction/format/KEGG/trait_table.tab -t ./genome_prediction/format/KEGG/pruned_tree.newick -o ./genome_prediction/asr/KEGG_asr_counts.tab -c ./genome_prediction/asr/asr_ci_KEGG.tab
# collect study sequence ids for predict_traits.py -l (greatly reduces runtime)
# convert biom to tsv using biom-format
biom convert -i sample_counts.biom -o sample_counts.tab --to-tsv
# predict traits
predict_traits.py -i ./genome_prediction/format/16S/trait_table.tab -t ./genome_prediction/format/16S/reference_tree.newick -r ./genome_prediction/asr/16S_asr_counts.tab -o ./genome_prediction/predict_traits/16S_precalculated.tab -a -c ./genome_prediction/asr/asr_ci_16S.tab -l sample_counts.tab
predict_traits.py -i ./genome_prediction/format/KEGG/trait_table.tab -t ./genome_prediction/format/KEGG/reference_tree.newick -r ./genome_prediction/asr/KEGG_asr_counts.tab -o ./genome_prediction/predict_traits/ko_precalculated.tab -a -c ./genome_prediction/asr/asr_ci_KEGG.tab -l sample_counts.tab
# add KEGG metadata
cat kegg_meta >> ./genome_prediction/predict_traits/ko_precalculated.tab

# run PICRUSt
normalize_by_copy_number.py -i sample_counts.biom -o norm_counts.biom -c ./genome_prediction/predict_traits/16S_precalculated.tab
predict_metagenomes.py -i norm_counts.biom -o meta_counts_asr.biom -c ./genome_prediction/predict_traits/ko_precalculated.tab
# optional: agglomerate counts by KEGG pathway level
categorize_by_function.py -i meta_counts_asr.biom -o cat_2_counts_asr.biom -c KEGG_Pathways -l 2



# Process PICRUST results -------------------------------------------------
library(biomformat)
library(phyloseq)
library(tidyverse)
library(Hmisc)
source("https://gist.githubusercontent.com/genomewalker/8abc47a044f0ac98b7392bbef8afcffa/raw/f4cc661bc3c7db5998b89d01006f710b1eb936d0/geom_flat_violin.R")

# Read results
# If downloaded file at osd2014_16S_asv/data/ use:
ko_picrust <- as.matrix(biom_data(read_biom("osd2014_16S_asv/data/meta_counts_asr.biom")))

# If remote use
download.file("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/meta_counts_asr.biom", destfile = "osd2014_16S_asv/data/meta_counts_asr.biom")
ko_picrust <- as.matrix(biom_data(read_biom("osd2014_16S_asv/data/meta_counts_asr.biom")))



ko_picrust <- phyloseq(otu_table = otu_table(ko_picrust, taxa_are_rows = T))

ko_picrust_prop <- transform_sample_counts(ko_picrust, function(x)x/sum(x))
ko_picrust_prop_long <- psmelt(ko_picrust_prop) %>% as_tibble() %>%
  arrange(Sample,OTU) %>%
  rename(ko = OTU, label = Sample, prop_picrust = Abundance)




# Run TAX4fun -------------------------------------------------------------

library(themetagenomics)

osd2014_t4fun <- list()
osd2014_t4fun$ABUND <- as(otu_table(osd2014_dada2_phyloseq_all, taxa_are_rows = TRUE), "matrix")
colnames(osd2014_t4fun$ABUND) <- paste("asv", 1:ncol(osd2014_t4fun$ABUND), sep = "_")
osd2014_t4fun$TAX <- as(tax_table(osd2014_dada2_phyloseq_all), "matrix")
rownames(osd2014_t4fun$TAX) <- paste("asv", 1:ncol(osd2014_t4fun$ABUND), sep = "_")

tmp <- tempdir()
download_ref(tmp,reference='silva_ko',overwrite=FALSE)

system.time(FUNCTIONS <- t4f(osd2014_t4fun$ABUND,rows_are_taxa=FALSE,tax_table=osd2014_t4fun$TAX[,1:6],
                             reference_path=tmp,type='uproc',short=TRUE,
                             cn_normalize=TRUE,sample_normalize=TRUE,drop=TRUE))


ko_tax4fun_prop_long <- FUNCTIONS$fxn_table %>% as.data.frame() %>% rownames_to_column("label") %>% gather(ko, prop_t4f, -label) %>% as_tibble()




osd2014_read_ko20140317_abun_prop <- osd2014_read_ko20140317_abun %>%
  group_by(label) %>%
  mutate(prop = abun/sum(abun)) %>%
  select(-abun) %>%
  arrange(label, ko)

# Calculate correlations and plot WGS vs PICRUST
osd2014_read_ko20140317_abun_prop %>%
  inner_join(ko_picrust_prop_long) %>%
  ungroup() %>%
  group_by(label) %>%
  do(rcorr(.$prop, .$prop_picrust, "spearman") %>% broom::tidy()) %>%
  ungroup() %>%
  ggplot(aes("A", y = estimate)) +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, errorbar.draw = TRUE, width = 0.2, jitter.colour = "black", fill = "#103C54", alpha = 0.7) +
  theme(legend.position = "none") +
  ylab(expression(Spearman~rho)) +
  xlab("") +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


# Calculate correlations and plot WGS vs TAX4FUN
osd2014_read_ko20140317_abun_prop %>%
  inner_join(ko_tax4fun_prop_long) %>%
  ungroup() %>%
  group_by(label) %>%
  do(rcorr(.$prop, .$prop_t4f, "spearman") %>% broom::tidy()) %>%
  ungroup() %>%
  ggplot(aes("A", y = estimate)) +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, errorbar.draw = TRUE, width = 0.2, jitter.colour = "black", fill = "#103C54", alpha = 0.7) +
  theme(legend.position = "none") +
  ylab(expression(Spearman~rho)) +
  xlab("") +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Calculate correlations and plot TAX4FUN vs PICRUST
ko_picrust_prop_long %>%
  inner_join(ko_tax4fun_prop_long) %>%
  ungroup() %>%
  group_by(label) %>%
  do(rcorr(.$prop_t4f, .$prop_picrust, "spearman") %>% broom::tidy()) %>%
  ungroup() %>%
  ggplot(aes("A", y = estimate)) +
  # geom_flat_violin(scale = "count", trim = FALSE) +
  # stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
  #              geom = "pointrange", position = position_nudge(0.05)) +
  # geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "down",
  #              position = position_nudge(-0.025)) +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, errorbar.draw = TRUE, width = 0.2, jitter.colour = "black", fill = "#103C54", alpha = 0.7) +
  theme(legend.position = "none") +
  ylab(expression(Spearman~rho)) +
  xlab("") +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
