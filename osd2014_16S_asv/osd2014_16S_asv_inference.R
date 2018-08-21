library(dada2); packageVersion("dada2")
library(tidyverse)

# File parsing
wdir <- "/Volumes/Extra/dada_2014/"
setwd(wdir)
pathF <-  file.path(wdir, "files")
pathR <- file.path(wdir, "files")

filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...

fastqFs <- sort(list.files(pathF, pattern="R1.fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="R2.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

sample.names <- sapply(strsplit(fastqFs, "_R1.fastq.gz"), `[`, 1)

filtpathF <- file.path(filtpathF, gsub("R1.fastq.gz","R1.filt.fastq.gz", fastqFs))
filtpathR <- file.path(filtpathR, gsub("R2.fastq.gz","R2.filt.fastq.gz", fastqRs))

fastqFs <- file.path(pathF, fastqFs)
fastqRs <- file.path(pathR, fastqRs)

qpF <- plotQualityProfile(fastqFs[1:5], aggregate = TRUE)
qpR <- plotQualityProfile(fastqRs[1:5], aggregate = TRUE)

# Filter and Trim reads (TruncLen and maxEE based on previous adjustments with reads.out mostly >= 70%)
ft <- filterAndTrim(fwd=fastqFs, filt=filtpathF,
                    rev=fastqRs, filt.rev=filtpathR,
                    truncLen=c(220,200), maxEE=c(8,8), maxN=0, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE, multithread=TRUE)

ft %>% tbl_df %>% mutate(prop = 100 * (reads.out/reads.in)) %>% .$prop %>% hist()

# Learn forward error rates
errF <- learnErrors(filtpathF, multithread = TRUE, randomize = TRUE)
# Learn reverse error rates
errR <- learnErrors(filtpathR, multithread = TRUE, randomize = TRUE)

ErrorPlotF <- plotErrors(errF, nominalQ = TRUE)
ErrorPlotR <- plotErrors(errR, nominalQ = TRUE)

#Dereplication
derepFs <- derepFastq(filtpathF, verbose=TRUE)
derepRs <- derepFastq(filtpathR, verbose=TRUE)


#Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread = TRUE, pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread = TRUE, pool="pseudo")

dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, returnRejects = TRUE)
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads throughout the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(ft, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

track_long <- as.data.frame(track) %>%
  rownames_to_column(var = "sample") %>%
  gather(variable, value, -sample) %>%
  tbl_df()

track_long$variable <- factor(track_long$variable, levels = colnames(track))

ggplot(track_long, aes(variable, value, fill = variable, color = variable)) +
  #geom_bar(stat = "identity", position = "dodge") +
  geom_jitter() +
  scale_fill_brewer(palette="Paired") +
  #facet_wrap(~sample, scales = "free_y") +
  xlab("DADA2 steps") +
  ylab("Number of sequences") +
  theme_bw()

taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)
taxa <- addSpecies(taxa, "~/Downloads/silva_species_assignment_v132.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

taxa_df <- as(taxa, "data.frame")

taxa_df <- taxa %>%
  as_tibble(rownames = "asv") %>%
  mutate(asv_name = paste("asv", row_number(), sep = "_"))

# read SILVA taxonomy
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

st_100_order_terrestrial <- tbl(my_db, "osd2014_st_order_terrestrial") %>%
  collect(n = Inf)

osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)

osd2014_cdata <- as.data.frame(osd2014_cdata)
rownames(osd2014_cdata) <- osd2014_cdata$label

osd2014_silva_dada2 <- tbl(my_db, "osd2014_silva_dada2") %>%
  collect(n = Inf)

osd2014_silva_dada2_blanks <- tbl(my_db, "osd2014_silva_dada2_blanks") %>%
  collect(n = Inf)

osd2014_silva_dada2_filt <- osd2014_silva_dada2 %>%
  filter(Kingdom != "Eukaryota" | is.na(Kingdom)) %>%
  filter(Family != "Mitochondria" | is.na(Family)) %>%
  filter(Order != "Chloroplast" | is.na(Order)) %>%
  filter(!(asv %in% osd2014_silva_dada2_blanks$asv)) %>%
  filter(!(is.na(Phylum)))

osd2014_taxa <- as.data.frame(osd2014_silva_dada2)
rownames(osd2014_taxa) <- osd2014_taxa$asv
osd2014_taxa$asv <- NULL

osd2014_taxa_filt <- as.data.frame(osd2014_silva_dada2_filt)
rownames(osd2014_taxa_filt) <- osd2014_taxa_filt$asv
osd2014_taxa_filt$asv <- NULL

row.names(seqtab.nochim) <- gsub("_R1.filt.fastq.gz", "", row.names(seqtab.nochim))

osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)

osd2014_dada2_phyloseq_all <- phyloseq(otu_table(seqtab.nochim[,osd2014_silva_dada2$asv], taxa_are_rows=FALSE), tax_table(as.matrix(osd2014_taxa)), sample_data(osd2014_cdata))
osd2014_dada2_phyloseq <- phyloseq(otu_table(seqtab.nochim[,osd2014_silva_dada2_filt$asv], taxa_are_rows=FALSE), tax_table(as.matrix(osd2014_taxa)), sample_data(osd2014_cdata))

sample_sums(osd2014_dada2_phyloseq) %>% sort

osd2014_dada2_phyloseq_all_prop <- transform_sample_counts(osd2014_dada2_phyloseq_all, function(X) X/sum(X))

euks <- subset_taxa(osd2014_dada2_phyloseq_all_prop, Kingdom == "Eukaryota")
chlr_prop <- subset_taxa(osd2014_dada2_phyloseq_all_prop, Order == "Chloroplast")
chlr <- subset_taxa(osd2014_dada2_phyloseq_all, Order == "Chloroplast")
mit <- subset_taxa(osd2014_dada2_phyloseq_all_prop, Family == "Mitochondria")

((sample_sums(chlr) %>% sum())/(sample_sums(osd2014_dada2_phyloseq) %>% sum()))*100

chlr_prop_long <- psmelt(chlr_prop) %>% tbl_df() %>% select(Sample, Abundance)

chlr_prop_long %>%
  group_by(Sample) %>%
  summarise(prop = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Sample = fct_relevel(Sample, st_100_order_terrestrial$label)) %>%
  ggplot(aes(Sample, prop)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("") +
  ylab("Proportion") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1))

write_tsv(osd2014_silva_dada2 %>% select(asv_name, asv), path = "osd2014_16S_asv/data/osd2014_16S_asv_sequences.tsv", col_names = FALSE)

save(osd2014_dada2_phyloseq, osd2014_dada2_phyloseq_all, file = "osd2014_16S_asv/data/osd2014_16S_asv_physeq.Rdata")
save(ft, errF, errR, derepFs, derepRs, dadaFs, dadaRs, mergers, seqtab, seqtab.nochim, taxa,file = "osd2014_16S_asv/data/osd2014_16S_asv_inference.Rdata")

osd2014_silva_dada2 <- osd2014_silva_dada2 %>% rename(asv = X7) %>%
  mutate(asv_name = paste("asv", row_number(), sep = "_"))

library(RPostgreSQL)  # loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_silva_dada2"), value=osd2014_silva_dada2,overwrite = TRUE, row.names = FALSE)



# Get non-merged reads ----------------------------------------------------
library(tidyverse)
library(phyloseq)
library(dada2)
library(purrr)
load("osd2014_16S_asv/data/osd2014_16S_asv_inference.Rdata")

concat <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, justConcatenate=TRUE)


get_nonmerged <- function(X, merg = merg, conc = conc){
m <- merg[[X]]
cn <- conc[[X]]
bind_cols(m, cn) %>% as_tibble() %>%
  filter(accept == FALSE, nmatch == 0, nmismatch == 0, nindel == 0) %>%
  mutate(sequence = sequence1, accept = TRUE) %>%
  select(sequence, abundance, forward, reverse, nmatch, nmismatch, nindel, prefer, accept) %>%
  as.data.frame()
}

get_nonmerged_fwd <- function(X, merg = merg, conc = conc){
  m <- merg[[X]]
  cn <- conc[[X]]
  bind_cols(m, cn) %>% as_tibble() %>%
    filter(accept == FALSE, nmatch == 0, nmismatch == 0, nindel == 0) %>%
    separate(sequence1, into = c("R1", "R2"), sep = "NNNNNNNNNN", remove = FALSE) %>%
    mutate(sequence = R1, accept = TRUE) %>%
    select(sequence, abundance, forward, reverse, nmatch, nmismatch, nindel, prefer, accept) %>%
    as.data.frame()
}

get_nonmerged_rev <- function(X, merg = merg, conc = conc){
  m <- merg[[X]]
  cn <- conc[[X]]
  bind_cols(m, cn) %>% as_tibble() %>%
    filter(accept == FALSE, nmatch == 0, nmismatch == 0, nindel == 0) %>%
    separate(sequence1, into = c("R1", "R2"), sep = "NNNNNNNNNN", remove = FALSE) %>%
    mutate(sequence = R2, accept = TRUE) %>%
    select(sequence, abundance, forward, reverse, nmatch, nmismatch, nindel, prefer, accept) %>%
    as.data.frame()
}

test_f <- lapply(names(mergers), get_nonmerged_fwd, merg = mergers, conc = concat)
test_r <- lapply(names(mergers), get_nonmerged_rev, merg = mergers, conc = concat)

names(test_r) <- names(mergers)
#Construct sequence table
seqtab_nm <- makeSequenceTable(test_r)
dim(seqtab_nm)
table(nchar(getSequences(seqtab_nm)))

#Remove chimeras
seqtab_nm.nochim <- removeBimeraDenovo(seqtab_nm, method = "consensus", multithread = 2, verbose = TRUE)
dim(seqtab_nm.nochim)
sum(seqtab_nm.nochim)/sum(seqtab_nm)

#Track reads throughout the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(ft, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

track_long <- as.data.frame(track) %>%
  rownames_to_column(var = "sample") %>%
  gather(variable, value, -sample) %>%
  tbl_df()

track_long$variable <- factor(track_long$variable, levels = colnames(track))

ggplot(track_long, aes(variable, value, fill = variable, color = variable)) +
  #geom_bar(stat = "identity", position = "dodge") +
  geom_jitter() +
  scale_fill_brewer(palette="Paired") +
  #facet_wrap(~sample, scales = "free_y") +
  xlab("DADA2 steps") +
  ylab("Number of sequences") +
  theme_bw()

taxa_nm <- assignTaxonomy(seqtab_nm.nochim, "~/Downloads/silva_nr_v132_train_set.fa.gz", multithread=2, tryRC = TRUE)
taxa_nm <- addSpecies(taxa_nm, "~/Downloads/silva_species_assignment_v132.fa.gz")
taxa_nm.print <- taxa_nm # Removing sequence rownames for display only
rownames(taxa_nm.print) <- NULL
head(taxa_nm.print)
row.names(seqtab_nm.nochim) <- gsub("_R1.filt.fastq.gz", "", row.names(seqtab_nm.nochim))

osd2014_dada2_phyloseq_nm <- phyloseq(otu_table(seqtab_nm.nochim, taxa_are_rows=FALSE), tax_table(as.matrix(taxa_nm)), sample_data(osd2014_cdata))

