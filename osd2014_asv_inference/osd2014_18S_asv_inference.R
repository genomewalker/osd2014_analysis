library(dada2); packageVersion("dada2")
library(tidyverse)

# File parsing
wdir <- "osd2014_asv_inference/data/dada2_2014"
wdir <- "/Volumes/Extra/dada2_18S/"

setwd(wdir)
pathF <-  file.path(wdir, "files")
pathR <- file.path(wdir, "files")

filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...

fastqFs <- sort(list.files(pathF, pattern="_R1_18S_raw.fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="_R2_18S_raw.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

sample.names <- sapply(strsplit(fastqFs, "_R1_18S_raw.fastq.gz"), `[`, 1)

filtpathF <- file.path(filtpathF, gsub("R1_18S_raw.fastq.gz","R1.filt.fastq.gz", fastqFs))
filtpathR <- file.path(filtpathR, gsub("R2_18S_raw.fastq.gz","R2.filt.fastq.gz", fastqRs))

fastqFs <- file.path(pathF, fastqFs)
fastqRs <- file.path(pathR, fastqRs)

qpF <- plotQualityProfile(fastqFs, aggregate = TRUE)
qpR <- plotQualityProfile(fastqRs, aggregate = TRUE)

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
  ggparl::geom_boxjitter(outlier.color = NA, jitter.shape = 21, errorbar.draw = TRUE, width = 0.2, alpha = 0.7, color = "#666666") +
  scale_fill_brewer(palette="Paired") +
  scale_color_brewer(palette="Paired") +
  #facet_wrap(~sample, scales = "free_y") +
  xlab("DADA2 steps") +
  ylab("Number of sequences") +
  scale_y_continuous(labels = scales::comma) +
  theme_bw()

taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads//pr2_version_4.10.0_dada2.fasta.gz", multithread=TRUE, tryRC = TRUE,
                       taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"))

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


osd2014_silva_dada2_filt <- taxa_df %>%
  filter(Kingdom == "Eukaryota")

osd2014_taxa <- as.data.frame(osd2014_silva_dada2_filt)
rownames(osd2014_taxa) <- osd2014_taxa$asv
osd2014_taxa$asv <- NULL

osd2014_taxa_filt <- as.data.frame(osd2014_silva_dada2_filt)
rownames(osd2014_taxa_filt) <- osd2014_taxa_filt$asv
osd2014_taxa_filt$asv <- NULL

row.names(seqtab.nochim) <- gsub("_R1.filt.fastq.gz", "", row.names(seqtab.nochim))

osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)

osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)

osd2014_cdata <- as.data.frame(osd2014_cdata)
rownames(osd2014_cdata) <- osd2014_cdata$label

osd2014_dada2_phyloseq <- phyloseq(otu_table(seqtab.nochim[,osd2014_silva_dada2_filt$asv], taxa_are_rows=FALSE), tax_table(as.matrix(osd2014_taxa)), sample_data(osd2014_cdata))

sample_sums(osd2014_dada2_phyloseq) %>% sort

write_tsv(osd2014_silva_dada2_filt %>% select(asv_name, asv), path = "osd2014_18S_asv/data/osd2014_18S_asv_sequences.tsv", col_names = FALSE)

save(osd2014_dada2_phyloseq, file = "osd2014_18S_asv/data/osd2014_18S_asv_physeq.Rdata")
save(ft, errF, errR, derepFs, derepRs, dadaFs, dadaRs, mergers, seqtab, seqtab.nochim, taxa,file = "osd2014_18S_asv/data/osd2014_18S_asv_inference.Rdata")



library(RPostgreSQL)  # loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")  # creates a connection to the postgres database  # note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "osd_analysis", host = "localhost", port = 5432)
dbWriteTable(con, c("osd_analysis", "osd2014_18S_pr2_dada2"), value=osd2014_silva_dada2_filt,overwrite = TRUE, row.names = FALSE)
