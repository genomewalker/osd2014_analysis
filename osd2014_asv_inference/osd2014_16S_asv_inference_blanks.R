library(dada2); packageVersion("dada2")
library(tidyverse)

# File parsing
wdir <- "osd2014_asv_inference/data/dada_2014_blanks/"
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

taxa <- assignTaxonomy(seqtab.nochim, "http://osd2014.metagenomics.eu/osd2014_asv_inference/data/silva_species_assignment_v132.fa.gz", multithread=TRUE, tryRC = TRUE)
taxa <- addSpecies(taxa, "http://osd2014.metagenomics.eu/osd2014_asv_inference/data/silva_species_assignment_v132.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

taxa_df <- as(taxa, "matrix") %>% as.data.frame()

taxa_df <- taxa %>%
  as_tibble(rownames = "asv") %>%
  mutate(asv_name = paste("asv", row_number(), sep = "_"))

cdata::

osd2014_dada2_phyloseq_blanks <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(as.matrix(taxa)), sample_data(cdata))
