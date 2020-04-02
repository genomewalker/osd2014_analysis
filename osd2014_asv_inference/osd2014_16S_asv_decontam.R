library(phyloseq)
library(decontam)
library(RPostgreSQL)
library(tidyverse)

# Load main samples phyloseq object
load(file = "osd2014_16S_asv/data/osd2014_16S_asv_physeq.Rdata")

# Load blank phyloseq object
load(file = "osd2014_16S_asv/data/osd2014_16S_asv_physeq_blanks.Rdata")


# Combine both objects
df <- bind_rows(as.data.frame(sample_data(osd2014_dada2_phyloseq)) %>% as_tibble() %>% select(label) %>% mutate(is.sample = TRUE, LibrarySize = sample_sums(osd2014_dada2_phyloseq)),
                as.data.frame(sample_data(osd2014_dada2_phyloseq_blanks)) %>% as_tibble() %>% select(label) %>% mutate(is.sample = FALSE, LibrarySize = sample_sums(osd2014_dada2_phyloseq_blanks)))
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=is.sample)) + geom_point()


intersect(colnames(seqtab.nochim), taxa_names(osd2014_dada2_phyloseq_blanks)) %>% length()

ps <- phyloseq::merge_phyloseq(osd2014_dada2_phyloseq_blanks, osd2014_dada2_phyloseq_all)

sample_data(ps)$is.neg <- ifelse(grepl("DNA", sample_data(ps)$label), TRUE, FALSE)
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

which(contamdf.prev$contaminant)

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

contamdf.prev05 %>% as_tibble(rownames = "asv") %>% inner_join(osd2014_silva_dada2_blanks)

ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$is.neg == TRUE, ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$is.neg == FALSE, ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

ps.prop <- transform_sample_counts(ps, function(x) x/sum(x))

ps.prop_contam <- psmelt(ps.prop) %>% as_tibble()

seqtab.nochim.prop <- vegan::decostand(seqtab.nochim, method = "total") %>%
  as_tibble(rownames = "Sample") %>%
  gather(key = OTU, value = Abundance, -Sample)

seqtab.nochim.prop %>% filter(OTU %in% osd2014_silva_dada2_blanks$asv) %>%
  ggplot(aes(Sample, Abundance)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("") +
  ylab("Proportion") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1))
