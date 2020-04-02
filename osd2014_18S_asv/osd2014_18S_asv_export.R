library(phyloseq)

load(file = "osd2014_18S_asv/data/osd2014_18S_asv_physeq.Rdata", verbose = TRUE)

asv_headers <- paste0(">",taxa_names(osd2014_dada2_phyloseq))
asv_seqs <- taxa_names(osd2014_dada2_phyloseq)

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "osd2014_18S_asv/data/osd2014_18S_ASVs.fasta")
