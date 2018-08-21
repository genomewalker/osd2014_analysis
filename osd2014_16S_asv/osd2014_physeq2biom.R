library(phyloseq)
library(biomformat)

load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_with_phylo.Rdata", verbose = TRUE)

osd2014_alpha_biom <- biomformat::make_biom(
  data = as(t(otu_table(osd2014_dada2_phyloseq_alpha, taxa_are_rows = FALSE)), "matrix"),
  matrix_element_type = "int"
)
write_biom(osd2014_alpha_biom, biom_file = "osd2014_16S_asv/data/osd2014_16S_asv_physeq_alpha_tmp.biom")
sdata <- as(sample_data(osd2014_dada2_phyloseq_alpha), "data.frame") %>% add_column(., `#SampleID` = sample_names(osd2014_dada2_phyloseq_alpha), .before = 1)
write.table(sdata, col.names = TRUE, quote = FALSE,  row.names = FALSE, file = "osd2014_16S_asv/data/osd2014_16S_asv_physeq_alpha_sdata.txt", sep = "\t")

tdata <- as(tax_table(osd2014_dada2_phyloseq_alpha), "matrix") %>% as.data.frame() %>% add_column(., `#OTUID` = taxa_names(osd2014_dada2_phyloseq_alpha), .before = 1)
write.table(tdata, col.names = TRUE, quote = FALSE,  row.names = FALSE, file = "osd2014_16S_asv/data/osd2014_16S_asv_physeq_alpha_tdata.txt", sep = "\t")

system("biom add-metadata -i osd2014_16S_asv/data/osd2014_16S_asv_physeq_alpha_tmp.biom -o osd2014_16S_asv/data/osd2014_16S_asv_physeq_alpha.biom --sample-metadata-fp osd2014_16S_asv/data/osd2014_16S_asv_physeq_alpha_sdata.txt --observation-metadata-fp osd2014_16S_asv/data/osd2014_16S_asv_physeq_alpha_tdata.txt")
system("rm osd2014_16S_asv/data/osd2014_16S_asv_physeq_alpha_tmp.biom osd2014_16S_asv/data/osd2014_16S_asv_physeq_alpha_sdata.txt osd2014_16S_asv/data/osd2014_16S_asv_physeq_alpha_tdata.txt")


osd2014_beta_biom <- biomformat::make_biom(
  data = as(t(otu_table(osd2014_dada2_phyloseq_beta, taxa_are_rows = FALSE)), "matrix"),
  matrix_element_type = "int"
)
write_biom(osd2014_beta_biom, biom_file = "osd2014_16S_asv/data/osd2014_16S_asv_physeq_beta_tmp.biom")
sdata <- as(sample_data(osd2014_dada2_phyloseq_beta), "data.frame") %>% add_column(., `#SampleID` = sample_names(osd2014_dada2_phyloseq_beta), .before = 1)
write.table(sdata, col.names = TRUE, quote = FALSE,  row.names = FALSE, file = "osd2014_16S_asv/data/osd2014_16S_asv_physeq_beta_sdata.txt", sep = "\t")

tdata <- as(tax_table(osd2014_dada2_phyloseq_beta), "matrix") %>% as.data.frame() %>% add_column(., `#OTUID` = taxa_names(osd2014_dada2_phyloseq_beta), .before = 1)
write.table(tdata, col.names = TRUE, quote = FALSE,  row.names = FALSE, file = "osd2014_16S_asv/data/osd2014_16S_asv_physeq_beta_tdata.txt", sep = "\t")

system("biom add-metadata -i osd2014_16S_asv/data/osd2014_16S_asv_physeq_beta_tmp.biom -o osd2014_16S_asv/data/osd2014_16S_asv_physeq_beta.biom --sample-metadata-fp osd2014_16S_asv/data/osd2014_16S_asv_physeq_beta_sdata.txt --observation-metadata-fp osd2014_16S_asv/data/osd2014_16S_asv_physeq_beta_tdata.txt")
system("rm osd2014_16S_asv/data/osd2014_16S_asv_physeq_beta_tmp.biom osd2014_16S_asv/data/osd2014_16S_asv_physeq_beta_sdata.txt osd2014_16S_asv/data/osd2014_16S_asv_physeq_beta_tdata.txt")

