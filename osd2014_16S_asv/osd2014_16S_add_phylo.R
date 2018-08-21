library(phyloseq)
library(ape)
library(picante)

load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata", verbose = TRUE)

# Read phylogenetic tree --------------------------------------------------
osd2014_dada2_tree <- read.tree("osd2014_16S_asv/data/qiime2_placement/exported-tree/tree.nwk")

osd2014_dada2_phyloseq_beta_comm <- match.phylo.comm(osd2014_dada2_tree, (otu_table(osd2014_dada2_phyloseq_beta_scaled)))
phy_tree(osd2014_dada2_phyloseq_beta) <- osd2014_dada2_phyloseq_beta_comm$phy
phy_tree(osd2014_dada2_phyloseq_beta_scaled) <- osd2014_dada2_phyloseq_beta_comm$phy
phy_tree(osd2014_dada2_phyloseq_beta_css) <- osd2014_dada2_phyloseq_beta_comm$phy
phy_tree(osd2014_dada2_phyloseq_beta_vst) <- osd2014_dada2_phyloseq_beta_comm$phy
phy_tree(osd2014_dada2_phyloseq_beta_norm) <- osd2014_dada2_phyloseq_beta_comm$phy


osd2014_dada2_phyloseq_alpha_comm <- match.phylo.comm(osd2014_dada2_tree, (otu_table(osd2014_dada2_phyloseq_alpha_scaled)))
phy_tree(osd2014_dada2_phyloseq_alpha) <- osd2014_dada2_phyloseq_alpha_comm$phy
phy_tree(osd2014_dada2_phyloseq_alpha_scaled) <- osd2014_dada2_phyloseq_alpha_comm$phy
phy_tree(osd2014_dada2_phyloseq_alpha_css) <- osd2014_dada2_phyloseq_alpha_comm$phy
phy_tree(osd2014_dada2_phyloseq_alpha_vst) <- osd2014_dada2_phyloseq_alpha_comm$phy
phy_tree(osd2014_dada2_phyloseq_alpha_norm) <- osd2014_dada2_phyloseq_alpha_comm$phy

save.image("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_with_phylo.Rdata")
