library(phyloseq)
library(ape)
library(picante)

# BEGIN: WARNING!!!! -------------------------------------------------------------
# You can access to the data used in this analysis in several ways:
# 1. You have a copy of the PostgreSQL DB
# 2. You downloaded the .Rdata files from http://osd2014.metagenomics.eu/ and placed them
#    in the data folder
# 3. You can load the files remotely, it might take a while when the file is very large
# END: WARNING!!!! -------------------------------------------------------------


# BEGIN: Load necessary data -----------------------------------------------------
# Use if you have the postgres DB in place
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")

# If downloaded file at osd2014_16S_asv/data/ use:
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata", verbose = TRUE)
osd2014_dada2_tree <- read.tree("osd2014_16S_asv/data/qiime2_placement/exported-tree/tree.nwk")

# If remote use
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata"), verbose = TRUE)
osd2014_dada2_tree <- read.tree(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/qiime2_placement/exported-tree/tree.nwk"))
# END: Load necessary data -----------------------------------------------------


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

# BEGIN: Save objects ------------------------------------------------------------
# WARNING!!! You might not want to run this code --------------------------
save.image("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_with_phylo.Rdata")
# END: Save objects ------------------------------------------------------------
