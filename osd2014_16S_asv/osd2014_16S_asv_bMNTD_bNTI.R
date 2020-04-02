## change to the directory on your computer that contains the OTU table and phylogeny
## note that the 'slash' needs to be changed to a forward slash like this /
setwd("C:/Users/steg815/Desktop/Stegen_PNNL/")
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_with_phylo.Rdata")
## load this library
## if not already installed, use install.packages('picante')
library(picante)

## read in OTU table

otu <- t(phyloseq:::veganifyOTU(osd2014_dada2_phyloseq_beta_vst))

phylo <- phy_tree(osd2014_dada2_phyloseq_beta_vst);

## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu <- match.phylo.data(phylo, otu);

## calculate empirical betaMNTD

beta.mntd.weighted <- as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));

beta.reps <- 999; # number of randomizations

rand.weighted.bMNTD.comp <- array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));

rand_mntd <- function(X) comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F)

zz <- pbmcapply::pbmclapply(1:beta.reps,rand_mntd, mc.cores = 4)

for (i in 1:beta.reps){
  rand.weighted.bMNTD.comp[,,i] <- as.matrix(zz[[i]])
}


weighted.bNTI <- matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    rand.vals <- rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] <- (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
  };
};

rownames(weighted.bNTI) <- colnames(match.phylo.otu$data);
colnames(weighted.bNTI) <- colnames(match.phylo.otu$data);

hist(weighted.bNTI)
save.image("osd2014_16S_asv/data/osd2014_16S_bMNTD_bNTI_results.Rdata")
