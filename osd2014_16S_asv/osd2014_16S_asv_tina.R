source("osd2014_16S_asv/lib/functions_com_sim.R")

# Load Packages
library("tidyverse")
library("foreach")
library("doMC")
library("iterators")
library("doParallel")
library("parallel")
library("Matrix")
library("bigmemory")
library("biganalytics")
library("gRbase")
library("gplots")
library("ggplot2")
library("grid")
library("gridExtra")
library("data.table")
library("plyr")
library("ape")
library("phyloseq")
library("vegan")
library("RColorBrewer")



# BEGIN: WARNING!!!! -------------------------------------------------------------
# You can access to the data used in this analysis in several ways:
# 1. You have a copy of the PostgreSQL DB
# 2. You downloaded the .Rdata files from http://osd2014.metagenomics.eu/ and placed them
#    in the data folder
# 3. You can load the files remotely, it might take a while when the file is very large
# END: WARNING!!!! -------------------------------------------------------------


# BEGIN: WARNING!!: This will load all the data and results for the analysis --------
# Uncomment if you want to use it. Some of the analysis step might require long
# computational times and you might want to use a computer with many cores/CPUs

# load("osd2014_16S_asv/data/osd2014_16S_asv_pina_tina.Rdata", verbose = TRUE)
# load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_asv_pina_tina.Rdata"), verbose = TRUE)

# END: WARNING!! ---------------------------------------------------------------



# BEGIN: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------

# Load necessary data -----------------------------------------------------
# If downloaded file at osd2014_16S_asv/data/ use:
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_with_phylo.Rdata")
load("osd2014_16S_asv/data/osd2014_sparcc_filtered.Rdata")

# If remote use
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_with_phylo.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_sparcc_filtered.Rdata"), verbose = TRUE)

# Load necessary data -----------------------------------------------------

# END: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------




PARAM <- list();
PARAM$cor.use <- "na.or.complete";
PARAM$p.adjust.method <- "hochberg";
PARAM$sample.steps <- c(0.01, 0.02, 0.05, 0.1, 0.15, seq(0.2, 0.9, by=0.1));
PARAM$use.cores <- 32;

#Set parameters
size.thresh <- 1;
pseudocount <- 10^-6;
nblocks <- 400;
use.cores <- 32;

########################
#Load packages for parallel processing
require("foreach");
require("bigmemory");
library("doMC", quietly=T);
#Register cluster
registerDoMC(cores=use.cores);
########################

########################
#Filter OTU table by removing all OTUs that are observed less then <size.thresh> times across all samples
#=> their correlations to all other OTUs will be (manually) set to 0 later on
#Add pseudocount to all observed counts (to avoid issues with zeroes in log-space)
cat(paste("SparCC preallocation steps =>", Sys.time()), sep="\n");
#o.t <- otu.table[1:1000, ] + pseudocount;
o.t <- t(as(object = otu_table(osd2014_dada2_phyloseq_beta), "matrix"))
ot <- t(as(object = otu_table(osd2014_dada2_phyloseq_beta), "matrix"))
o.t <- o.t + pseudocount;
otus <- rownames(o.t);
n.otu <- length(otus);
########################
#Preallocate blocks for parallel processing & Aitchinson's T matrix
#=> based on https://gist.github.com/bobthecat/5024079
size.split <- floor(n.otu / nblocks);
if (size.split < 1) {size.split <- 1}
my.split <- list(); length(my.split) <- nblocks;
my.split[1:(nblocks-1)] <- split(1:(size.split*(nblocks-1)), rep(1:(nblocks-1), each = size.split));
my.split[[nblocks]] <- (size.split*(nblocks-1)):n.otu;
dat.split <- mclapply(my.split, function(g) {o.t[g,]}, mc.cores=use.cores);
#Get combinations of splits
my.combs <- expand.grid(1:length(my.split), 1:length(my.split));
my.combs <- t(apply(my.combs, 1, sort));
my.combs <- unique(my.combs);
#Preallocate Aitchinson's T matrix as big.matrix ("shared" in memory, so accessible from w/in foreach loop)
mat.T <- big.matrix(nrow=n.otu, ncol=n.otu, dimnames=list(otus, otus), shared=T);
mat.T.desc <- describe(mat.T);
cat(paste("Done with preallocations =>", Sys.time()), sep="\n");
########################

########################
#Compute Aitchinson's T matrix
#=> iterate through each block combination, calculate matrix
#		between blocks and store them in the preallocated matrix on both
#		symmetric sides of the diagonal
cat(paste("Starting parallel T matrix calculation =>", Sys.time()), sep="\n");
results <- foreach(i = 1:nrow(my.combs)) %dopar% {
  #Get current combination
  curr.comb <- my.combs[i, ];
  #Get current data
  g.1 <- my.split[[curr.comb[1]]];
  g.2 <- my.split[[curr.comb[2]]];
  dat.1 <- dat.split[[curr.comb[1]]];
  dat.2 <- dat.split[[curr.comb[2]]];
  #Get current part of Aitchinson's matrix
  curr.T <- apply(dat.1, 1, function(x) {apply(dat.2, 1, function(y) {var(log(x/y))})});
  #Store
  curr.mat.T <- attach.big.matrix(mat.T.desc);
  curr.mat.T[g.2, g.1] <- curr.T;
  curr.mat.T[g.1, g.2] <- t(curr.T);
  #Return
  TRUE
}
cat(paste("Done with parallel T matrix calculation =>", Sys.time()), sep="\n");
########################

########################
#Compute component variations ("t_i")
cat(paste("Computing component variations =>", Sys.time()), sep="\n");
var.t <- colsum(mat.T);
cat(paste("Done with component variation calculations =>", Sys.time()), sep="\n");
#Estimate component variances ("omega_i") from t_i by solving a linear equation system
cat(paste("Estimating component variances from linear equation system =>", Sys.time()), sep="\n");
mat.a <- matrix(data=1, nrow=n.otu, ncol=n.otu); diag(mat.a) <- n.otu-1;
omega <- sqrt(solve(a=mat.a, b=var.t));
cat(paste("Done with component variation estimation =>", Sys.time()), sep="\n");
#Estimate pairwise correlations based on these values
cat(paste("Estimating correlations =>", Sys.time()), sep="\n");
global.sparcc <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %dopar% {(omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)}
rownames(global.sparcc) <- colnames(global.sparcc) <- otus;
cat(paste("Done with correlation estimation; returning data matrix =>", Sys.time()), sep="\n");
########################

########################
#Plot histogram
curr.data <- global.sparcc[upper.tri(global.sparcc, diag=F)];
ggplot(data.frame(data=sample(curr.data, 100000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + ggtitle("Distribution of Pairwise SparCC Correlations") + ylab("Density");
########################


############################
#Calculate derived OTU similarity matrix S
cat(paste("Calculating correlations of SparCC correlations =>", Sys.time()), sep="\n");
tmp.S <- cor.par(global.sparcc, method="pearson", use=PARAM$cor.use, use.cores=32);
cat(paste("Done =>", Sys.time()), sep="\n");
S.sparcc <- 0.5 * (tmp.S + 1);
#Plot histogram
curr.data <- S.sparcc[upper.tri(S.sparcc, diag=F)];
ggplot(data.frame(data=sample(curr.data, 10000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + ggtitle("Distribution of Pairwise SparCC Correlations") + ylab("Density");


#Calculate derived OTU similarity matrix S
cat(paste("Calculating correlations of SparCC correlations =>", Sys.time()), sep="\n");
tmp.S <- osd2014_sparcc_filtered
cat(paste("Done =>", Sys.time()), sep="\n");
S.sparcc.filt <- 0.5 * (tmp.S + 1);
#Plot histogram
curr.data <- S.sparcc.filt[upper.tri(S.sparcc.filt, diag=F)];
ggplot(data.frame(data=sample(curr.data, 10000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + ggtitle("Distribution of Pairwise SparCC Correlations") + ylab("Density");


############################
#Get cophenetic distance matrix for current tree
#=> based on the "cophenetic.phylo" function from the ape package
my.tree <- phy_tree(osd2014_dada2_phyloseq_beta);
my.tree$node.label <- c("Root", paste("Node", 2:my.tree$Nnode, sep = "_"));
global.cophenetic_tree <- fast.cophenetic.phylo(my.tree, use.cores=32);
#Reorder to match order in OTU table
global.cophenetic_tree <- global.cophenetic_tree[rownames(ot), rownames(ot)];
#Save
############################
#Calculate derived OTU similarity matrix S
cat(paste("Calculating correlations of cophenetic phylogenetic distances =>", Sys.time()), sep="\n");
tmp.S <- cor.par(global.cophenetic_tree, method="pearson", use=PARAM$cor.use, use.cores=32);
cat(paste("Done =>", Sys.time()), sep="\n");
S.phylo <- 0.5 * (tmp.S + 1);
#Plot histogram
curr.data <- S.phylo[upper.tri(S.phylo, diag=F)];
curr.plot <- ggplot(data.frame(data=sample(curr.data, 10000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + ggtitle("Distribution of Pairwise Phylogenetic Correlations") + ylab("Density");
#Save correlations & tidy up
rm(tmp.S, curr.data);
############################


############################
#Preallocate classical indices to calculate
get.cs <- list(); t <- 1;

############################
#Preallocate corrected indices to calculate
t <- start.t <- length(get.cs) + 1;
#Jaccard index, SparCC-corrected, unweighted, normalized
#=> unweighted TINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "TINA, unweighted";
get.cs[[t]]$call <- "jaccard.corr.uw.norm";
t <- t + 1;
#Jaccard index, SparCC-corrected, weighted, normalized
#=> weighted TINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "TINA, weighted";
get.cs[[t]]$call <- "jaccard.corr.w.norm";
t <- t + 1;
#Jaccard index, SparCC-corrected-filtered, unweighted, normalized
#=> unweighted TINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "TINA, unweighted filtered";
get.cs[[t]]$call <- "jaccard.corr.uw.norm";
t <- t + 1;
#Jaccard index, SparCC-corrected-filtered, weighted, normalized
#=> weighted TINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "TINA, weighted filtered";
get.cs[[t]]$call <- "jaccard.corr.w.norm";
t <- t + 1;
#Jaccard index, phylo-corrected, unweighted, normalized
#=> unweighted PINA
get.cs[[t]] <- list();
get.cs[[t]]$name <- "PINA, unweighted";
get.cs[[t]]$call <- "jaccard.corr.uw.norm";
t <- t + 1;
#Jaccard index, phylo-corrected, weighted, normalized
get.cs[[t]] <- list();
get.cs[[t]]$name <- "PINA, weighted";
get.cs[[t]]$call <- "jaccard.corr.w.norm";
t <- t + 1;



############################
#Iterate through (corrected) indices and calculate pairwise community similarities
for (t in seq(start.t, length(get.cs))) {
  print(paste(Sys.time(), "Starting", get.cs[[t]]$name));
  #Calculate all pairwise similarities
  if (get.cs[[t]]$name %in% c("TINA, unweighted", "TINA, weighted")) {
    curr.cs <- community.similarity.corr.par(ot, S=S.sparcc, distance=get.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
  } else if (get.cs[[t]]$name %in% c("PINA, unweighted", "PINA, weighted")) {
    curr.cs <- community.similarity.corr.par(ot, S=S.phylo, distance=get.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
  } else if (get.cs[[t]]$name %in% c("TINA, unweighted filtered", "TINA, weighted filtered")) {
    curr.cs <- community.similarity.corr.par(ot, S=S.sparcc.filt, distance=get.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
  }

  #Correct for rounding errors
  curr.cs[curr.cs < 0] <- 0;
  #Export histogram
  curr.data <- curr.cs[upper.tri(curr.cs)];
  curr.plot <- ggplot(data.frame(data=sample(curr.data, 6000)), aes(x=data)) + geom_density(alpha=0.2, fill="green") + xlim(0,1) + ggtitle(paste("Distribution of Pairwise", get.cs[[t]]$name, "Distances")) + ylab("Density");
  #ggsave(curr.plot, width=10, height=10, filename=paste(PARAM$folder.output, "hmp_samples.similarity_hist.", get.cs[[t]]$name, ".pdf", sep = ""), useDingbats=F);
  #Store current similarities
  get.cs[[t]]$cs <- curr.cs;
  get.cs[[t]]$plot <- curr.plot
  #Tidy up
  rm(curr.data, curr.cs);
  print(paste(Sys.time(), "Done with", get.cs[[t]]$name));
  gc();
}
############################

pina_tina_results <- get.cs

# BEGIN: Save objects ------------------------------------------------------------
# WARNING!!! You might not want to run this code --------------------------
save(pina_tina_results, file = "osd2014_16S_asv/data/osd2014_16S_asv_pina_tina_results.Rdata")
save.image("osd2014_16S_asv/data/osd2014_16S_asv_pina_tina.Rdata")
# END: Save objects ------------------------------------------------------------



