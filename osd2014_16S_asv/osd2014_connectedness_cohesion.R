# Online script to generate cohesion metrics for a set of samples
# CMH 06Dec17; cherren@wisc.edu

# User instructions: read in a sample table (in absolute or relative abundance) as object "b".
# If using a custom correlation matrix, read in that matrix at the designated line.
# Run the entire script, and the 4 vectors (2 of connectedness and 2 of cohesion) are generated for each sample at the end.
# Parameters that can be adjusted include pers.cutoff (persistence cutoff for retaining taxa in analysis),
# iter (number of iterations for the null model), tax.shuffle (whether to use taxon shuffle or row shuffle randomization),
# and use.custom.cors (whether to use a pre-determined correlation matrix)

library(tidyverse)
library(ggridges)
library(phyloseq)

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

# load("osd2014_16S_asv/data/osd2014_connectedness_cohesion.Rdata", verbose = TRUE)
# load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_connectedness_cohesion.Rdata"), verbose = TRUE)

# END: WARNING!! ---------------------------------------------------------------


# BEGIN: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------

# Load necessary data -----------------------------------------------------
# Use if you have the postgres DB in place
my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)
st_100_order_terrestrial <- tbl(my_db, "osd2014_st_order_coastal") %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label)
osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)
osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf)
osd2014_meow_regions <- tbl(my_db, "osd2014_meow_regions") %>%
  collect(n = Inf)

# If downloaded file at osd2014_16S_asv/data/ use:

load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_with_phylo.Rdata", verbose = TRUE)
load("osd2014_16S_asv/data/osd2014_sparcc_filtered.Rdata", verbose = TRUE)

# Basic contextual data
load("osd2014_16S_asv/data/osd2014_basic_cdata.Rdata", verbose = TRUE)

# If remote use
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects_with_phylo.Rdata"), verbose = TRUE)
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_sparcc_filtered.Rdata"), verbose = TRUE)

# Basic contextual data
load(url("http://osd2014.metagenomics.eu/osd2014_16S_asv/data/osd2014_basic_cdata.Rdata"), verbose = TRUE)
# Load necessary data -----------------------------------------------------

# END: SKIP THIS IF YOU ALREADY LOADED ALL RESULTS AND DATA --------------------





####################create necessary functions######################

#find the number of zeroes in a vector
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}

#create function that averages only negative values in a vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}

#create function that averages only positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}

###################################################################
###################################################################
### Workflow options ####
###################################################################
###################################################################

## Choose a persistence cutoff (min. fraction of taxon presence) for retaining taxa in the analysis
pers.cutoff <- 0.10
## Decide the number of iterations to run for each taxon. (>= 200 is recommended)
# Larger values of iter mean the script takes longer to run
iter <- 200
## Decide whether to use taxon/column shuffle (tax.shuffle = T) or row shuffle algorithm (tax.shuffle = F)
tax.shuffle <- T
## Option to input your own correlation table
# Note that your correlation table MUST have the same number of taxa as the abundance table.
# There should be no empty (all zero) taxon vectors in the abundance table.
# Even if you input your own correlation table, the persistence cutoff will be applied
use.custom.cors <- T

###################################################################
###################################################################

# Read in dataset
## Data should be in a matrix where each row is a sample.

b <- as(otu_table(osd2014_dada2_phyloseq_beta), "matrix")
custom.cor.mat <- as.matrix(osd2014_sparcc_filtered)

# Read in custom correlation matrix, if desired. Must set "use.custom.cors" to TRUE
if(use.custom.cors == T) {
  #custom.cor.mat <- read.csv("your_path_here.csv", header = T, row.names = 1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  #Check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2] == dim(custom.cor.mat)[2])
}


# Suggested steps to re-format data. At the end of these steps, the data should be in a matrix "c" where there are no
# empty samples or blank taxon columns.
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]

# Optionally re-order dataset to be in chronological order. Change date format for your data.
#c <- c[order(as.Date(rownames(c), format = "%m/%d/%Y")), ]

# Save total number of individuals in each sample in the original matrix. This will be 1 if data are in relative abundance,
# but not if matrix c is count data
rowsums.orig <- rowSums(c)

# Based on persistence cutoff, define a cutoff for the number of zeroes allowed in a taxon's distribution
zero.cutoff <- ceiling(pers.cutoff * dim(c)[1])

# Remove taxa that are below the persistence cutoff
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
# Remove any samples that no longer have any individuals, due to removing taxa
d <- d[rowSums(d) > 0, ]

#If using custom correlation matrix, need to remove rows/columns corresponding to the taxa below persistence cutoff
if(use.custom.cors == T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c, 2, zero) < (dim(c)[1]-zero.cutoff), apply(c, 2, zero) < (dim(c)[1]-zero.cutoff)]
}

# Create relative abundance matrix.
rel.d <- d / rowsums.orig
# Optionally, check to see what proportion of the community is retained after cutting out taxa
hist(rowSums(rel.d))

# Create observed correlation matrix
cor.mat.true <- cor(rel.d)

# Create vector to hold median otu-otu correlations for initial otu
med.tax.cors <- vector()

# Run this loop for the null model to get expected pairwise correlations
# Bypass null model if the option to input custom correlation matrix is TRUE
if(use.custom.cors == F) {
  if(tax.shuffle) {
    for(which.taxon in 1:dim(rel.d)[2]){

      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()

      for(i in 1:iter){
        #Create empty matrix of same dimension as rel.d
        perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)

        #For each otu
        for(j in 1:dim(rel.d)[2]){
          # Replace the original taxon vector with a permuted taxon vector
          perm.rel.d[, j ] <- sample(rel.d[ ,j ])
        }

        # Do not randomize focal column
        perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]

        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)

        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])

      }
      # Save the median correlations between the focal taxon and all other taxa
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))

      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  } else {
    for(which.taxon in 1:dim(rel.d)[2]){

      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()

      for(i in 1:iter){
        #Create duplicate matrix to shuffle abundances
        perm.rel.d <- rel.d

        #For each taxon
        for(j in 1:dim(rel.d)[1]){
          which.replace <- which(rel.d[j, ] > 0 )
          # if the focal taxon is greater than zero, take it out of the replacement vector, so the focal abundance stays the same
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]

          #Replace the original taxon vector with a vector where the values greater than 0 have been randomly permuted
          perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal])
        }

        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)

        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])

      }
      # Save the median correlations between the focal taxon and all other taxa
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))

      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  }
}

# Save observed minus expected correlations. Use custom correlations if use.custom.cors = TRUE
if(use.custom.cors == T) {
  obs.exp.cors.mat <- custom.cor.mat.sub
} else {
  obs.exp.cors.mat <- cor.mat.true - med.tax.cors
}

diag(obs.exp.cors.mat) <- 0

####
#### Produce desired vectors of connectedness and cohesion

# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

####
#### Combine vectors into one list and print
output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)

names(output) <- c("Negative Connectedness", "Positive Connectedness", "Negative Cohesion", "Positive Cohesion")

cohesion_df <- bind_rows(tibble(label = rownames(output$`Negative Cohesion`), cohesion = output$`Negative Cohesion`[,1], class = "negative"),
                         tibble(label = rownames(output$`Positive Cohesion`), cohesion = output$`Positive Cohesion`[,1], class = "positive")) %>%
  mutate(label = fct_relevel(label, st_100_order_terrestrial$label)) %>%
  left_join(osd2014_cdata) %>%
  filter(meow_province %in% osd2014_meow_regions$meow_province)

ggplot(cohesion_df, aes(y = meow_province, x = cohesion, fill = class)) +
  geom_density_ridges(aes(point_fill = class),
                      jittered_points = TRUE, point_color = "grey20", point_shape = 21, point_size = 0.8, scale = 0.7, alpha = .7,
                      rel_min_height = 0.0 ,panel_scaling = T, size = 0.2, color = "black")+
  theme_light()


cohesion_df <- tibble(label = rownames(output$`Negative Cohesion`), cohesion_negative = output$`Negative Cohesion`[,1]) %>%
  left_join(tibble(label = rownames(output$`Positive Cohesion`), cohesion_positive = output$`Positive Cohesion`[,1]))

connectedness_df <- tibble(asv = names(connectedness.neg), connectedness_negative = connectedness.neg) %>%
  left_join(tibble(asv = names(connectedness.pos), connectedness_positive = connectedness.pos))

# BEGIN: Save objects ------------------------------------------------------------
# WARNING!!! You might not want to run this code --------------------------
save.image("osd2014_16S_asv/data/osd2014_connectedness_cohesion.Rdata")
# END: Save objects ------------------------------------------------------------

