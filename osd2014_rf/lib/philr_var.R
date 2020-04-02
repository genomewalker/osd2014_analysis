# Here I use permutation of the leaf labels on the tree as the null model.

library(philr); packageVersion("philr")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(matrixStats); packageVersion('matrixStats')
library(parallel); packageVersion("parallel")
library(magrittr); packageVersion("magrittr")

set.seed(4)

#load('~/Downloads/HMP.RData')

# Helpful fxns ------------------------------------------------------------

calc.gm.coords <- function(site, cs){

  node_labels <-  BYSITE[[site]] %>%
    phy_tree() %>%
    .$node.label

  if (taxa_are_rows(BYSITE[[site]])){
    df <- t(otu_table(BYSITE[[site]]))
  }else{
    df <- (otu_table(BYSITE[[site]]))
  }

  cs <- intersect(node_labels, cs)

  df.clo <- df %>%
    #t() %>%
    as('matrix') %>%
    compositions::clo()

  sbp <- BYSITE[[site]] %>%
    phy_tree() %>%
    phylo2sbp()

  cs.l <- array_branch(cs)
  names(cs.l) <- cs

  map_df(cs.l, var.truncate.zeroes.df, df = df.clo, sbp = sbp)
}

var.truncate.zeroes <- function(coord, df, sbp, return.gm.only=FALSE,
                                return.balance.only=FALSE,
                                plot=FALSE, plot.type='pair', plot.log=TRUE){
  df <- as(df, 'matrix')
  up <- which(sbp[,coord]==1)
  down <- which(sbp[,coord]==-1)
  df.up <- df[,up, drop=F]
  df.down <- df[,down, drop=F]
  df.up[abs(df.up-0) < 1e-14] <- NA
  df.down[abs(df.down-0) < 1e-14] <- NA
  log.gm.up <- rowMeans(log(df.up), na.rm=TRUE)
  log.gm.down <- rowMeans(log(df.down), na.rm=TRUE)
  log.ratio <- log.gm.up - log.gm.down
  # Now calculate scaling constant
  n.up <- ncol(df.up)
  n.down <- ncol(df.down)
  sc <- sqrt(n.up*n.down/(n.up+n.down))
  # Now calculate variance
  v <- var(sc*log.ratio, na.rm=T)
  support <- sum(!is.na(log.ratio))

  if (return.gm.only==TRUE){
    return(list(gm.up=exp(log.gm.up),
                gm.down=exp(log.gm.down)))
  } else if (return.balance.only==TRUE){
    return(list(balance=sc*log.ratio))
  }

  if (plot==TRUE){
    subtitle <- paste('var:',v, 'support:',support, sep=' ')
    if (plot.type=='pair'){
      gm.up <- exp(log.gm.up)
      gm.down <- exp(log.gm.down)
      if (plot.log==TRUE) plot(x=gm.up, y=gm.down, main=coord, log='xy')
      else plot(x=gm.up, y=gm.down, main=coord)
    } else if(plot.type == 'hist'){
      hist(na.exclude(sc*log.ratio), main=coord)
    }
    mtext(subtitle)
    return(NULL)
  }
  return(c(v, support))
}

var.truncate.zeroes.df <- function(coord, df, sbp, return.gm.only=FALSE,
                                return.balance.only=FALSE,
                                plot=FALSE, plot.type='pair', plot.log=TRUE){
  df <- as(df, 'matrix')
  up <- which(sbp[,coord]==1)
  down <- which(sbp[,coord]==-1)
  df.up <- df[,up, drop=F]
  df.down <- df[,down, drop=F]
  df.up[abs(df.up-0) < 1e-14] <- NA
  df.down[abs(df.down-0) < 1e-14] <- NA
  log.gm.up <- rowMeans(log(df.up), na.rm=TRUE)
  log.gm.down <- rowMeans(log(df.down), na.rm=TRUE)
  log.ratio <- log.gm.up - log.gm.down
  # Now calculate scaling constant
  n.up <- ncol(df.up)
  n.down <- ncol(df.down)
  sc <- sqrt(n.up*n.down/(n.up+n.down))
  # Now calculate variance
  v <- var(sc*log.ratio, na.rm=T)
  support <- sum(!is.na(log.ratio))

  if (return.gm.only==TRUE){
    return(list(gm.up=exp(log.gm.up),
                gm.down=exp(log.gm.down)))
  } else if (return.balance.only==TRUE){
    return(list(balance=sc*log.ratio))
  }

  if (plot==TRUE){
    subtitle <- paste('var:',v, 'support:',support, sep=' ')
    if (plot.type=='pair'){
      gm.up <- exp(log.gm.up)
      gm.down <- exp(log.gm.down)
      if (plot.log==TRUE) plot(x=gm.up, y=gm.down, main=coord, log='xy')
      else plot(x=gm.up, y=gm.down, main=coord)
    } else if(plot.type == 'hist'){
      hist(na.exclude(sc*log.ratio), main=coord)
    }
    mtext(subtitle)
    return(NULL)
  }
  return(tibble(coord = coord, var = v, support = support, gm_up = exp(log.gm.up), gm_down = exp(log.gm.down), label = names(exp(log.gm.up))))
}

# Calls var.truncate.zeroes
calc.var <- function(sbp, bmd, df, n.support=10, return.var=FALSE){
  # Var Truncate Zeroes
  df <- as(df, 'matrix')
  df <- compositions::clo(df)
  var.tz <- sapply(colnames(sbp), var.truncate.zeroes, df, sbp)
  var.tz <- t(var.tz)
  colnames(var.tz) <- c('var.tz', 'support')
  var.tz <- add_rownames(as.data.frame(var.tz), var = 'coord')

  # Add distance to tips
  var.tz$mean.dist.to.tips <- bmd[var.tz$coord]
  if (return.var==TRUE)return(var.tz)

  # Set threshold below which can't calculate variance and fit model
  fit <- var.tz %>%
    filter(var.tz !=0, mean.dist.to.tips !=0) %>%
    filter(support >= n.support) %>%
    lm(log(var.tz) ~log(mean.dist.to.tips), data=.) %>%
    broom::tidy()
  fit$term  <- factor(c('intercept', 'log_mdtt'))

  return(fit)
}

null.model <- function(i, sbp, bmd, df, n.support){
  i <- NULL # Just for mclapply
  sbp.perm <- sbp
  colnames(sbp.perm) <- sample(colnames(sbp.perm),
                               length(colnames(sbp.perm)),
                               replace=FALSE)
  return(calc.var(sbp.perm, bmd, df, n.support, return.var=FALSE))
}

# High level
# n - number of permutations
# nmr - precomputed null model results
# n.support - number of non-zero/missing points/samples needed to calculate var
run.var.analysis <- function(phyloseq.filt, n, ncores, n.support, nmr=NULL){
  if (taxa_are_rows(phyloseq.filt)){
    df <- t(otu_table(phyloseq.filt))
  }else{
    df <- (otu_table(phyloseq.filt))
  }
  tr <- phy_tree(phyloseq.filt)
  bmd <- mean_dist_to_tips(tr)
  sbp <- phylo2sbp(tr)

  # Now calculate var
  var.tz <- calc.var(sbp, bmd, df, return.var=TRUE, n.support = 10)

  # Fit linear model
  lm <- var.tz %>%
    filter(var.tz !=0, mean.dist.to.tips !=0) %>%
    filter(support >= n.support) %>%
    lm(log(var.tz) ~log(mean.dist.to.tips), data=.) %>%
    summary()
  print(lm)

  if (is.null(nmr)){
    # Now create null model/permute
    null.model.results <- pbmcapply::pbmclapply(1:n, null.model, sbp, bmd, df, n.support, mc.cores=ncores)
    nmr <- bind_rows(null.model.results, .id='id')
    nmr %<>%
      select(id, term, estimate) %>%
      spread(term,estimate)
  }

  # Display Null distribution
  p.null <- ggplot(nmr, aes(x=log_mdtt)) +
    geom_density() +
    geom_vline(aes(xintercept=lm$coefficients[2,'Estimate']), color='red') +
    ggtitle(paste('Null Distribution for Slope',
                  ggtitle(deparse(substitute(phyloseq.filt))),
                  sep=': '))

  # Plot var vs. depth with null models plotted as well
  p.pvd <- var.tz %>%
    filter(var.tz !=0, mean.dist.to.tips !=0) %>%
    filter(support >= n.support) %>%
    ggplot(aes(x=log(mean.dist.to.tips), y=log(var.tz), color=support), alpha=0.5) +
    geom_point() +
    geom_abline(aes(intercept=intercept, slope=log_mdtt), data=nmr, color='darkgreen', alpha=0.01) +
    geom_smooth(method='lm', se=FALSE)+
    ggtitle(deparse(substitute(phyloseq.filt)))

  # Is there a relationship between support and var?
  p.mp <- ggplot(var.tz, aes(x=support, y=log(var.tz), color=log(mean.dist.to.tips))) +
    geom_point() +
    ggtitle(deparse(substitute(phyloseq.filt)))

  result <- list(var.tz=var.tz,
                 lm = lm,
                 nmr=nmr,
                 p.null=p.null,
                 p.pvd=p.pvd,
                 #p.tree=p.tree,
                 p.mp=p.mp)
  return(result)
}


# From https://github.com/DomBennett/MoreTreeTools/blob/master/R/g --------
getChildren <- function(tree, node, display = FALSE) {
  if (!is.numeric (node)) {
    stop("Node is not numeric!")
  }
  if (node > tree$Nnode + length (tree$tip.label)) {
    stop("Node is greater than the number of nodes in tree!")
  }
  if (node <= length (tree$tip.label)) {
    term.nodes <- node
  } else {
    term.nodes <- vector ()
    temp.nodes <- node
    while (length (temp.nodes) > 0) {
      connecting.nodes <- tree$edge[tree$edge[,1] %in% temp.nodes, 2]
      term.nodes <- c(term.nodes, connecting.nodes[connecting.nodes <=
                                                     length (tree$tip.label)])
      temp.nodes <- connecting.nodes[connecting.nodes > length(tree$tip.label)]
    }
  }
  children <- tree$tip.label[term.nodes]
  if (display) {
    tip.cols <- ifelse(tree$tip.label %in% children, "black", "grey")
    plot.phylo(tree, tip.color = tip.cols, show.tip.label = TRUE)
    nodelabels("node", node)
  }
  return(children)
}

tree_subset_mod <- function (tree, node, levels_back = 5, group_node = TRUE)
{
  if (!is.numeric(levels_back)) {
    levels_back <- as.numeric(levels_back)
    if (is.na(levels_back))
      stop("'levels_back' must be of class numeric")
  }
  tree_df <- tidytree::as_tibble(tree)
  selected_node <- node
  is_tip <- tree_df %>% dplyr::mutate(isTip = !.data$node %in%
                                        .data$parent) %>% dplyr::filter(.data$node == selected_node |
                                                                          .data$label == selected_node) %>% dplyr::pull(.data$isTip)
  if (is_tip & levels_back == 0) {
    stop("The selected node (", selected_node, ") is a tip. 'levels_back' must be > 0",
         call. = FALSE)
  }
  if (is_tip) {
    group_labels <- tree_df %>% dplyr::filter(.data$node ==
                                                selected_node | .data$label == selected_node) %>%
      dplyr::pull(.data$label)
  }
  else {
    group_labels <- tree_df %>% tidytree::offspring(selected_node) %>%
      dplyr::filter(!.data$node %in% .data$parent) %>%
      dplyr::pull(.data$label)
  }
  if (levels_back == 0) {
    subset_labels <- tidytree::offspring(tree_df, selected_node) %>%
      dplyr::filter(!.data$node %in% .data$parent) %>%
      dplyr::pull(.data$label)
  }
  else {
    subset_labels <- tidytree::ancestor(tree_df, selected_node) %>%
      tail(levels_back) %>% head(1) %>% dplyr::pull(.data$node) %>%
      tidytree::offspring(tree_df, .) %>% dplyr::filter(!.data$node %in%
                                                          .data$parent) %>% dplyr::pull(.data$label)
  }
  subset_nodes <- which(tree$tip.label %in% subset_labels)
  subtree <- drop.tip(tree, tree$tip.label[-subset_nodes],
                      rooted = TRUE)
  if (group_node)
    subtree <- tidytree::groupOTU(subtree, .node = group_labels)
  return(subtree)
}


simplify.phylum <- function(tax, thresh=50){
  # Define color presets
  phylum.colors <- c('Actinobacteria'='#BFBBFF', 'Bacteroidetes'='#CB4D42',
                     'Cyanobacteria'='#66CDAA', 'Fusobacteria'='#EEFF9A',
                     'Firmicutes'='#FFB756','Proteobacteria'='#00569D',
                     'Spirochaetes'='#9B00A6', 'Synergistetes'='#CDAF95',
                     'Tenericutes'='#FFEA95', 'Verrucomicrobia'='#68688E',
                     'Other'='#FFAAAA')

  # Retain only phyla that have more than thresh representatives for the given tax
  # table - lump everythin else in with "Other"
  phylum.to.keep <- c(names(table(tax)[table(tax) > thresh]), 'Other')
  phylum.colors <- phylum.colors[phylum.to.keep]

  # Simplify Tax vector by only keeping certain ones - combine the rest
  # to "Others"
  if(!is.factor(tax))tax <- factor(tax) # Make sure we are dealing with a factor
  lt <- levels(tax)
  lt[!(lt %in% names(phylum.colors))] <- 'Other'
  levels(tax) <- lt

  # Reorder levels so Other is at the end
  lt <- c(lt[lt!='Other'], 'Other')
  tax <- factor(tax, levels=lt)

  list(tax=tax, phylum.colors=phylum.colors)
}

plot.pvd <- function(var.tz, labels=NULL, n.support=40, blw.byrank=NULL,
                     plot.rug=F, plot.rank.lines=T, rank.species.only=F,
                     point.size=1){
  df <- var.tz %>%
    filter(var.tz !=0, mean.dist.to.tips !=0) %>%
    filter(support >= n.support)

  # Plot basic scatterplot
  p <- ggplot(df, aes(x=mean.dist.to.tips, y=var.tz)) +
    geom_point(size=point.size, shape = 21, color = "black", fill = "grey", alpha = 0.8)

  # Add layer of taxonomic rank info if needed
  if (!is.null(blw.byrank)){
    if (plot.rug){
      p <- p + geom_rug(data=(blw.byrank %>% filter(rank %in% c('species'))),
                        aes(x=mean.dist.to.tips, y=NULL, color=rank), lwd=1, show.legend = F)
    }
    if (plot.rank.lines){
      if (rank.species.only){
        blw.byrank %<>% filter(rank %in% c('species'))
      }
      blw.byrank  %<>%  group_by(rank)  %>% summarise(median=median(mean.dist.to.tips))
      p <- p + geom_vline(data=blw.byrank, aes(xintercept=median), linetype='dashed') #+
      #geom_text(data=blw.byrank, aes(x=median, y=0.25, label=substr(rank, 1, 1)))
    }
  }

  # Now add loess regression and linear regression then labels
  p <- p + geom_smooth(method='loess', color='darkgreen', lwd = 0.5) +
    geom_smooth(method='lm', se=FALSE,  lwd = 0.5)
  if (!is.null(labels)){
    data <- var.tz %>% filter(coord %in% labels)
    #p <- p + geom_label_repel(aes(label=coord), data=data)
  }
  # Some final stylizing
  p +
    theme_bw()
}

figure.pvd <- function(site, to.plot, plot.rug=F, plot.rank.lines=T, rank.species.only=F,
                       point.size=1){
  tr <- tax2tree(phy_tree(BYSITE[[site]]), BYSITE[[site]])
  blw.byrank <- extract.rank.data(tr)
  p <- plot.pvd(bysite[[site]]$var.tz, labels=to.plot, blw.byrank=blw.byrank,
                plot.rug=plot.rug, plot.rank.lines=plot.rank.lines,
                rank.species.only=rank.species.only,
                point.size=point.size, n.support = 20)
  p
}


figure.pvd_18S <- function(site, to.plot, plot.rug=F, plot.rank.lines=T, rank.species.only=F,
                       point.size=1){
  tr <- tax2tree_18S(phy_tree(BYSITE[[site]]), BYSITE[[site]])
  blw.byrank <- extract.rank.data_18S(tr)
  p <- plot.pvd(bysite[[site]]$var.tz, labels=to.plot, blw.byrank=blw.byrank,
                plot.rug=plot.rug, plot.rank.lines=plot.rank.lines,
                rank.species.only=rank.species.only,
                point.size=point.size, n.support = 20)
  p
}





require(phyloseq)
require(magrittr)
require(ape)
require(tidyr)
require(dplyr)

# Given a tree and a tax table (corresponding to that tree) runs
# tax2tree (the python scripts) and outputs a tree decorated with the
# taxonomy - will strip existing node names on the tree.
# tax should be a filename - should not have headers
# Cleanup removes temporary tax2tree output files.
# Note path to nlevel is hardcoded currently TODO:
tax2tree <- function(tr, physeq, cleanup=TRUE){
  tr <- phy_tree(physeq)
  df <- as(tax_table(physeq), "matrix") %>%
    as_tibble(rownames = "asv") %>%
    mutate(Kingdom = paste0("k__", Kingdom),
           Phylum = paste0("p__", Phylum),
           Class = paste0("c__", Class),
           Order = paste0("o__", Order),
           Family = paste0("f__", Family),
           Genus = paste0("g__", Genus),
           Species = paste0("s__", Species)) %>%
    unite(Kingdom, Phylum, Class, Order, Family, Genus, Species, col = "tax", sep = "; ", remove = TRUE) %>%
    mutate(tax = gsub("__NA", "__", tax))
  # Write cleaned to file
  write.table(df, file='/tmp/tmp.tax.cleaned.txt', sep='\t',
              quote=F, col.names=F, row.names=F)

  # Write tree with labels removed
  tr$node.label <- NULL
  write.tree(tr, file='/tmp/tmp.tree')

  # Run tax2tree
  system("pyenv local miniconda2-latest; /Users/ufo/.pyenv/shims/nlevel -t /tmp/tmp.tree -m /tmp/tmp.tax.cleaned.txt -o /tmp/tax2tree_output", wait = T)
  while(!file.exists('/tmp/tax2tree_output')){Sys.sleep(1)}

  # Note Documentation of problem
  # https://groups.google.com/forum/#!msg/qiime-forum/v-TfjP20uws/IGkeWpUO6rwJ
  # This is the reason I am using phyloseq's tree reader for greengenes style trees
  # That apparently tax2tree also outputs in.
  tr <- phyloseq::read_tree_greengenes('/tmp/tax2tree_output')

  # remove weirdness
  tr$node.label <- gsub('><-><','.',tr$node.label)

  # Remove extra quotes around names:
  # http://joey711.github.io/phyloseq-demo/HMP_import_example.html
  tr$node.label <- gsub("'","", tr$node.label)
  tr$tip.label <- gsub("'","", tr$tip.label)

  # Remove tax2tree output, and tmp.tax, tmp.tree
  system('rm /tmp/tmp.tax.cleaned.txt')
  system('rm /tmp/tmp.tree')
  if (cleanup)system('rm /tmp/tax2tree_output*')

  # return annotated tree
  tr
}

# # called from tax2tree - input is data.frame to be 'cleaned' for
# # tax2table - e.g., make sure everything has the same number of ranks
# # (even if they are empty) and ensure that stuff uses the "Unclassified"
# # keyword.
# clean.tax <- function(df) {
#   df %<>% separate(V2, paste('c', 1:7,sep=''), sep=';')
#
#   # Fix the NAs
#   df[['c7']][is.na(df[['c7']])] <- ' s__'
#   df[['c6']][is.na(df[['c6']])] <- ' g__'
#   df[['c5']][is.na(df[['c5']])] <- ' f__'
#   df[['c4']][is.na(df[['c4']])] <- ' o__'
#   df[['c3']][is.na(df[['c3']])] <- ' c__'
#   df[['c2']][is.na(df[['c2']])] <- ' p__'
#   df[['c1']][is.na(df[['c1']])] <- ' k__'
#
#   # Fix "Unassigned"
#   if (any(df$c1=='Unassigned')){
#     df[df$c1=='Unassigned',c(3,8)] <- rep(NA, 6)
#     df[df$c1=='Unassigned','c1'] <- 'Unclassified'
#   }
#
#   # Strip Whitespace
#   # returns string w/o leading or trailing whitespace
#   trim <- function (x) gsub("^\\s+|\\s+$", "", x)
#   df %<>% sapply(trim) %>% as.data.frame(stringsAsFactors=F)
#
#   # Collapse taxonomy to a single column
#   df %<>% unite(c, c1, c2, c3, c4, c5, c6, c7, sep='; ')
#   df
# }
#
#
# extract.rank.data <- function(tr) {
#   # calculate mean.dist.to.tips for each node
#   blw <- mean_dist_to_tips(tr)
#
#   # If a given node has multiple ranks assigned take the lowest.
#   blw.lowest.names <- strsplit(names(blw), '.', fixed=TRUE) %>%
#     map(~(ifelse(length(.x)>0, .x[length(.x)], ""))) %>%
#     as_vector()
#   names(blw) <- blw.lowest.names
#
#   # Depth of species
#   blw.byrank <- list()
#   blw.byrank[['species']] <- blw[grep('s__[A-Za-z]+', names(blw))]
#   blw.byrank[['genus']] <- blw[grep('g__[A-Za-z]+', names(blw))]
#   blw.byrank[['family']] <- blw[grep('f__[A-Za-z]+', names(blw))]
#   blw.byrank[['order']] <- blw[grep('o__[A-Za-z]+', names(blw))]
#   blw.byrank[['phylum']]<- blw[grep('p__[A-Za-z]+', names(blw))]
#   blw.byrank[['kingdom']]<- blw[grep('k__[A-Za-z]+', names(blw))]
#
#   blw.byrank %<>%
#     map(~data.frame(mean.dist.to.tips=.x)) %>%
#     bind_rows(.id='rank')
#
#   # Reorder factors
#   blw.byrank$rank <- factor(blw.byrank$rank,
#                             levels=c('species', 'genus','family', 'order',
#                                      'phylum', 'kingdom'))
#   blw.byrank
# }


# By Site - Analyze -------------------------------------------------------
# BYSITE <- list()
# for (site in levels(get_variable(HMP, 'groupedsites'))){
#   SITE <- prune_samples(sample_data(HMP)$groupedsites == site, HMP)
#   SITE <- filter_taxa(SITE, function(x) sum(x > 1) > (0.2*length(x)) , TRUE)
#   SITE <- prune_samples(colSums(otu_table(SITE)) > 50, SITE)
#   BYSITE[[site]] <- SITE
# }
#
# bysite = list()
# for (site in levels(get_variable(HMP, 'groupedsites'))){
#   site <- "GU:Posterior_fornix"
#   bysite[[site]] <- run.var.analysis(BYSITE[[site]], n=20000, n.support = 40, ncores=3)
# }
# save(bysite, file='~/Downloads/bysite.20000.40.RData')
