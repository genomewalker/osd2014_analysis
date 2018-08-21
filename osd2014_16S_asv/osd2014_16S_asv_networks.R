library(igraph)
library(phyloseq)
#library(SpiecEasi)
library(tidyverse)
library(tidygraph)
library(ggraph)
source("osd2014_16S_asv/lib/graph_lib.R")
load("osd2014_16S_asv/data/osd2014_16S_asv_networks.Rdata", verbose = TRUE)
load("osd2014_16S_asv/data/osd2014_16S_niche_breadth.Rdata", verbose = TRUE)
load("osd2014_16S_asv/data/osd2014_16S_asv_physeq_filt_objects.Rdata", verbose = TRUE)

rm(osd2014_16S_asv_se_gl_minus2)
rm(osd2014_16S_asv_se_mb_minus3)

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")
osd2014_silva_dada2_names <- tbl(my_db, "osd2014_silva_dada2") %>%
  collect(n = Inf) %>%
  select(asv, asv_name)

st_100_order_terrestrial <- tbl(my_db, "osd2014_st_order_terrestrial") %>%
  collect(n = Inf)

osd2014_amp_mg_intersect <- tbl(my_db, "osd2014_amp_mg_intersect_2018") %>%
  collect(n = Inf)

osd2014_sample_cohesion <- tbl(my_db, "osd2014_sample_cohesion") %>%
  collect(n = Inf)

osd2014_cdata <- tbl(my_db, "osd2014_cdata") %>%
  collect(n = Inf) %>%
  filter(label %in% osd2014_amp_mg_intersect$label) %>%
  left_join(osd2014_sample_cohesion) %>%
  arrange(match(label, st_100_order_terrestrial$label))

osd2014_sparcc_g <- osd2014_sparcc %>%
  rename(weight = correlation) %>%
  mutate(weight_orig = weight, weight = abs(weight)) %>%
  graph_from_data_frame(directed = FALSE)


# library(Matrix)
# secor <- cov2cor(forceSymmetric(getOptCov(osd2014_16S_asv_se_gl_minus2), ifelse(sum(Matrix::tril(getOptCov(osd2014_16S_asv_se_gl_minus2)))>sum(Matrix::triu(getOptCov(osd2014_16S_asv_se_gl_minus2))), 'L', 'U')))
# refit <- forceSymmetric(getRefit(osd2014_16S_asv_se_gl_minus2), ifelse(sum(Matrix::tril(getRefit(osd2014_16S_asv_se_gl_minus2)))>sum(Matrix::triu(getRefit(osd2014_16S_asv_se_gl_minus2))), 'L', 'U'))
# osd2014_asv_g_gl <- adj2igraph(as.matrix(secor*refit),  vertex.attr=list(name=taxa_names(osd2014_16s_otuXsample_physeq_filt_se)))
# E(osd2014_asv_g_gl)$weight %>% hist
#
#
# osd2014_asv_g_mb <- symBeta(getOptBeta(osd2014_16S_asv_se_mb_minus3), mode='maxabs')
# osd2014_asv_g_mb <- adj2igraph(osd2014_asv_g_mb, vertex.attr=list(name=taxa_names(osd2014_16s_otuXsample_physeq_filt_se)))
# E(osd2014_asv_g_mb)$weight %>% hist
#
binary_search_filter <- function(g, low, high) {
  x <- g

  if ( high < low ) {

    return(NULL)

  }else {
    mid <- floor((low + high) / 2)
    x<-igraph::delete.edges(g, which(E(g)$weight <= weights[mid]))
    cat("low:", low, " mid:", mid, " high:", high, " mid_weight:", weights[mid], " components:", count_components(x), " edges:", ecount(x), "\n")

    if ((mid - low) == 0){
      x<-igraph::delete.edges(g, which(E(g)$weight <= weights[mid-1]))
      cat("Final: low:", low - 1, " mid:", mid - 1, " high:", high - 1, " mid_weight:", weights[mid - 1], " components:",count_components(x), " edges:", ecount(x), "\n")
      return(list(weight=mid - 1, graph=x))
      exit
    }

    if (((high - low) == 0)){
      x<-igraph::delete.edges(g, which(E(g)$weight <= weights[mid+1]))
      cat("Final: low:", low, " mid:", mid, " high:", high, " mid_weight:", weights[mid], " components:",count_components(x), " edges:", ecount(x), "\n")
      return(list(weight=mid , graph=x))
      exit
    }


    if (!(is.connected(x))){
      binary_search_filter(g, low, mid - 1)
    }else if (is.connected(x)){
      binary_search_filter(g, mid + 1, high)
    }
  }
}

G <- osd2014_sparcc_g

weights <- unique(sort(E(G)$weight, decreasing = FALSE))
g3 <- binary_search_filter(g = G, low = 1, high = length(weights))$graph


E(g3)$pvalue %>% hist
E(g3)$weight_orig %>% hist

# Invert p-values
E(g3)$pvalue_inv <- 1 - E(g3)$pvalue
E(g3)$pvalue_inv %>% hist
# Assign inverted p-values to weight
E(g3)$weight <- E(g3)$pvalue_inv

weights <- unique(sort(E(g3)$weight, decreasing = FALSE))
g4 <- binary_search_filter(g = g3, low = 1, high = length(weights))$graph


E(g4)$weight <- abs(E(g4)$weight_orig)
E(g4)$pvalue %>% skimr::skim()
E(g4)$weight_orig %>% skimr::skim()






#
#
# a_names <- paste("asv", 1:vcount(g3), sep = "_")
#
# g_taxonomy <- tax_table(osd2014_16s_otuXsample_physeq_filt_prev_beta) %>%
#   as.matrix() %>%
#   base::as.data.frame( stringsAsFactors = FALSE) %>%
#   rownames_to_column(var = "asv") %>%
#   tbl_df %>%
#   mutate(id = paste("asv", as.character(row_number()), sep = "_")) %>% select(id, asv)
#
#
#
# # Transform grapj
#
# g <- adj2igraph(symBeta(getOptBeta(osd2014_se_mb_beta), mode='ave'))
#
# hist(summary(symBeta(getOptBeta(osd2014_se_mb_beta), mode='ave')))
#
# g<-adj2igraph(se_asv_gl_1e4$opt.cov * se_asv_gl_1e4$refit)
# a_names <- paste("asv", 1:ntaxa(osd2014_16s_otuXsample_physeq_filt_prev_beta), sep = "_")
# #r_names <- tax_table(osd2014_resfam_physeq_se) %>% as.data.frame() %>% .$resfam_id %>% as.character()
# #g_names <- c(a_names, r_names)
#
#
# g_taxonomy <- tax_table(osd2014_16s_otuXsample_physeq_filt_prev_beta) %>%
#   as.matrix() %>%
#   base::as.data.frame( stringsAsFactors = FALSE) %>%
#   rownames_to_column(var = "asv") %>%
#   tbl_df %>%
#   mutate(id = paste("asv", as.character(row_number()), sep = "_")) %>% select(id, asv)
#


# g3 <- osd2014_asv_g_gl
g <- as_tbl_graph(g4) %>%
  activate(nodes) %>%
  mutate(asv = name) %>%
  left_join(osd2014_silva_dada2_names) %>%
  mutate(name = asv_name) %>%
  select(-asv_name) %>%
  #filter(!node_is_isolated()) %>%
  activate(edges) %>%
  select(-pvalue_inv) %>%
  #mutate(weight_orig = weight, weight = abs(weight)) %>%
  #filter(pvalue < 0.01) %>%
  #mutate(sign = ifelse(weight_orig < 0, "NEG", "POS")) %>%
  activate(nodes) %>%
  filter(!node_is_isolated())


#gx<-g
# Louvain
community_bin <- file.path("~/Desktop/BiG-SCAPE/bin/community ")
convert_bin <- file.path("~/Desktop/BiG-SCAPE/bin/convert")
hierarchy_bin <- file.path("~/Desktop/BiG-SCAPE/bin/hierarchy -n")

community_options <- "637268 -l -1 -v"
g_louvain <- run_louvain(X = "g", graph = g, community_bin = community_bin, convert_bin = convert_bin, hierarchy_bin = hierarchy_bin, community_options = community_options)


#com_l_0 <- make_clusters(g, as.numeric(g_louvain[[1]]), algorithm = 'Louvain', modularity = TRUE)
#com_l_1 <- make_clusters(g, as.numeric(g_louvain[[2]]), algorithm = 'Louvain', modularity = TRUE)
com_l_2 <- make_clusters(g, as.numeric(g_louvain[[3]]), algorithm = 'Louvain', modularity = TRUE)
#com_l_3 <- make_clusters(g, as.numeric(g_louvain[[4]]), algorithm = 'Louvain', modularity = TRUE)
# V(gx)$com_louvain_0 <- as.character(membership(com_l_0))
# V(gx)$com_louvain_1 <- as.character(membership(com_l_1))
# V(gx)$com_louvain_2 <- as.character(membership(com_l_2))
#V(g_mock)$com_louvain_3 <- as.character(membership(com_l_3))

V(g)$com <-  as.character(membership(com_l_2))
#V(g)$com <- as.character(membership(cluster_louvain(g, gamma = 1)))
# g_plot <- ggraph(g)  +
#   geom_edge_link(aes(color = sign, alpha = weight), show.legend = FALSE, width = 0.5) +
#   geom_node_point(aes(fill = com),shape = 21, color = "black" , show.legend = FALSE)  +
#   theme_graph()


counts <- dplyr::as_data_frame(g) %>%
  group_by(com) %>%
  dplyr::count() %>%
  dplyr::rename(n_mem = n)

g <- g %>%
  left_join(counts)

g1 <- g %>%
  select(-asv)

g.c <- contract.vertices(as.igraph(g1), V(as.igraph(g1))$com,
                         vertex.attr.comb=list(name= function(x) paste(x, collapse="|"),
                                               com= function(x) unique(x),
                                               n_mem = function(x) unique(x))
)

g.c <- igraph::simplify(g.c, edge.attr.comb = list(weight = function(x) median(x),
                                                   weight_orig =  function(x) median(x),
                                                   pvalue = function(x) median(x)), remove.loops = TRUE)

g <- g %>% activate(nodes) %>% mutate(com = as.character(paste("com", com, sep = "_")))
g.c <- as_tbl_graph(g.c) %>% activate(nodes) %>% mutate(com = as.character(paste("com", com, sep = "_")))

g.c <- as_tbl_graph(g.c) %>%
  activate(edges) %>%
  mutate(sign = ifelse(weight_orig < 0, "NEG", "POS"))

g <- as_tbl_graph(g) %>%
  activate(edges) %>%
  mutate(sign = ifelse(weight_orig < 0, "NEG", "POS"))

#write.graph(as.igraph(g), file = "~/Downloads/test_resfam_asv.graphml", format = "graphml")
write.graph(as.igraph(g.c), file = "osd2014_16S_asv/data/osd2014_dada_sparcc_modules.graphml", format = "graphml")



# Plot distribution of the communities in our samples ---------------------

my_db <- src_postgres(host = "localhost", port = 5432, dbname = "osd_analysis", options = "-c search_path=osd_analysis")


l_p <- transform_sample_counts(osd2014_dada2_phyloseq_beta, function (x) x/sum(x))


qpsmelt <- function(X) {
  if (taxa_are_rows(X)) {
    count_table <- as(otu_table(X), "matrix") %>%
      as_tibble(rownames = "OTU") %>%
      gather(label, Abundance, -OTU)
  }else{
    count_table <-as(otu_table(X), "matrix") %>%
      as_tibble(rownames = "label") %>%
      gather(OTU, Abundance, -label)
  }
  sample_table <- as(sample_data(X), "matrix") %>%
    as_tibble()
  taxa_table <- as(tax_table(X), "matrix") %>%
    as_tibble(rownames = "OTU")

  count_table %>%
    left_join(sample_table) %>%
    left_join(taxa_table)
}


l <- dplyr::as_data_frame(g %>% activate(nodes)) %>%
  left_join(qpsmelt(l_p) %>% rename(asv = OTU))

l$label <- factor(l$label, levels = st_100_order_terrestrial$label)

l_pl <- l %>%
  select(label, com, Abundance, Phylum) %>%
  #filter(counts >=5) %>%
  ggplot(aes(label, Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~com) +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Proportion") +
  xlab("OSD samples") +
  scale_fill_manual(values = randomcoloR::distinctColorPalette(k = l$Phylum %>% unique() %>% length()))

textcol <- "grey40"
l1 <- l %>%
  select(label, com, Abundance) %>%
  group_by(label, com) %>%
  summarise(N = sum(Abundance)) %>%
  ungroup()


l2 <- dplyr::as_data_frame(g %>% activate(nodes)) %>%
  left_join(qpsmelt(osd2014_dada2_phyloseq_beta) %>% rename(asv = OTU))



l3 <- l2 %>%
  ungroup() %>%
  select(com, Abundance, Phylum, Class, Order) %>%
  mutate(Phylum = paste("P:",Phylum,";O:", Order, sep = '')) %>%
  group_by(com, Phylum) %>%
  summarise( N = sum(Abundance)) %>%
  group_by(com) %>%
  mutate(prop = N/sum(N)) %>%
  ungroup() %>%
  mutate(Phylum = ifelse(prop > 0.005 & !(grepl('O:NA', Phylum)), Phylum, "Other")) %>%
  group_by(com, Phylum) %>%
  summarise(prop = sum(prop)) %>%
  ungroup() %>%
  tidyr::complete(com, Phylum)



com_colors <- c("#4D4D4D", #732
                "#F15854", #700
                "#60BD68", #548
                "#F17CB0", #410
                "#FAA43A", #431
                "#5DA5DA", #379
                "#B276B2", #351
                "#DECF3F", #311
                "#B2912F"  #29
)


#colors <- sample(com_colors)
#colors <- c("#DECF3F", "#60BD68", "#5DA5DA", "#4D4D4D", "#FAA43A", "#B2912F", "#F17CB0", "#F15854", "#B276B2")

l1 <- l1 %>% mutate(label = fct_relevel(label, st_100_order_terrestrial$label)) %>% group_by(label) %>% mutate(com = fct_reorder(com, desc(N))) %>% ungroup()
com_colors_order <- data_frame(r = seq(1,9,1), color = com_colors)

com_colors_order <- counts %>% ungroup %>% arrange(desc(n_mem)) %>% mutate(r = row_number(), com = paste0("com_", com)) %>%
  left_join(com_colors_order) %>% arrange(match(com, levels(l1$com)))

ggplot(l1, aes(label, N, fill = com)) +
  geom_col(color = "grey20", size = 0.2, width = 1) +
  scale_fill_manual(values = com_colors_order$color)  +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
        axis.ticks.x = element_blank(),
        #axis.title= element_text(size=8, face="bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        #legend.text = element_text(size = 8),
        #legend.key.size = unit(0.5, "cm"),
        legend.position = "top") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  ylab("Mixing proportions") +
  xlab("Samples") +
  labs(fill="")

ggsave(plot = last_plot(), filename = "osd2014_16S_asv/figures/osd2014_sparcc_com_barplot.pdf", width = 11.69, height = 8.27)

l_d <-  l2 %>%
  ungroup() %>%
  select(com, Abundance, Phylum, Class, Order) %>%
  mutate(Phylum = paste("P:",Phylum,";O:", Order, sep = '')) %>%
  group_by(com, Phylum) %>%
  summarise( N = sum(Abundance)) %>%
  group_by(com) %>%
  mutate(prop = N/sum(N)) %>%
  ungroup() %>% select(-N) %>% unique() %>% spread(com, prop, fill = 0) %>%
  column_to_rownames(var = "Phylum") %>%
  as.data.frame()


# Plot the occurrence of taxa in our communities --------------------------

l_hc <- hclust(vegan::vegdist(l_d, method = "jaccard"), method = "ward.D2")
l_hc_l <- hclust(vegan::vegdist(t(l_d), method = "jaccard"), method = "ward.D2")

# l_hc <- hclust(as.dist(cor(t(l_d),method = "spearman")), method = "average")
# l_hc_l <- hclust(1 - as.dist(cor((l_d),method = "pearson")), method = "complete")


#o <- l3 %>% group_by(Phylum) %>% summarise(L=sum(prop, na.rm = TRUE)) %>% arrange((L)) %>% .$Phylum
# o <- l3 %>% ungroup() %>% .$Phylum %>% unique() %>% as.character() %>% sort(decreasing = F)
# p <- l3 %>% group_by(com) %>% summarise(L=sum(prop, na.rm = TRUE)) %>% arrange(desc(L)) %>% .$com
#
# l3$Phylum <- factor(l3$Phylum, levels = o)
#l3$com <- factor(l3$com, levels = p)

l3$com <- factor(l3$com, levels = l_hc_l$labels[l_hc_l$order])
l3$Phylum <- factor(l3$Phylum, levels =l_hc$labels[l_hc$order])
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}
cols <- c(colorRampPalette(c("#e7f0fa", "#c9e2f6", "#95cbee", "#0099dc", "#4ab04a", "#ffd73e"))(10),
          colorRampPalette(c("#eec73a", "#e29421", "#e29421", "#f05336","#ce472e"), bias=9)(90))

p1 <- ggplot(l3 %>% mutate(prop = ifelse(is.na(prop), 0.00001, prop)), aes(x=Phylum,y=com, fill = prop))+
  geom_tile()+
  #redrawing tiles to remove cross lines from legend
  geom_tile(colour="gray30",size=0.3, show.legend = FALSE, na.rm = TRUE)+
  #remove axis labels, add title
  labs(x="",y="",title="")+
  #remove extra space
  scale_y_discrete(expand=c(0,0), position = "right")+
  #custom breaks on x-axis
  scale_x_discrete(expand=c(0,0))+
  #custom colours for cut levels and na values
  #scale_fill_gradientn(colours =rev(c("#d53e4f","#f46d43","#fdae61",
  # "#fee08b","#e6f598","#abdda4","#ddf1da")), na.value="grey90", trans = "log10", labels = scales::percent) +
  viridis::scale_fill_viridis(option = "B", direction = 1,trans = scales::log10_trans(), labels = scales::percent, breaks = base_breaks()) +
  #scale_fill_gradientn(colours = c(rev(viridis::viridis(100, begin = 0)), viridis::magma(100, begin = 0)), labels = scales::percent, trans = "log10") + #, name="Number of restaurant", guide = guide_legend( keyheight = unit(3, units= "mm"), keywidth=unit(12, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) ) +

  #labels=c("0k", "1k", "2k", "3k", "4k")) +
  #scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = scales::percent) +
  #mark year of vaccination
  #geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
  #equal aspect ratio x and y axis
  coord_fixed() +
  #set base size for all font elements
  theme_grey(base_size=10)+
  #theme options
  theme(
    legend.position = "top",
    #remove legend title
    legend.title=element_blank(),
    #remove legend margin
    legend.spacing = grid::unit(0,"cm"),
    #change legend text properties
    legend.text=element_text(colour=textcol,size=7,face="bold"),
    #change legend key height
    #legend.key.height=grid::unit(0.8,"cm"),
    #set a slim legend
    #legend.key.width=grid::unit(0.2,"cm"),
    #set x axis text size and colour
    axis.text.x=element_text(hjust = 1, vjust = 1, colour=textcol, angle = 45, size = 6),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.text.y = element_blank(),
    #set y axis text colour and adjust vertical justification
    #axis.text.y=element_text(vjust = 0.2,colour=textcol, size = 6),
    #change axis ticks thickness
    axis.ticks=element_line(size=0.4),
    #change title font, size, colour and justification
    plot.title=element_blank(),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank())


p2 <- g %>%
  tidygraph::activate(nodes) %>%
  left_join(results %>% select(OTU, sign) %>% unique() %>% rename(asv = OTU)) %>%
  dplyr::as_data_frame() %>%
  dplyr::select(com, sign) %>%
  group_by(com,  sign) %>%
  mutate(n = n()) %>%
  group_by(com) %>%
  unique() %>%
  mutate(N= sum(n), prop = n/N) %>%
  arrange(com) %>%
  ungroup() %>%
  mutate(sign = fct_relevel(sign, (c("Broad", "Non significant", "Narrow"))), com = fct_relevel(com, l_hc_l$labels[l_hc_l$order])) %>%
  ggplot(aes(com, prop, fill = sign)) +
  geom_col(width = 1, color = "grey20", size = 0.2) +
  ggpubr::rotate() +
  scale_fill_manual(values = c("#588157", "#9A9A9A", "#3891A6")) +
  theme_grey(base_size=10)+
  #theme options
  theme(
    legend.position = "right",
    #remove legend title
    legend.title=element_blank(),
    #remove legend margin
    legend.spacing = grid::unit(0,"cm"),
    #change legend text properties
    #legend.text=element_text(colour=textcol,size=7,face="bold"),
    #change legend key height
    #legend.key.height=grid::unit(0.8,"cm"),
    #set a slim legend
    #legend.key.width=grid::unit(0.2,"cm"),
    #set x axis text size and colour
    axis.text.x=element_text(colour=textcol, size = 6),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.text.y = element_blank(),
    #set y axis text colour and adjust vertical justification
    #axis.text.y=element_text(vjust = 0.2,colour=textcol, size = 6),
    #change axis ticks thickness
    axis.ticks=element_line(size=0.4),
    #change title font, size, colour and justification
    plot.title=element_blank(),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank()) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  xlab("")


ggpubr::ggarrange(p1, p2, nrow = 1, widths = c(0.7, 0.3), common.legend = F, ncol = 2)

ggsave(plot = last_plot(), filename = "osd2014_16S_asv/figures/osd2014_sparcc_com_heatmap.pdf", width = 11.69, height = 8.27)


# Plot the relationship between communities -------------------------------

#gc.1 <- g.c %>% activate(nodes) %>% mutate(wdegree = centrality_degree()/local_ave_degree(), name = fct_reorder(name, n_mem), com = fct_reorder(com, desc(n_mem)))

gc.1 <- g.c %>% activate(nodes) %>% mutate(wdegree = centrality_degree()/local_ave_degree(), name = fct_reorder(name, n_mem), com = fct_relevel(com, l_hc_l$labels[l_hc_l$order]))

gc.1 %>% dplyr::as_data_frame() %>% .$com %>% levels()

gc1_nodes <- gc.1 %>% activate(nodes) %>% dplyr::as_data_frame() %>% mutate(id = row_number())
gc1_edges <- gc.1 %>% activate(edges) %>% dplyr::as_data_frame()

gc1_edges$from <- plyr::mapvalues(gc1_edges$from, from = as.character(gc1_nodes$id), to = as.character(gc1_nodes$com))
gc1_edges$to <- plyr::mapvalues(gc1_edges$to, from = as.character(gc1_nodes$id), to = as.character(gc1_nodes$com))


gc1_edges$from <- factor(gc1_edges$from, levels = levels(gc1_nodes$com))
gc1_edges$to <- factor(gc1_edges$to, levels = levels(gc1_nodes$com))

gc1_edges_pos <- gc1_edges %>% filter(sign == "POS")
gc1_edges_neg <- gc1_edges %>% filter(sign == "NEG")


pg1 <- ggplot() +
  geom_curve(data = subset(gc1_edges_pos %>% arrange(weight) %>% head(1000), as.numeric(from) > as.numeric(to)),
             aes(x = from, xend = to, y = 0, yend = 0, alpha = weight),
             curvature = 0.5, ncp = 1000, lineend = 'butt', color = "black") +
  geom_curve(data = subset(gc1_edges_pos %>% arrange(weight) %>% head(1000), as.numeric(from) < as.numeric(to)),
             aes(x = from, xend = to, y = 0, yend = 0, alpha = weight),
             curvature = -0.5, ncp = 1000, lineend = 'butt', color = "black") +
  geom_point(data = gc1_nodes, aes(y = 0, x = com, fill = com), size = 10, shape = 21) +
  #scale_size_continuous(range = c(2,10)) +
  #geom_point(data = gc1_nodes, aes(y = 0, x = as.numeric(com))) +
  theme_bw() +
  scale_x_discrete(limits = levels(gc1_nodes$com)) +
  ylim(c(-0.5,1)) +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.ticks.x = element_blank(),
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4.5),
    axis.text.x = element_blank(),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top") +
  scale_fill_manual(values = com_colors_order %>% arrange(match(com, levels(gc1_nodes$com))) %>% .$color)

pg2 <- ggplot() +
  geom_curve(data = subset(gc1_edges_neg %>% arrange(weight) %>% head(1000), as.numeric(from) > as.numeric(to)),
             aes(x = from, xend = to, y = 0, yend = 0, alpha = weight),
             curvature = 0.5, ncp = 1000, lineend = 'butt',  color = "red") +
  geom_curve(data = subset(gc1_edges_neg %>% arrange(weight) %>% head(1000), as.numeric(from) < as.numeric(to)),
             aes(x = from, xend = to, y = 0, yend = 0, alpha = weight),
             curvature = -0.5, ncp = 1000, lineend = 'butt',  color = "red") +
  geom_point(data = gc1_nodes, aes(y = 0, x = com, fill = com), size = 10, shape = 21) +
  scale_size_continuous(range = c(2,10)) +
  #geom_point(data = gc1_nodes, aes(y = 0, x = as.numeric(com))) +
  theme_bw() +
  scale_x_discrete(limits = levels(gc1_nodes$com)) +
  ylim(c(-0.5,1)) +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=8),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=4.5),
    axis.ticks.x = element_blank(),
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4.5),
    axis.text.x = element_blank(),
    axis.title= element_text(size=8, face="bold"),
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "top") +
  scale_fill_manual(values = com_colors_order %>% arrange(match(com, levels(gc1_nodes$com))) %>% .$color)

ggpubr::ggarrange(pg1,pg2, nrow = 2, ncol = 1, common.legend = T)
ggsave(plot = last_plot(), filename = "osd2014_16S_asv/figures/osd2014_sparcc_com_arcplot.pdf", width = 11.69, height = 8.27)

# Calculate eigengenes
library(WGCNA)
mat <- (otu_table(osd2014_dada2_phyloseq_beta))
df_nodes <- dplyr::as_data_frame(g %>% activate(nodes)) %>%
  select(name, asv, com) %>%
  slice(match(colnames(mat), asv))
colnames(mat) <- df_nodes$name
PCs <- WGCNA::moduleEigengenes(mat, colors=df_nodes$com)
ME <- PCs$eigengenes
rownames(ME) <- rownames(mat)


osd2014_metadata <- as(sample_data(osd2014_dada2_phyloseq_beta), "matrix") %>% tbl_df()

osd2014_mld_adata <- tbl(my_db, "osd2014_mld_data") %>%
  collect(n = Inf) %>%
  rename(label = osd_id) %>%
  select(label, mld) %>%
  mutate_if(is.numeric, vegan::decostand, method = "standardize", na = na.omit)

osd2014_phenology_adata <- tbl(my_db, "osd2014_phenology_data") %>%
  collect(n = Inf) %>%
  rename(label = osd_id) %>%
  select(label, pp_8d_0.5) %>%
  mutate_if(is.numeric, vegan::decostand, method = "standardize", na = na.omit)

osd2014_satellite_adata <- tbl(my_db, "osd2014_satellite_data") %>%
  collect(n = Inf) %>%
  rename(label = osd_id) %>%
  select(label, chlor_a_8d, par_8d, poc_8d, KD490_8d) %>%
  mutate_if(is.numeric, vegan::decostand, method = "standardize", na = na.omit)

osd2014_woa13_adata <- tbl(my_db, "osd2014_woa13_data") %>%
  collect(n = Inf) %>%
  rename(label = osd_id) %>%
  select(-long_woa13, -lat_woa13) %>%
  mutate_if(is.numeric, vegan::decostand, method = "standardize", na = na.omit)

osd2014_iron_adata <- tbl(my_db, "osd2014_iron_data") %>%
  collect(n = Inf) %>%
  rename(label = osd_id) %>%
  select(label, iron) %>%
  mutate_if(is.numeric, vegan::decostand, method = "standardize", na = na.omit)

osd2014_rescaled_2013_median_long <- tbl(my_db, "osd2014_halpern_scaled_median") %>%
  collect(n = Inf) %>%
  filter(buffer == "1km") %>%
  select(-buffer) %>%
  spread(ohi_variable, median, fill =0) %>%
  select(-global_cumul_impact, -global_cumul_impact_diff_2008)


library(tidyverse)
library(broom)

osd2014_adata <- osd2014_mld_adata %>%
  left_join(osd2014_phenology_adata) %>%
  left_join(osd2014_satellite_adata) %>%
  left_join(osd2014_woa13_adata) %>%
  left_join(osd2014_iron_adata) %>%
  #left_join(osd2014_halpern_adata) %>%
  left_join(osd2014_rescaled_2013_median_long) %>%
  left_join(osd2014_sample_cohesion)


l <- osd2014_metadata %>% select(label, start_lat, start_lon) %>%
  #left_join(osd2014_halpern_adata)  %>%
  left_join(osd2014_adata) %>%
  gather(variable, value, -label) %>%
  right_join(ME %>% dplyr::as_data_frame() %>% rownames_to_column(var = "label") %>% gather(com, eigengene, -label)
  ) %>%
  mutate(value =  as.numeric(value), com = gsub("ME", "", com)) %>%
  group_by(com, variable) %>%
  do(tidy(Hmisc::rcorr(.$eigengene, .$value))) %>%
  ungroup() %>%
  filter(n > 110) %>%
  mutate(p.value.adj.bh = p.adjust(p.value, method = 'BH'),
         p.value.adj.holm = p.adjust(p.value, method = 'holm'),
         p.value.adj.fdr = p.adjust(p.value, method = 'fdr'))



l_pl <- l %>%
  ungroup %>%
  select(com, variable, estimate, contains("p.val")) %>%
  complete(com, variable) %>%
  mutate(estimate = ifelse(p.value.adj.bh < 0.01, estimate, NA), com = fct_relevel(com,  rev(l_hc_l$labels[l_hc_l$order]))) %>%
  group_by(variable) %>% mutate(nas = sum(is.na(estimate))) %>% filter(nas != 9) %>% ungroup()

l_pl$value<-cut(l_pl$estimate,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),include.lowest=TRUE,label=c("(-1,-0.75)","(-0.75,-0.5)","(-0.5,-0.25)","(-0.25,0)","(0,0.25)","(0.25,0.5)","(0.5,0.75)","(0.75,1)")) # this can be customized to put the correlations in categories using the "cut" function with appropriate labels to show them in the legend, this column now would be discrete and not continuous

l_pl_variable <- c(
  "start_lat",
  "start_lon",

  "temp", #WOD13
  "sal", #WOD13
  "silicate", #WOD13
  "oxy", #WOD13
  "phos", #WOD13

  "chlor_a_8d", #Satellite MODIS AQUA
  "KD490_8d", #Satellite MODIS AQUA
  "par_8d", #Satellite MODIS AQUA
  "poc_8d", #Satellite MODIS AQUA

  "pp_8d_0.5", #phenology SeaWiFs Sea-Viewing Wide Field-of-View Sensor

  "inorganic", #OHI
  "ocean_acidification", #OHI
  "ocean_pollution", #OHI
  "plumes_fert", #OHI
  "shipping", #OHI
  "slr", #OHI
  "sst", #OHI
  "uv",#OHI
  "cohesion_negative",
  "cohesion_positive"
)

  ggplot(l_pl %>% mutate(variable = fct_relevel(variable, rev(l_pl_variable))), aes(com, variable, fill = value)) +
    geom_tile()+
    #redrawing tiles to remove cross lines from legend
    geom_tile(colour="grey60",size=0.3, show.legend = FALSE, na.rm = TRUE)+
    #remove axis labels, add title
    labs(x="",y="",title="")+
    #remove extra space
    scale_y_discrete(expand=c(0,0), position = "left")+
    #custom breaks on x-axis
    scale_x_discrete(expand=c(0,0))+
    #custom colours for cut levels and na values
    #scale_fill_gradientn(colours = (c("#d53e4f","#f46d43",
    # "#fee08b","#e6f598","#ddf1da","#abdda4")), na.value="#000004", limits = c(-0.5,0.75)) +
    scale_fill_brewer(palette = "RdYlGn",name="Correlation",  na.value="#000004") +
    #scale_fill_gradientn(colours = c("#D61818","#FFAE63","#FFFFBD","#B5E384"), na.value="#000004", limits = c(-0.8,0.8)) +
    #scale_fill_gradient2(low = "#d53e4f", high = "#1F629A", mid = "ivory2", midpoint = 0,  na.value="#000004", limits = c(-0.5,0.75), ) +
    #scale_fill_gradientn(colours = c(rev(viridis::viridis(100, begin = 0)), viridis::magma(100, begin = 0)), labels = scales::percent, trans = "log10") + #, name="Number of restaurant", guide = guide_legend( keyheight = unit(3, units= "mm"), keywidth=unit(12, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) ) +

    #labels=c("0k", "1k", "2k", "3k", "4k")) +
    #scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,"YlGnBu")),na.value="grey90", trans = scales::log_trans(), breaks = base_breaks(), labels = scales::percent) +
    #mark year of vaccination
    #geom_vline(aes(xintercept = 36),size=3.4,alpha=0.24)+
    #equal aspect ratio x and y axis
    coord_fixed() +
    #set base size for all font elements
    theme_grey(base_size=10)+
    #theme options
    theme(
      legend.position = "top",
      #remove legend title
      legend.title=element_blank(),
      #remove legend margin
      legend.spacing = grid::unit(0,"cm"),
      #change legend text properties
      legend.text=element_text(colour=textcol,size=7,face="bold"),
      #change legend key height
      #legend.key.height=grid::unit(0.8,"cm"),
      #set a slim legend
      #legend.key.width=grid::unit(0.2,"cm"),
      #set x axis text size and colour
      axis.text.x=element_text(hjust = 1, vjust = 1, colour=textcol, angle = 45, size = 6),
      #axis.text.x = element_blank(),
      #axis.ticks.x = element_blank(),
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      #set y axis text colour and adjust vertical justification
      #axis.text.y=element_text(vjust = 0.2,colour=textcol, size = 6),
      #change axis ticks thickness
      axis.ticks=element_line(size=0.4),
      #change title font, size, colour and justification
      plot.title=element_blank(),
      #remove plot background
      plot.background=element_blank(),
      #remove plot border
      panel.border=element_blank())

ggsave(plot = last_plot(), filename = "osd2014_16S_asv/figures/osd2014_sparcc_com_cor.pdf", width = 11.69, height = 8.27)

save.image("osd2014_16S_asv/data/osd2014_16S_asv_networks_results.Rdata")
