
# Libraries for the fuzzyforests ------------------------------------------

brStick <- function (x) {
  m <- 0
  out <- matrix(NA, ncol = 2, nrow = length(x))
  colnames(out) <- c("% of Variability", "B-Stick Threshold")
  for (i in 1:length(x)) {
    for (k in i:length(x)) {
      m <- m + ((1 / length(x)) * (1 / k))
    }
    out[i, ] <- c((x[i] / sum(x)) * 100, m * 100)
    m <- 0
  }
  return(list("Use PCs:" = which(out[, 1] > out[, 2]), Table = out))
}


get_ftable_bs <- function(X){
  f <- ff_fit_r[[X]]$feature_list %>% as_tibble %>% mutate(run = names(ff_fit_r[X]))
  f1 <- tibble(pos = brStick(f$variable_importance)[["Use PCs:"]]) %>%
    mutate(pos_difference = pos- lag(pos), pos_difference = ifelse(is.na(pos_difference), 1, pos_difference)) %>%
    filter(pos_difference > 1)
  if (nrow(f1) == 0) {
    f[brStick(f$variable_importance)[["Use PCs:"]], ]
  }else{
    head(f, f1[1,]$pos - 1)
  }
}

get_ftable <- function(X){
  f <- ff_fit_r[[X]]$feature_list %>% as_tibble %>% mutate(run = names(ff_fit_r[X]))
}

get_predict <- function(X){
  ff_predict[[X]] %>% as_tibble(rownames = "label") %>% mutate(run = names(ff_predict[X]), label =  names(ff_predict[[X]]))
}

get_accuracy <- function(X){
  ff_accuracy[[X]] %>% as_tibble() %>% mutate(run = names(ff_predict[X]))
}


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


get_g <- function(X){
  g <- igraph::induced_subgraph(g, vids = top_features %>% filter(com == X) %>% .$name)
  if (vcount(g) > 0){
    #weights <- unique(sort(E(g)$weight, decreasing = FALSE))
    #g <- binary_search_filter(g = g, low = 1, high = length(weights))$graph
    g %>%
      as_tbl_graph() %>%
      activate(nodes) %>%
      mutate(weighted_degree = centrality_degree() / local_ave_degree()) %>%
      inner_join(as(tax_table(osd2014_dada2_phyloseq_beta_filt), "matrix") %>% as_tibble(rownames = "asv"))
  }
}

get_g_ohi <- function(X){
  g <- igraph::induced_subgraph(g, vids = top_features %>% filter(com == X) %>% .$name)
  if (vcount(g) > 0){
    #weights <- unique(sort(E(g)$weight, decreasing = FALSE))
    #g <- binary_search_filter(g = g, low = 1, high = length(weights))$graph
    g %>%
      as_tbl_graph() %>%
      activate(nodes) %>%
      mutate(weighted_degree = centrality_degree() / local_ave_degree()) %>%
      inner_join(as(tax_table(Y), "matrix") %>% as_tibble(rownames = "asv"))
  }
}

plot_com <- function(X){
  all_comps %>%
    inner_join(osd2014_cdata) %>%
    group_by(label, Order, com, meow_province) %>%
    dplyr::summarise(prop = sum(prop)) %>%
    ungroup() %>%
    filter(prop >0) %>%
    group_by(Order) %>%
    mutate(agg_prop = sum(prop)) %>% ungroup() %>% mutate(Order = ifelse(agg_prop > 0.01, Order, "Other")) %>%
    ungroup() %>%
    filter(!(com %in% c("com_8")))%>%#, meow_province %in% c("Mediterranean Sea", "Northern European Seas")) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial $label)) %>%
    mutate(meow_province = fct_relevel(meow_province, meow_provinces)) %>%
    select(label, prop, Order, meow_province, com) %>%
    filter(meow_province == X) %>%
    droplevels() %>%
    complete(label,com, meow_province, Order, fill = list(prop = 0)) %>%
    filter(!is.na(com)) %>% #filter(meow_province == "Cold Temperate Northwest Atlantic") %>% View
    inner_join(all_comps_order %>% rename(Order = order_mod)) %>%
    ggplot(aes(label, prop, fill = Order)) +
    geom_col(color = "grey20", size = 0.2, width = 1) +
    facet_wrap(meow_province~com , scales = "free_x", ncol = 5) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    xlab("") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(limits = all_comps_order$order_mod, values = all_comps_order$colour)
}

plot_com_ohi <- function(X){
  all_comps %>%
    inner_join(osd2014_cdata) %>%
    group_by(label, Order, com, meow_province) %>%
    dplyr::summarise(prop = sum(prop)) %>%
    ungroup() %>%
    filter(prop >0) %>%
    group_by(Order) %>%
    mutate(agg_prop = sum(prop)) %>% ungroup() %>% mutate(Order = ifelse(agg_prop > 0.1, Order, "Other")) %>%
    ungroup() %>%
    filter(!(com %in% c("com_8")))%>%#, meow_province %in% c("Mediterranean Sea", "Northern European Seas")) %>%
    mutate(label = fct_relevel(label, osd2014_order_terrestrial$label)) %>%
    mutate(meow_province = fct_relevel(meow_province, meow_provinces)) %>%
    select(label, prop, Order, meow_province, com) %>%
    filter(meow_province == X) %>%
    droplevels() %>%
    complete(label,com, meow_province, Order, fill = list(prop = 0)) %>%
    filter(!is.na(com)) %>% #filter(meow_province == "Cold Temperate Northwest Atlantic") %>% View
    inner_join(all_comps_order %>% rename(Order = order_mod)) %>%
    ggplot(aes(label, prop, fill = Order)) +
    geom_col(color = "grey20", size = 0.2, width = 1) +
    facet_wrap(meow_province~com , scales = "free_x", ncol = 5) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    xlab("") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(limits = all_comps_order$order_mod, values = all_comps_order$colour)
}
# Plotting ----------------------------------------------------------------

# somewhat hackish solution to:
# https://twitter.com/EamonCaddigan/status/646759751242620928
# based mostly on copy/pasting from ggplot2 geom_violin source:
# https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r

library(ggplot2)
library(dplyr)


"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)

            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
          },

  draw_group = function(data, panel_scales, coord) {
    # Find the points for the line to go all the way around
    data <- transform(data, xminv = x,
                      xmaxv = x + violinwidth * (xmax - x))

    # Make sure it's sorted properly to draw the outline
    newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                     plyr::arrange(transform(data, x = xmaxv), -y))

    # Close the polygon: set first and last point the same
    # Needed for coord_polar and such
    newdata <- rbind(newdata, newdata[1,])

    ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
  },

  draw_key = draw_key_polygon,

  default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                    alpha = NA, linetype = "solid"),

  required_aes = c("x", "y")
  )


### Example:
ggplot(diamonds, aes(cut, carat)) +
  geom_flat_violin() +
  coord_flip()


# Fast phyloseq melt ------------------------------------------------------

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
    left_join(taxa_table) %>%
    filter(Abundance > 0)
}
