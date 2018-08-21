run_louvain <- function(X, graph = graph, convert_options = NULL, community_options = NULL,
         community_bin = community_bin, convert_bin = convert_bin, hierarchy_bin = hierarchy_bin){

  # convert names to numeric
  V(graph)$id <- V(graph)$name
  V(graph)$name <- as.factor(V(graph)$name)

  # Get convert options

  if (is.null(convert_options)){
    convert_options <-""
  }else{
    convert_options <- convert_options
  }

  # Get community options
  if (is.null(community_options)){
    community_options <-"12345 -l -1 -v"
  }else{
    community_options <- community_options
  }

  # convert files
  bin_file <- paste("~/Downloads/", X,"-louvain.bin", sep = "")
  w_file <- paste("~/Downloads/", X,"-louvain_w.bin", sep = "")

  # community files
  tree_file <- paste("~/Downloads/", X,"-louvain.tree", sep = "")
  # hierarchy files
  g_name <- paste("~/Downloads/", X,".tsv", sep = "")

  convert_command <- paste(convert_bin, "-i", g_name, "-o", bin_file, "-w", w_file ,convert_options, sep = " ")
  community_command <- paste(community_bin, bin_file, community_options,"-w", w_file, ">", tree_file, sep = " ")

  # export file
  write.graph(graph, g_name, format = "ncol")

  cat("Command called: ",convert_command,"\n")
  system(convert_command)
  cat("Command called: ", community_command,"\n")
  system(community_command)

  # guess how many levels do we have:
  hierarchy_command <- paste(hierarchy_bin, tree_file, "-n", sep = " ")
  cat("Command called: ", hierarchy_command,"\n")
  louvain_hier <- system(hierarchy_command, intern = TRUE)
  louvain_levels <- as.numeric(gsub(pattern = "Number of levels: ", "", louvain_hier[[1]]))

  modules <- vector(mode = "list")
  mod_louvain <- vector(mode = "list")

  for (i in 1:louvain_levels){
    com_file <- paste("~/Downloads/", X,"-louvain_", i, ".com", sep = "")

    hierarchy_command <- paste(hierarchy_bin, tree_file, "-l", i - 1, ">", com_file, sep = " ")
    system(hierarchy_command)
    modules[[i]] <- read.table(file = com_file, fill = TRUE, stringsAsFactors = F, sep = " ", header = T) %>%
      magrittr::set_colnames(c("node", "com"))

    cat(paste("Parsing Louvain communities output for level", i - 1, ":", length(unique(modules[[i]]$com)), "communities detected.\n"))
    mod_louvain[[i]] <- modules[[i]][match(get.vertex.attribute(graph, "name"), modules[[i]]$node),] %>%
      dplyr::select(com) %>%
      .$com
  }

  #map_file <- paste(out_folder, "/", X, ".map", sep = "")

  return(mod_louvain)
}
