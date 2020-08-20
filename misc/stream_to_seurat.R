get_node_coordinates <- function(node = NULL, file.h5 = NULL,all_tree,root){
  if(is.null(file.h5)){
    stop("No H5 file were provided")
  }
  if(is.null(node)){
    stop("No node were provided, cannot extract coordinate")
  }
  
  coordinate_accessor <- glue::glue(all_tree[stringr::str_which(all_tree,root)],"/nodes/",node)
  coordinate <- file.h5[[coordinate_accessor]][]
  names(coordinate) <- c("x","y")
  return(coordinate)
}

compute_trajectory <- function(i,sobj,name,timeseq){
  `%>%` <- magrittr::`%>%`
  cells <- tibble::rownames_to_column(sobj[[]],var = "barcode") %>% dplyr::filter({{name}} < timeseq[i+1] & {{name}} >= timeseq[i]) %>% dplyr::pull(barcode)
  sobj@reductions$umap@cell.embeddings %>% tibble::as_tibble(rownames = "barcode") %>% dplyr::filter(barcode %in% cells) %>% dplyr::select(-barcode) %>% dplyr::summarise_all(mean)
}

add_trajectory_umap <- function(p,data_traj,my.seed){
  p <- p+ggplot2::geom_path(data = data_traj, mapping = ggplot2::aes(x=UMAP_1,y=UMAP_2,group = edge),arrow = ggplot2::arrow(ends = 'last', type = "closed", angle = 15))+
    ggrepel::geom_text_repel(data = data_traj, mapping = ggplot2::aes(x=UMAP_1,y=UMAP_2,label =labels, group = edge),seed = my.seed) + ggplot2::ggtitle("")
  return(p)
}

traj_by_edge <- function(edge,node_time,sobj, tree_path, pseudotime_accessor_name,nbins = 10){
  `%>%` <- magrittr::`%>%`
  
  edge_num <- edge %>% stringr::str_remove_all("\\(|\\)") %>% 
    stringr::str_split(',') %>% 
    magrittr::extract2(1) %>%
    stringr::str_remove_all(" ")
  node_time_sub <- node_time %>% 
    dplyr::filter(node == edge_num[1] | node == edge_num[2])
  
  node <- FetchData(object = sobj, vars = "nodes") %>% 
    dplyr::pull(nodes)

  sobj_sub <- sobj[, which(x = node %in% tree_path[[edge]])]
  
  expr <- FetchData(object = sobj, vars = pseudotime_accessor_name)
  sobj_sub <- sobj_sub[, which(x = expr >= node_time_sub$pseudotime[1] & expr < node_time_sub$pseudotime[2])]
  
  min_time <- min(sobj_sub[[pseudotime_accessor_name]],na.rm = TRUE)
  max_time <- max(sobj_sub[[pseudotime_accessor_name]],na.rm = TRUE)
  
  timeseq <- seq(from=min_time, to = max_time, length.out = nbins)
  edge_name <- node_time_sub %>% 
    dplyr::pull(label) %>% 
    glue::glue_collapse(sep='_')
  
  traj <- purrr::map_df(seq_along(timeseq)[-nbins],compute_trajectory, sobj = sobj_sub, name = !!rlang::sym(pseudotime_accessor_name), timeseq = timeseq) %>% 
    dplyr::mutate(edge = edge_name, pseudotime = timeseq[-nbins]) %>% 
    dplyr::mutate(labels = dplyr::case_when(pseudotime == max(pseudotime) ~ dplyr::pull(node_time_sub,label)[2],                                                                         pseudotime == min(pseudotime) ~ dplyr::pull(node_time_sub,label)[1]))
}

plot_traj_umap <- function(root = 'S0',sobj = NULL ,features = NULL, my.seed = 1337,combine = TRUE, path = NULL, name = NULL){
  `%>%` <- magrittr::`%>%`
  
  if(is.null(sobj)){
    stop("Please provide a Seurat object")
  }
  if(is.null(path)){
    stop("No path provided do not know where to lookup for files")
  }
  if(!dir.exists(path)){
    stop("Specified directory does not exist")
  }
  if(is.null(name)){
    stop("No name provided, cannot search for trajectory files")
  }
  if(any(!file.exists(c(paste0(path,'/',name,'.h5ad'),paste0(path,'/tree_',name,'.json'), paste0(path,'/',root,'_tree_path_',name,'.json'))))){
    stop("One of the requested files does not exist")
  }
  
  file.h5 <- hdf5r::H5File$new(paste0(path,'/',name,'.h5ad'), mode="r+")
  tree_node_label <- jsonlite::fromJSON(txt = paste0(path,'/tree_',name,'.json')) %>% 
    tibble::as_tibble() %>% 
    tidyr::gather(key = 'node',value="label")
  
  tree_path <-  jsonlite::fromJSON(txt = paste0(path,'/',root,'_tree_path_',name,'.json')) 
  # GEt all pseudotime frome different root (allow user to change root)
  all_pseudotime <- stringr::str_subset(file.h5$ls(recursive=TRUE)$name,'pseudotime')
  all_tree <- stringr::str_subset(file.h5$ls(recursive=TRUE)$name,'uns/subwaymap_S\\d+$')
  
  if(!any(stringr::str_detect(all_tree,root))){
    stop(paste("Tree with root",root,"is not available"))
  }
  
  if(!any(stringr::str_detect(all_pseudotime,root))){
    stop(paste("Pseudotime from root",root,"is not available"))
  }
  
  nodes_accessor <- glue::glue(all_tree[stringr::str_which(all_tree,root)],"/nodes")
  edges_accessor <- glue::glue(all_tree[stringr::str_which(all_tree,root)],"/edges")
  pseudotime_accessor <-glue::glue("obs/",root,"_pseudotime")
  pseudotime_accessor_name <- make.names(pseudotime_accessor)
  
  if(is.null(features)){
    features = pseudotime_accessor_name
  }
  
  # Extract nodes information from the tree 
  nodes <- file.h5[[nodes_accessor]]$ls()$name
  names(nodes) <- nodes
  
  edges <- file.h5[[edges_accessor]]$ls()$name
  
  node_time <- purrr::map(nodes,get_node_coordinates,file.h5 = file.h5,all_tree = all_tree, root = root) %>% 
    purrr::map('x') %>% tibble::as_tibble() %>% 
    tidyr::gather(key = 'node', value = 'pseudotime') %>% 
    dplyr::arrange(pseudotime) %>% 
    dplyr::inner_join(tree_node_label,by='node')
  
  # Add pseudotime from root in sobj metadata
  cell_with_pseudotime_node <- tibble::tibble(barcode = file.h5[["obs/_index"]][], pseudotime = file.h5[[pseudotime_accessor]][],node =file.h5[["obs/node"]][])
  
  if(!all(dplyr::pull(cell_with_pseudotime_node,barcode) %in% rownames(sobj[[]]))){
    stop("Cell barcodes does not match. Did you provide the Seurat object corresponding to this trajectory")
  }
  
  cell_with_pseudotime_node <- tibble::rownames_to_column(sobj[[]],var = "barcode") %>% dplyr::select(barcode) %>% dplyr::left_join(cell_with_pseudotime_node) 
  
  file.h5$close_all()
  
  sobj[[pseudotime_accessor_name]] <- dplyr::pull(cell_with_pseudotime_node,pseudotime)
  sobj[["nodes"]] <- dplyr::pull(cell_with_pseudotime_node,node)
  
  traj <- purrr::map_df(edges,traj_by_edge, node_time = node_time, sobj = sobj, tree_path = tree_path,pseudotime_accessor_name = pseudotime_accessor_name)
  
  p <- Seurat::FeaturePlot(sobj,features = features, cols=c("lightgrey","blue"),combine = FALSE)
  p <- purrr::map(p,add_trajectory_umap, data_traj = traj, my.seed = my.seed)
  
  if(combine & length(features) > 1){
    p <- Seurat::CombinePlots(p)
  }
  
  return(p)
}

# load("/mnt/Data/Project/NF1-MPNST/NF1-MPNST_20190304/ANALYSIS/DEFAULT/NF1-MPNST_20190304_SCT_DEFAULT_pca60_resolution0.7.rda")
load("/mnt/Data/Project/NF1-MPNST/NF1-MPNST_20190212/ANALYSIS/DEFAULT/NF1-MPNST_20190212_SCT_DEFAULT_pca100.rda")

plot_traj_umap(root = 'S0', sobj = sobj, features = NULL, my.seed = 1337, combine = FALSE, path = "~/stream/stream_result_grafted/",name = 'grafted')
  ## lets look at the content
#file.h5$ls(recursive=TRUE)$name





# Cell pseudotime
# temp <- file.h5[["obs/S0_pseudotime"]]
# temp[]

# Get nodes from a tree (only important one : start, end, branching)
# file.h5[["uns/subwaymap_S3/nodes"]]$ls()$name

# Get edges from a tree 
# file.h5[["uns/subwaymap_S3/edges"]]$ls()$name

# Get coordinate of edges : Line 1 = pseudotime 
# file.h5[["uns/subwaymap_S3/edges/(81, 2)"]][1,]


# Cell barcode
# file.h5[["obs/_index"]][]

# load("/mnt/Data/Project/NF1-MPNST/NF1-MPNST_20190304/ANALYSIS/DEFAULT/NF1-MPNST_20190304_SCT_DEFAULT_pca60_resolution0.7.rda")


