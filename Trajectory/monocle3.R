ti_monocle3 <- function(cds, parameters){
  suppressPackageStartupMessages(require(SingleCellExperiment))
  method <- names(parameters)
  if(parameters$preprocess_method != 'seurat'){
    cds <- monocle3::preprocess_cds(
      cds,
      method = parameters$preprocess_method,
      num_dim = parameters$num_dim,
      norm_method = parameters$norm_method
    )
  } else {
    parameters$preprocess_method <- 'PCA'
  }
  if(parameters$reduction_method != 'seurat'){
    # perform dimensionality reduction
    cds <- monocle3::reduce_dimension(
      cds,
      preprocess_method = parameters$preprocess_method,
      max_components = parameters$max_components,
      reduction_method = parameters$reduction_method,
      cores = parameters$cores
    )
  } else {
    parameters$reduction_method <- "UMAP"
  }
  
  # perform clustering
  cds <- monocle3::cluster_cells(
    cds,
    k = parameters$k,
    cluster_method = parameters$cluster_method,
    reduction_method = parameters$reduction_method,
    random_seed = parameters$random_seed
  )
  
  # calculate trajectory
  cds <- monocle3::learn_graph(cds, verbose = FALSE)
  
  # Only UMAP
  if(parameters$reduction_method == 'UMAP'){
    cds <- monocle3::order_cells(cds, reduction_method = parameters$reduction_method ,root_pr_nodes = colnames(cds@principal_graph_aux$UMAP$dp_mst))
  }
  
  S4Vectors::metadata(cds)$start_cell <- parameters[['start_cell']]
  return(cds)
}

wrapper_monocle3 <- function(cds){
  dimred <- SingleCellExperiment::reducedDim(cds,"UMAP") %>%
    magrittr::set_colnames(c("comp_1", "comp_2"))
  dimred_milestones <- t(cds@principal_graph_aux$UMAP$dp_mst) %>%
    magrittr::set_colnames(colnames(dimred))
  milestone_network <-
    igraph::as_data_frame(cds@principal_graph$UMAP) %>%
    dplyr::transmute(
      from,
      to,
      length = sqrt(rowSums((dimred_milestones[from, ] - dimred_milestones[to, ])^2)),
      directed = FALSE
    )
  
  milestone_percentages <-
    cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex %>%
    magrittr::set_colnames("index") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::transmute(
      cell_id,
      milestone_id = rownames(dimred_milestones)[index],
      percentage = 1
    )
  
  dimred_segment_progressions <-
    milestone_network %>%
    dplyr::select(from, to) %>%
    dplyr::mutate(percentage = purrr::map(seq_len(dplyr::n()), ~ c(0, 1))) %>%
    tidyr::unnest(percentage)
  
  dsp_names <-
    dimred_segment_progressions %>%
    {ifelse(.$percentage == 0, .$from, .$to)}
  dimred_segment_points <- dimred_milestones[dsp_names, , drop = FALSE]
  
  test <-
    dynwrap::wrap_data(
      cell_ids = milestone_percentages$cell_id
    ) %>%
    dynwrap::add_trajectory(
      milestone_ids = rownames(dimred_milestones),
      milestone_network = milestone_network,
      milestone_percentages = milestone_percentages
    ) %>%
    dynwrap::add_dimred(
      dimred = dimred,
      dimred_milestones = dimred_milestones,
      dimred_segment_progressions = dimred_segment_progressions,
      dimred_segment_points = dimred_segment_points
    ) %>% 
    dynwrap::add_expression(counts = Matrix::t(SingleCellExperiment::counts(cds)) %>% as("dgCMatrix") ,
                            expression = Matrix::t(SingleCellExperiment::logcounts(cds)) %>% as("dgCMatrix")) %>%
    dynwrap::add_root(root_cell_id = S4Vectors::metadata(cds)$start_cell) %>%
    dynwrap::simplify_trajectory()
  
    return(test)
}


find_nearest_edge <- function(cds){
  Z <- t(SingleCellExperiment::reducedDims(cds)[['UMAP']])
  Y <- cds@principal_graph_aux[['UMAP']]$dp_mst
  dp_mst <- monocle3::principal_graph(cds)[['UMAP']]
  closest_vertex <- cds@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_closest_vertex"]]
  closest_vertex_names <- colnames(Y)[closest_vertex[, 1]]
  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))
  nearest_edges <- matrix(rep(0, length(Z[1:2, ])), ncol = 2)
  
  for(i in 1:length(closest_vertex)) { # This loop is going to be slow
    neighbors <- names(igraph::neighborhood(dp_mst,
                                            nodes = closest_vertex_names[i],
                                            mode = 'all')[[1]])[-1]
    projection <- NULL
    distance <- NULL
    Z_i <- Z[, i]
    for(neighbor in neighbors) {
      tmp <- project_point_to_line_segment(Z_i, Y[, c(closest_vertex_names[i],
                                                      neighbor)])
      if(any(is.na(tmp))) {
        tmp <- Y[, neighbor]
      }
      projection <- rbind(projection, tmp)
      distance <- c(distance, stats::dist(rbind(Z_i, tmp)))
    }
    if(class(projection) != 'matrix') {
      projection <- as.matrix(projection)
    }
    
    which_min <- which.min(distance)
    
    if(length(which_min) == 0)
      browser()
    
    nearest_edges[i, ] <- c(closest_vertex_names[i], neighbors[which_min])
  }
  rownames(nearest_edges) <- rownames(closest_vertex)
  
  return(nearest_edges)
}


project_point_to_line_segment <- function(p, df){
  # returns q the closest point to p on the line segment from A to B
  A <- df[, 1]
  B <- df[, 2]
  # vector from A to B
  AB <- (B-A)
  # squared distance from A to B
  AB_squared = sum(AB^2)
  if(AB_squared == 0) {
    # A and B are the same point
    q <- A
  }
  else {
    # vector from A to p
    Ap <- (p-A)
    # from http://stackoverflow.com/questions/849211/
    # Consider the line extending the segment, parameterized as A + t (B - A)
    # We find projection of point p onto the line.
    # It falls where t = [(p-A) . (B-A)] / |B-A|^2
    # t <- max(0, min(1, sum(Ap * AB) / AB_squared))
    t <- sum(Ap * AB) / AB_squared
    
    if (t < 0.0) {
      # "Before" A on the line, just return A
      q <- A
    }
    else if (t > 1.0) {
      # "After" B on the line, just return B
      q <- B
    }
    else {
      # projection lines "inbetween" A and B on the line
      q <- A + t * AB#
    }
  }
  return(q)
}


data_prep_monocle3 <- function(sobj, path, filename, reductions, seed, cell_identities, start_cell){
  suppressPackageStartupMessages(require(SingleCellExperiment))
  fs::path(path,"monocle3") %>% fs::dir_create()
  expression <- sobj[['RNA']]@counts
  # expression <- expression[,Matrix::colSums(expression) != 0]
  
  gene_metadata <- Matrix::rowSums(expression) %>% 
    tibble::enframe(name = "gene_short_name", value = "num_cells_expressed") %>% 
    dplyr::mutate(rowname = gene_short_name) %>% 
    tibble::column_to_rownames()
  
  # Create a cds object
    cds <- monocle3::new_cell_data_set(expression_data = expression, cell_metadata = sobj@meta.data, gene_metadata = gene_metadata) 
    SingleCellExperiment::logcounts(cds) <- log(expression+1)
  monocle3_default_parameters <- list("num_dim" = 50L,
                                      "preprocess_method" = 'PCA',
                                      "norm_method" = 'log',
                                      "max_components" = 2L,
                                      "reduction_method" = 'UMAP',
                                      "cores" = 1L,
                                      "k" = 20L,
                                      "cluster_method" = 'leiden',
                                      "random_seed" = seed,
                                      "start_cell" = start_cell)
  
  # Add the different dimensions reduction
  monocle3_params <- list('default' = monocle3_default_parameters)
  available.reductions <- purrr::map_chr(reductions,~stringr::str_extract(names(sobj@reductions), .x) %>% purrr::discard(is.na))
  for(red in available.reductions){
    if(stringr::str_detect(red,'(?i)pca')){
      redname_monocle3 <- "PCA"
      temp <- monocle3_default_parameters
      temp[['preprocess_method']] <- 'seurat'
      monocle3_params[['seurat_pca']] <- temp
    } 
    if(stringr::str_detect(red,'(?i)umap')){
      redname_monocle3 <- "UMAP"
      temp <- monocle3_default_parameters
      temp[['preprocess_method']] <- 'seurat'
      temp[['reduction_method']] <- 'seurat'
      monocle3_params[['seurat_umap']] <- temp
    }
    SingleCellExperiment::reducedDim(cds, redname_monocle3) <- sobj@reductions[[red]]@cell.embeddings
  }
  
  
  # SingleCellExperiment::reducedDim(cds, "UMAP") <- sobj@reductions[['umap']]@cell.embeddings 
  # SingleCellExperiment::reducedDim(cds, "PCA") <- sobj@reductions[['pca']]@cell.embeddings 
  saveRDS(object = cds, file = fs::path(path,"monocle3",paste0(filename,".rds")))
  # write parameters list
  yaml::write_yaml(monocle3_params,file = fs::path(path,"monocle3","monocle3_params.yaml"))
  #grouping
  sobj[[cell_identities]] %>% 
    tibble::as_tibble(rownames = "cell_id") %>% 
    dplyr::rename(group_id = {{cell_identities}}) %>%
    saveRDS(file = fs::path(path,"monocle3",paste0(filename,"_grouping.rds")))
}

run_monocle3 <- function(path,filename){
  path %<>% fs::path("monocle3")
  cds <- readRDS(fs::path(path,paste0(filename,".rds")))
  monocle3_params <- yaml::read_yaml(fs::path(path,"monocle3_params.yaml")) %>%
    filter_computed_traj(path)
  
  create_trajectory_dir(monocle3_params,path)
  groups <- readRDS(fs::path(path,paste0(filename,"_grouping.rds")))
  
  trajectory_cds <- purrr::map(monocle3_params, ~ ti_monocle3(cds = cds, parameters = .))
  
  # purrr::iwalk(trajectory_cds, ~ saveRDS(.x,file = fs::path(path,.y, "trajectory_result.rds")))
  output <- purrr::map(trajectory_cds,wrapper_monocle3) %>%
    purrr::map(dynwrap::add_grouping, grouping = groups)
  return(output)
}

