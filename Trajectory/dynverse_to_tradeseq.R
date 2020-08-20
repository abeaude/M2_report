dynverse_to_tradeseq <- function(trajectory) {
  # Retrieve graph from milestone network
  gr <- igraph::graph_from_data_frame(
    d = trajectory$milestone_network %>% dplyr::rename(weight = length),
    directed = any(trajectory$milestone_network$directed),
    vertices = trajectory$milestone_ids
  )

  # To which milestone each cell belongs to
  y_to_cells <-
    trajectory$milestone_percentages %>%
    dplyr::group_by(cell_id) %>%
    dplyr::slice(which.max(percentage)) %>%
    dplyr::mutate(percentage = 1) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::transmute(Y = which(milestone_id == igraph::V(gr)$name), cells = cell_id) %>%
    dplyr::ungroup()

  root <- trajectory$root_milestone_id

  endpoints <- names(which(igraph::degree(gr) == 1))
  endpoints <- endpoints[!endpoints %in% root]
  if (rlang::is_empty(endpoints)) {
    endpoints_w_root <- names(which(igraph::degree(gr) == 1))
    no_endpoints <- rlang::is_empty(endpoints_w_root)
    is_root <- root == endpoints_w_root
    mess <- dplyr::case_when(
      no_endpoints ~ "No vertex have a degree 1. Can't find any endpoints",
      is_root ~ stringr::str_c("The only endpoints is ", root, " which is your root.\n Please consider changing your root")
    )
    rlang::abort(message = mess)
  }

  # map over waypoints
  cellWeights <- purrr::map_dfc(endpoints, function(endpoint) {
    # We find the path between the endpoint and the root
    path <- igraph::shortest_paths(gr, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    # We find the cells that map along that path
    cells_in_path <- y_to_cells %>%
      dplyr::filter(Y %in% path) %>%
      dplyr::pull(cells)
    tibble::tibble({{ endpoint }} := as.numeric(trajectory$cell_ids %in% cells_in_path))
  }) %>% magrittr::set_rownames(trajectory$cell_ids) # warning messages : setting rownames on a tibble


  pseudotime <- matrix(trajectory$pseudotime,
    ncol = ncol(cellWeights),
    nrow = length(trajectory$cell_ids), byrow = FALSE
  )

  return(list(cellWeights = cellWeights, pseudotime = pseudotime))
}

# genes <- trajectory$expression %>% Matrix::colSums() %>% tibble::enframe() %>% dplyr::filter(value > 10) %>% dplyr::pull(name)
# trajectory$expression[,genes]
# 
# test <- dynverse_to_tradeseq(trajectory)
