wrapper_scorpius <- function(traj_space)  {
  output <- dynwrap::wrap_data(cell_ids = names(traj_space[['traj']][['time']])) %>% 
    dynwrap::add_linear_trajectory(pseudotime = traj_space[['traj']][['time']])
  dimred_segment_points <- traj_space[['traj']][['path']]
  dimred_segment_progressions <- output$progressions %>% dplyr::select("from", "to", "percentage")
  output <- output %>% dynwrap::add_dimred(dimred = traj_space[['space']],
                                           dimred_segment_points = dimred_segment_points,
                                           dimred_segment_progressions = dimred_segment_progressions,
                                           connect_segments = TRUE) %>%
    dynwrap::label_milestones(labelling = c("milestone_begin" = "start", "milestone_end" = "end" ) ) %>%
    dynwrap::add_root(root_milestone_id = "milestone_begin")
  output
}

get_space_scorpius <- function(params, path, filename){
  if(length(params) > 1){
    log_counts <- fs::path(path,paste0(filename,"_log_counts.rds")) %>% readRDS()
    space <- SCORPIUS::reduce_dimensionality(log_counts, dist = params$dist, ndim = params$ndim, num_landmarks = params$num_landmarks) 
  } else {
    space <- fs::path(path,paste0(filename,"_redim_",params,".rds")) %>% readRDS()
  }
  return(space)
}

run_scorpius <- function(path,filename){
  path %<>% fs::path("scorpius")
  # from https://github.com/rcannood/SCORPIUS/blob/devel/vignettes/seurat.md 
  scorpius_params <- yaml::read_yaml(fs::path(path,"scorpius_params.yaml")) %>%
    filter_computed_traj(path)
  
  create_trajectory_dir(scorpius_params,path)
  
  groups <- fs::path(path,paste0(filename,"_grouping.rds")) %>% readRDS()
  log_counts <- fs::path(path,paste0(filename,"_log_counts.rds")) %>% readRDS()
  raw_counts <- fs::path(path,paste0(filename,"_raw_counts.rds")) %>% readRDS()
  
  output <- purrr::map(scorpius_params, 'space') %>%
    purrr::map(get_space_scorpius, path = path, filename = filename) %>%
    purrr::map2(scorpius_params, ~ list("traj" = SCORPIUS::infer_trajectory(.x, 
                                                              k = .y$k,
                                                              thresh = .y$thresh,
                                                              maxit = .y$maxit,
                                                              stretch = .y$stretch,
                                                              smoother = .y$smoother,
                                                              approx_points = .y$approx_points),
                                        "space" = .x)
                ) %>%
    purrr::map(wrapper_scorpius) %>%
    purrr::map(dynwrap::add_expression, counts = raw_counts, expression = log_counts) %>%
    purrr::map(dynwrap::add_grouping, grouping = groups)
  return(output)
}

data_prep_scorpius <- function(sobj,filename,reductions,path, cell_identities){
  fs::path(path,"scorpius") %>% fs::dir_create()
  Matrix::t(sobj[["RNA"]]@counts) %>% saveRDS(file = fs::path(path,"scorpius",paste0(filename,"_raw_counts.rds")))
  Matrix::t(sobj[["RNA"]]@data) %>% saveRDS(file = fs::path(path,"scorpius",paste0(filename,"_log_counts.rds")))
  available.reductions <- purrr::map_chr(reductions,~stringr::str_extract(names(sobj@reductions), .x) %>% purrr::discard(is.na))
  purrr::walk(available.reductions, ~ sobj[[.]]@cell.embeddings %>% saveRDS(file = fs::path(path,"scorpius",paste0(filename,"_redim_",.x,".rds"))))
  
  # Create params files
  scorpius_default <- list('space' = list('dist' = "spearman", "ndim" = 2L, "num_landmarks" = 1000),
                           'k' = 4L,
                           'thresh' = 0.001,
                           'maxit' = 10L, 
                           "stretch" = 0,
                           "smoother" = "smooth_spline",
                           "approx_points" = 100L)
  scorpius_params <- purrr::map(available.reductions, function(x) {scorpius_default['space'] <- x
  return(scorpius_default)
  } ) %>% purrr::set_names(purrr::map_chr(.,"space") %>% stringr::str_c("seurat_",.))
  scorpius_params[['default']] <- scorpius_default
  scorpius_default[['k']] <- nlevels(sobj@meta.data[[cell_identities]])
  scorpius_params[['default_ncluster']] <- scorpius_default
  yaml::write_yaml(scorpius_params, file = fs::path(path,"scorpius","scorpius_params.yaml"))
  #grouping
  sobj[[cell_identities]] %>% 
    tibble::as_tibble(rownames = "cell_id") %>% 
    dplyr::rename(group_id = {{cell_identities}}) %>%
    saveRDS(file = fs::path(path,"scorpius",paste0(filename,"_grouping.rds")))
}
