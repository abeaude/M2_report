# if (any(names(a) == "cells")){
#   sobj <- subset(sobj, cells = a[['cells']])
# } else if (any(names(a) == "idents")){
#   sobj <- subset(sobj, idents = a[["idents"]])
# } else {
#   # NOT SAFE
#   # conditions <- purrr::imap_chr(a, ~ stringr::str_c(.y ," %in% c(", stringr::str_c(.x, collapse = ","),")" )) %>% stringr::str_c(collapse = " & ")
#   # subset_cmd <- stringr::str_c("subset(sobj,subset = ", conditions, ")")
#   # sobj <- eval(parse(text = subset_cmd))
#   # SAFE
#   sobj <- subset_sobj(sobj = sobj, conditions = a)
#
# }

# subset a Seurat object based on metadata
subset_sobj <- function(sobj, conditions) {
  assertthat::assert_that(is_seurat(sobj))
  assertthat::assert_that(rlang::is_bare_list(conditions))
  assertthat::assert_that(is_named(conditions))

  metadata <- sobj@meta.data %>%
    tibble::as_tibble(rownames = "cells")

  if (any(!names(conditions) %in% colnames(metadata))) {
    names_cond <- names(conditions) %in% colnames(metadata)
    names_not_in_list <- names(conditions)[!names_cond]
    conditions <- purrr::keep(conditions, names_cond)
    warning(stringr::str_c(names_not_in_list, collapse = ", "), " are not in the seurat object metadata, they will be ignored")
  }

  cells <- purrr::imap(conditions, function(x, y) {
    dplyr::filter(metadata, .data[[y]] %in% x) %>% dplyr::pull(cells)
  }) %>%
    purrr::reduce(intersect)
  if (length(cells) == 0) {
    stop("No cells were found")
  } else {
    sobj[, cells]
  }
}

data_preparation_trajectories <- function(sobj = NULL, method = c("STREAM", "slingshot", "PAGA", "scorpius", "monocle3"), subset = NULL, cell_groups = "seurat_clusters", path = NULL, dimred = NULL, start_cell = Seurat::Cells(sobj)[1], seed = NULL, filename = "data_for_trajectories") {
  assertthat::assert_that(is_seurat(sobj))
  assertthat::assert_that(is.character(method))
  assertthat::assert_that(has_metadata(sobj, cell_groups))

  if (is.null(path)) {
    stop("Please specify a path to save all the datasets for trajectories")
  }
  assertthat::assert_that(assertthat::is.dir(path))
  assertthat::assert_that(assertthat::is.writeable(path))
  assertthat::assert_that(has_cell(sobj, start_cell))
  assertthat::assert_that(is.character(filename))

  seed <- seed %||% sobj@misc$params$seed %||% 1337L

  if (is.null(dimred)) {
    dimred <- stringr::str_c("(?i)^", names(sobj@reductions) , "$")
  } else {
    assertthat::assert_that(has_dimred(sobj, dimred))
    dimred <- stringr::str_c("(?i)^", dimred, "$")
  }

  if (any(stringr::str_detect(dimred, "(?i)UMAP"))) {
    warning("You will run trajecotry inference on UMAP, this is not recommended : see ",
      "https://www.biorxiv.org/content/10.1101/689851v3",
      call. = FALSE
    )
  }

  assay <- Seurat::DefaultAssay(sobj)
  if (sum(dim(sobj@assays[[assay]]@scale.data)) == 0) {
    sobj <- Seurat::ScaleData(object = sobj, assay = assay)
  }
  
  if (is.null(subset)) {
    subdir_name <- "all_object"
  } else if (any(names(subset) == "cells")) {
    subdir_name <- name_cell_directory(fs::path(path,"TRAJECTORIES"))
    sobj_sub <- subset(sobj, cells = subset[["cells"]])
  } else if (any(names(subset) == "idents")) {
    subdir_name <- stringr::str_c("idents", stringr::str_c(subset[["idents"]], collapse = "-"), sep = "_")
    sobj_sub <- subset(sobj, idents = subset[["idents"]])
  } else {
    subdir_name <- purrr::imap_chr(subset, ~ stringr::str_c(.y, stringr::str_c(.x, collapse = "-"), sep = "_")) %>% stringr::str_c(collapse = "-_-")
    sobj_sub <- subset_sobj(sobj = sobj, conditions = subset)
  }
  # Experimental
  # display a message
  ## subsetting the graph
  if(!is.null(subset)){
    cells <- colnames(sobj_sub)
    sobj_sub@graphs <- purrr::map(sobj@graphs, ~ .[cells, cells])
    sobj <- sobj_sub
    rm(sobj_sub)
  }
  
  # DROP unused factor when subsetting
  sobj@meta.data %<>%
    tibble::rownames_to_column("cells") %>%
    dplyr::mutate(dplyr::across(where(is.factor), forcats::fct_drop)) %>%
    tibble::column_to_rownames("cells")

  # Create directories
  fs::path(path, "TRAJECTORIES", subdir_name) %>%
    fs::dir_exists() %>%
    tibble::enframe() %>%
    dplyr::filter(!value) %>%
    dplyr::pull(name) %>%
    fs::dir_create()

  # if cells subsetting write cells inside directory
  if (any(names(subset) == "cells")) {
    subset[["cells"]] %>% readr::write_lines(fs::path(path, "TRAJECTORIES", subdir_name, "cells", ext = "txt"))
  }

  # check function
  funs <- stringr::str_c("data_prep_", stringr::str_to_lower(method)) %>% check_functions()

  args <- funs %>%
    purrr::set_names(.) %>%
    purrr::map(get, envir = globalenv()) %>%
    purrr::map(rlang::fn_fmls) %>%
    purrr::map(names) %>%
    purrr::map(function(l, sobj, path, filename, redim, cell_identities) {
      out <- list(
        "sobj" = sobj,
        "path" = fs::path(path, "TRAJECTORIES", subdir_name),
        "filename" = filename,
        "cell_identities" = cell_identities
      )
      if (any(l == "reductions")) {
        out[["reduction"]] <- dimred
      }
      if (any(l == "seed")) {
        out[["seed"]] <- seed
      }
      if (any(l == "start_cell")) {
        out[["start_cell"]] <- start_cell
      }
      return(out)
    }, sobj = sobj, path = path, filename = filename, redim = redim, cell_identities = cell_groups)

  funs %>%
    purrr::walk2(args, ~ rlang::exec(.x, !!!.y))
}

run_ti <- function(path = NULL, filename = "data_for_trajectories", force = FALSE) {
  assertthat::assert_that(assertthat::is.dir(path))
  assertthat::assert_that(assertthat::is.writeable(path))
  assertthat::assert_that(is.character(filename))

  silent_exec <- purrr::quietly(rlang::exec)
  subsets <- fs::dir_ls(fs::path(path, "TRAJECTORIES"), type = "directory") %>% fs::path_file()

  if (force) {
    fs::dir_ls(fs::path(path, "TRAJECTORIES"), type = "directory") %>% # subsets
      purrr::map(~ fs::dir_ls(path = .x, type = "directory")) %>% # method
      purrr::map(~ fs::dir_ls(path = .x, type = "directory")) %>% # params
      unlist() %>%
      purrr::map(~ fs::dir_delete(.x)) # delete previous analysis
  }

  output <- list()
  n_method <- fs::dir_ls(fs::path(path, "TRAJECTORIES"), type = "directory") %>%
    purrr::map(~ fs::dir_ls(path = .x, type = "directory")) %>%
    unlist() %>%
    length()
  pb_ti <- progress::progress_bar$new(format = "TI :subset - :method : |:bar| :percent (:elapsedfull)", total = n_method, complete = "█", incomplete = " ", current = " ", clear = FALSE)

  for (subset in subsets) {
    method <- fs::path(path, "TRAJECTORIES", subset) %>%
      fs::dir_ls(type = "directory") %>%
      fs::path_file()

    args <- list("filename" = filename, "path" = fs::path(path, "TRAJECTORIES", subset))
    trajectory <- stringr::str_c("run_", stringr::str_to_lower(method)) %>%
      check_functions() %>%
      purrr::set_names(stringr::str_extract(., "_.*") %>% stringr::str_remove("_")) %>%
      purrr::imap(function(fun, method.) {
        res <- silent_exec(fun, !!!args)
        pb_ti$tick(tokens = list(subset = subset, method = method.))
        res
      })
    output[[subset]] <- purrr::map(trajectory, "result") %>% purrr::compact()
  }
  pb_ti$terminate()
  output %<>% purrr::compact()

  n_traj <- get_number_trajecotry(output)

  pb <- progress::progress_bar$new(format = "Adding Top Features : |:bar| :percent (:elapsedfull)", total = n_traj, complete = "█", incomplete = " ", current = " ", clear = FALSE)
  pb_time <- progress::progress_bar$new(format = "Computing pseudotime : |:bar| :percent (:elapsedfull)", total = n_traj, complete = "█", incomplete = " ", current = " ", clear = FALSE)
  pb_genes <- progress::progress_bar$new(format = "Finding DE genes along trajectory : |:bar| :percent (:elapsedfull)", total = n_traj, complete = "█", incomplete = " ", current = " ", clear = FALSE)

  output %<>%
    add_feature_oi(n_features = 20, pb = pb) %>%
    add_pseudotime_traj(pb = pb_time) %>%
    add_gene_importances(pb = pb_genes)

  pb$terminate()
  pb_time$terminate()
  pb_genes$terminate()
  return(output)
}

check_functions <- function(fun_names) {
  assertthat::assert_that(is.character(fun_names))

  fun_exists <- fun_names %>%
    purrr::set_names(.) %>%
    purrr::map_lgl(exists)

  if (any(!fun_exists)) {
    not_def <- fun_exists[!fun_exists] %>% names()
    fun_prefix <- Biobase::lcPrefix(fun_names)
    methods <- stringr::str_split(not_def, "_") %>%
      purrr::map_chr(~ rev(.)[1])
    stringr::str_c("No functions :\n", stringr::str_c("- ", not_def, collapse = "\n"), "\nTo enable support for methods : ", stringr::str_c(methods, collapse = ", "), ", consider adding `", fun_prefix, "` functions") %>% warning(call. = FALSE)
  }
  return(fun_exists[fun_exists] %>% names())
}

filter_computed_traj <- function(params, path) {
  assertthat::assert_that(assertthat::is.dir(path))
  assertthat::assert_that(assertthat::is.readable(path))
  assertthat::assert_that(is_named(params))

  not_done <- params %>%
    names() %>%
    fs::path(path, .) %>%
    fs::dir_exists() %>%
    tibble::enframe(name = "Path", value = "exists") %>%
    dplyr::mutate(Name = fs::path_file(Path)) %>%
    dplyr::filter(!exists) %>%
    dplyr::pull(Name)

  return(magrittr::extract(params, not_done))
}

create_trajectory_dir <- function(params, path) {
  assertthat::assert_that(assertthat::is.dir(path))
  assertthat::assert_that(assertthat::is.writeable(path))
  assertthat::assert_that(is_named(params))

  params %>%
    names() %>%
    fs::path(path, .) %>%
    fs::dir_create()
}

save_trajectories <- function(trajectories, path) {
  assertthat::assert_that(assertthat::is.dir(path))
  assertthat::assert_that(assertthat::is.writeable(path))
  assertthat::assert_that(list_depth(trajectories) >= 3)

  trajectories <- squashed_name(trajectories)

  purrr::iwalk(trajectories, function(x, nm) {
    rel_path <- stringr::str_split(nm, "_-_") %>%
      unlist() %>%
      stringr::str_c(collapse = "/")
    if (ggplot2::is.ggplot(x)) {
      ggplot2::ggsave(plot = x, filename = fs::path(path, "TRAJECTORIES", paste0(rel_path, "_plot"), ext = "png"), width = 7, height = 7)
    } else {
      saveRDS(x, fs::path(path, "TRAJECTORIES", rel_path, "trajectory_result.rds"))
    }
  })
}

plot_trajectories <- function(trajectories, pb = if (get_number_trajecotry(trajectories) > 1) progress::progress_bar$new(format = "Plotting : |:bar| :percent (:elapsedfull)", total = get_number_trajecotry(trajectories), complete = "█", incomplete = " ", current = " ", clear = FALSE) else NULL, extra_grouping = NULL) {
  if (list_depth(trajectories) >= 1) {
    purrr::map(trajectories, plot_trajectories, pb = pb, extra_grouping = extra_grouping)
  } else if (list_depth(trajectories) == 0 & dynwrap::is_data_wrapper(trajectories)) {
    dim_test <- purrr::possibly(dynplot::plot_dimred, otherwise = NULL)(trajectories, color_cells = "pseudotime", plot_trajectory = TRUE, label_milestones = TRUE, color_trajectory = "none")

    plot_traj <- ifelse(is.null(dim_test), FALSE, TRUE)

    dim_pseudotime <- dynplot::plot_dimred(trajectories, color_cells = "pseudotime", plot_trajectory = plot_traj, label_milestones = TRUE, color_trajectory = "none")
    dim_grouping <- dynplot::plot_dimred(trajectories, color_cells = "grouping", plot_trajectory = plot_traj, label_milestones = TRUE, color_trajectory = "none")

    pbuild <- ggplot2::ggplot_build(dim_pseudotime)
    x_range <- pbuild$layout$panel_params[[1]]$x.range %>% diff()
    y_range <- pbuild$layout$panel_params[[1]]$y.range %>% diff()

    if (x_range > y_range) { # align horizontaly : /
      dimred_plot <- dim_pseudotime / dim_grouping
      dimred_plot[[1]] <- dimred_plot[[1]] + ggplot2::theme(legend.position = "right", legend.direction = "vertical")
    } else { # align verticaly : |
      dimred_plot <- dim_pseudotime | dim_grouping
      dimred_plot[[1]] <- dimred_plot[[1]] + ggplot2::theme(legend.position = "left", legend.direction = "vertical")
    }

    ## Extra grouping plot :
    if (!is.null(extra_grouping)) {
      assertthat::assert_that(is.list(extra_grouping), is_named(extra_grouping))

      extra_grouping_plot <- purrr::map(
        extra_grouping,
        ~ dynplot::plot_dimred(trajectories, color_cells = "grouping", plot_trajectory = plot_traj, label_milestones = TRUE, color_trajectory = "none", grouping = .x)
      ) %>%
        purrr::set_names(nm = paste0("dimred_", names(.), "_grouping"))
    } else {
      extra_grouping_plot <- NULL
    }

    dendro_plot <- purrr::possibly(dynplot::plot_dendro, NULL)(trajectories, color_cells = "pseudotime") /
      purrr::possibly(dynplot::plot_dendro, NULL)(trajectories, color_cells = "grouping")
    if (ggplot2::is.ggplot(dendro_plot)) {
      dendro_plot[[1]] <- dendro_plot[[1]] + ggplot2::theme(legend.position = "right", legend.direction = "vertical")
    } else {
      dendro_plot <- NULL
    }


    # graph_plot <- dynplot::plot_graph(wrap_traj, color_cells = "pseudotime", label_milestones = TRUE, plot_milestones = TRUE) /
    # dynplot::plot_graph(wrap_traj, color_cells = "grouping", label_milestones = TRUE, plot_milestones = TRUE)
    #
    # graph_plot[[1]] <- graph_plot[[1]] + ggplot2::theme(legend.position = "right",legend.direction = 'vertical')

    topo_plot <- dynplot::plot_topology(trajectories)

    onedim_plot <- dynplot::plot_onedim(trajectories, color_cells = "pseudotime", label_milestones = TRUE) /
      dynplot::plot_onedim(trajectories, color_cells = "grouping", label_milestones = TRUE, orientation = -1)

    onedim_plot[[1]] <- onedim_plot[[1]] + ggplot2::theme(legend.position = "right", legend.direction = "vertical")

    heatmap_plot <- dynplot::plot_heatmap(trajectories, color_cells = "grouping", label_milestones = TRUE, features_oi = trajectories[["features_oi"]])

    output <- list(
      "dimred" = dimred_plot,
      "dendro" = dendro_plot,
      "top" = topo_plot,
      "onedim" = onedim_plot,
      "heatmap" = heatmap_plot
    ) %>%
      append(extra_grouping_plot) %>%
      purrr::compact()

    if (rlang::has_name(trajectories, "gene_importances")) {
      gene_sel <- trajectories[["gene_importances"]]
      expr_sel <- trajectories$expression[, gene_sel]
      expression <- as.matrix(expr_sel)
      time <- trajectories$pseudotime

      # label
      dynwrap::get_milestone_labelling(trajectories)
      start <- trajectories[["root_milestone_id"]]
      end <- trajectories[["pseudotime"]] %>%
        {
          names(.)[which.max(.)]
        } %>%
        {
          dplyr::filter(trajectories[["milestone_percentages"]], cell_id == .) %>%
            dplyr::arrange(desc(percentage)) %>%
            dplyr::slice(1) %>%
            dplyr::pull(milestone_id)
        }

      label <- trajectories$milestone_ids %>%
        purrr::set_names(trajectories[["milestone_labelling"]]) %>%
        .[which(. %in% c(start, end))]
      label <- {
        if (rlang::is_named(label)) names(label) else label
      }

      heatmap_traj <- heatmap_trajectory(expression, time, label = label) +
        patchwork::plot_annotation(title = "Trajectory heatmap")
      output[["heatmap_trajectory"]] <- heatmap_traj
    }

    if (!is.null(pb)) pb$tick(1)

    return(output)
  } else {
    rlang::abort("Change ME")
  }
}

add_feature_oi <- function(trajectories, n_features = 20, pb = NULL) {
  assertthat::assert_that(is.numeric(n_features))

  if (list_depth(trajectories) >= 1) {
    purrr::map(trajectories, add_feature_oi, n_features = n_features, pb = pb)
  } else if (list_depth(trajectories) == 0 & dynwrap::is_data_wrapper(trajectories)) {
    expression <- dynwrap::get_expression(trajectories, "expression")

    features_oi <- dynfeature::calculate_overall_feature_importance(trajectories, expression = expression) %>%
      dplyr::top_n(n_features, importance) %>%
      dplyr::pull(feature_id) %>%
      as.character()

    trajectories[["features_oi"]] <- features_oi
    if (!is.null(pb)) pb$tick(1)

    return(trajectories)
  } else {
    rlang::abort("trajectories is not an accepted object")
  }
}

list_depth <- function(x) {
  if (rlang::is_null(x) | ggplot2::is.ggplot(x) | dynwrap::is_data_wrapper(x)) {
    0L
  } else if (rlang::is_atomic(x)) {
    1L
  } else if (rlang::is_list(x)) {
    depths <- purrr::map_int(x, list_depth)
    1L + max(depths, 0L)
  } else {
    rlang::abort("`x` must be a vector")
  }
}

add_pseudotime_traj <- function(trajectories, pb = NULL) {
  if (list_depth(trajectories) >= 1) {
    purrr::map(trajectories, add_pseudotime_traj, pb = pb)
  } else if (list_depth(trajectories) == 0 & dynwrap::is_data_wrapper(trajectories)) {
    if (!is.null(pb)) pb$tick(1)
    dynwrap::add_pseudotime(trajectories)
  } else {
    rlang::abort("trajectories is not an accepted object")
  }
}

# Seurat::Reductions(sobj,slot = "umap")@cell.embeddings
project_traj <- function(trajectories, space, space_name) {
  assertthat::assert_that(is.matrix(space))
  assertthat::assert_that(is.character(space_name))

  if (list_depth(trajectories) >= 1) {
    purrr::map(trajectories, project_traj, space = space, space_name = space_name) %>%
      purrr::compact()
  } else if (list_depth(trajectories) == 0 & dynwrap::is_data_wrapper(trajectories)) {
    # subset the space to only contains cells in dynverse object
    space <- space[trajectories$cell_ids, ]
    if (!identical(unname(space), unname(trajectories$dimred))) {
      plot_name <- stringr::str_c("projected_on_", space_name)
      dimred_projected <- dynwrap::project_trajectory(trajectories, space)
      trajectories %>%
        dynwrap::add_dimred(
          dimred = space,
          dimred_milestones = dimred_projected$dimred_milestones,
          dimred_segment_progressions = dimred_projected$dimred_segment_progressions,
          dimred_segment_points = dimred_projected$dimred_segment_points
        ) %>%
        dynplot::plot_dimred(color_cells = "grouping", plot_trajectory = TRUE, label_milestones = TRUE, color_trajectory = "none") %>%
        list(.) %>%
        magrittr::set_names(plot_name)
    } else {
      warning("Original and projection space are identical", call. = FALSE)
      NULL
    }
  } else {
    rlang::abort("trajectories is not an accepted object")
  }
}

get_number_trajecotry <- function(trajectories) {
  rlang::squash(trajectories) %>% length()
}

name_cell_directory <- function(path){
  fs::dir_ls(path, type = "directory") %>% 
    fs::path_file() %>% 
    stringr::str_subset("cells_subsetting") %>%
    c("cells_subsetting") %>%
    make.names(unique = TRUE) %>%
    tail(1)
}