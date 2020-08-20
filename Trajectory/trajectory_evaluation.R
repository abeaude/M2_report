eval_trajectory <- function(trajectories, metrics_id = c("correlation","rf_nmse","lm_rsq","featureimp_cor","featureimp_wilcox")) {
  # add cell waypooints (required for trajectory evaluation)
  traj_to_eval <- purrr::map(trajectories, ~ purrr::map(., ~ purrr::map(.x, dynwrap::add_cell_waypoints))) %>%
    # prepare for evaluation
    prepare_eval()

  pb <- progress::progress_bar$new(format = "Trajectory evaluation : |:bar| :percent (:elapsedfull)", total = length(traj_to_eval), complete = "█", incomplete = " ", current = " ", clear = FALSE)
  pb$tick(0)
  res <- purrr::map_dfr(traj_to_eval, metrics_trajectory, metrics_id = metrics_id, pb = pb)
  pb$terminate()

  return(res)
}

# Compute metrics to assess the quality of the trajectory based on the metrics available
# in dyneval package
# Input : a list of 2 elements (the trajectory)
metrics_trajectory <- function(list_traj, metrics_id = c("correlation","rf_nmse","lm_rsq","featureimp_cor","featureimp_wilcox"), pb = NULL) {
  if (is.null(metrics_id)) {
    metrics_id <- dyneval::metrics %>% dplyr::pull(metric_id)
  } else if (any(!metrics_id %in% dyneval::metrics$metric_id)) {
    metrics_id <- metrics_id[metrics_id %in% dyneval::metrics$metric_id]
    unkown_metrics <- metrics_id[!metrics_id %in% dyneval::metrics$metric_id]
    warning("Unknown metrics : ", stringr::str_c(unkown_metrics, collapse = ", "))
  }

  split_1 <- names(list_traj[1]) %>%
    stringr::str_split("_-_") %>%
    unlist()
  split_2 <- names(list_traj[2]) %>%
    stringr::str_split("_-_") %>%
    unlist()
  metrics <- dyneval::calculate_metrics(dataset = list_traj[[1]], model = list_traj[[2]], metrics = metrics_id) %>%
    dplyr::select(!!!metrics_id) %>%
    dplyr::mutate(subset_1 = split_1[1], subset_2 = split_2[1], method_1 = split_1[2], method_2 = split_2[2], params_1 = split_1[3], params_2 = split_2[3])
  if (!is.null(pb)) {
    pb$tick(1)
  }
  return(metrics)
}

plot_trajectory_metrics <- function(metrics_res, metrics_id = c("correlation","rf_nmse","lm_rsq","featureimp_cor","featureimp_wilcox")) {
  eval_plot <- dplyr::select(metrics_res, method_1, subset_1) %>%
    dplyr::distinct() %>%
    dplyr::group_by(subset_1) %>%
    dplyr::summarise(method = list(method_1), .groups = "drop") %>%
    purrr::map(as.list) %>%
    purrr::pmap(~ purrr::cross2(..1, ..2)) %>%
    rlang::flatten() %>%
    purrr::set_names(nm = purrr::map(., ~ unlist(.) %>% stringr::str_c(collapse = "_-_"))) %>%
    purrr::map(~ tidyr::pivot_longer(metrics_res, !!metrics_id, names_to = "metric_id") %>%
      dplyr::filter(subset_1 == .x[[1]], subset_2 == .x[[1]]) %>%
      dplyr::inner_join(dyneval::metrics %>% dplyr::select(metric_id, plotmath, perfect, worst), by = "metric_id") %>%
      dplyr::filter(method_1 == .x[[2]]) %>%
      dplyr::mutate(label = dplyr::case_when(
        abs(value - worst) < 0.1 ~ "bad",
        abs(value - perfect) <= 0.1 ~ "good"
      ),
      group = 1:nrow(.)) %>%
      ggplot2::ggplot(ggplot2::aes(y = paste(params_1, " vs ", params_2), x = plotmath, color = value)) +
      ggplot2::scale_x_discrete(label = label_parse) +
      ggplot2::geom_point(size = 5) +
      ggforce::geom_mark_rect(ggplot2::aes(group = group, filter = label == "bad"), show.legend = TRUE, expand = grid::unit(3, "mm"), radius = grid::unit(2, "mm"), color = "red") +
      ggforce::geom_mark_rect(ggplot2::aes(group = group, filter = label == "good"), show.legend = TRUE, expand = grid::unit(3, "mm"), radius = grid::unit(2, "mm"), color = "green") +
        ggplot2::facet_grid(~method_2) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank(), plot.caption = ggtext::element_markdown(size = 10, halign = 0), panel.background = ggplot2::element_blank()) +
      ggplot2::labs(caption = paste("<span style='color:#ff0000;'>\u274c Bad score</span>", "<span style='color:#00ff00;'>\u2714 Good score</span>", sep = "<br>"), y = NULL, x = "Metrics", title = stringr::str_to_title(.x[[2]]) %>% stringr::str_c(" vs"), color = "Metric score", subtitle = stringr::str_c("Subset : ", .x[[1]])))

  return(eval_plot)
}

label_parse <- function(breaks) {
  parse(text = breaks)
}

save_eval_plot <- function(eval_plot, path) {
  
  assertthat::assert_that(rlang::is_named(eval_plot))
    purrr::iwalk(eval_plot, function(plot, name) {
      pbuild <- ggplot2::ggplot_build(plot)
      panel_number <- pbuild$layout$panel_params %>% length()
      x_range <- pbuild$layout$panel_params[[1]]$x.range %>% diff()
      y_range <- pbuild$layout$panel_params[[1]]$y.range %>% diff()
      n_char <- pbuild$layout$panel_params[[1]]$y$breaks %>%
        stringr::str_count() %>%
        max()
      
      # To ensure a good dimension of the plot and readability
      width <- 0.37 * x_range * panel_number + 0.08 * (panel_number - 1) + 1 + (n_char + 1) * 0.06
      height <- 0.28 * y_range + 0.5 + 1
      
      # Get where to save the plot based on list name 
      name_split <- name %>%
        stringr::str_split("_-_") %>%
        unlist()
      subset <- name_split[1]
      method <- name_split[2]
      path_plot <- fs::path(path, "TRAJECTORIES", subset, method)
      
      # saving
      ggplot2::ggsave(filename = "trajectory_evaluation.png", plot = plot, width = width, height = height, path = path_plot)
    })
}

load_all_traj <- function(path, subset = NULL, method = NULL) {
  # Load previous trajectory from the disk
  # Return a list with the same structure as the output of run_ti
  # USe of recursivity
  # Possibility to only return some subset/method
  
  assertthat::assert_that(assertthat::is.dir(path))
  assertthat::assert_that(assertthat::is.readable(path))

  if (fs::path_file(path) != "TRAJECTORIES") {
    path %<>% fs::path("TRAJECTORIES")
  }
  
  # Recursive function to create the different level of the list
  # .l a named list : the different level of the list are represented by / in the name
  # depth : the current level of depth in the list
  # max_depth : stopping condition, what is the depth maximum of the list (# words between /)
  create_depth_list <- function(.l, depth, max_depth) {
    if (depth == max_depth) {
      new_name <- .l %>%
        names() %>%
        stringr::str_split("/") %>%
        purrr::map(depth) %>%
        unlist()
      purrr::set_names(.l, new_name) %>%
        purrr::map(readRDS)
    } else {
      sub_name <- names(.l) %>%
        stringr::str_split("/") %>%
        purrr::map(depth) %>%
        unlist() %>%
        unique() %>%
        purrr::set_names(., nm = .)
      purrr::map(sub_name, ~ purrr::keep(.l, stringr::str_detect(names(.l), .))) %>%
        purrr::map(create_depth_list, depth = depth + 1, max_depth = max_depth)
    }
  }
  # Construct the list by first finding all the rds correponding to a trjaectory result (dynverse object)
  output <- fs::dir_ls(path, regexp = "trajectory_result.rds", recurse = TRUE) %>%
    purrr::set_names(purrr::map_chr(., ~ fs::path_rel(., start = path) %>% fs::path_dir())) %>%
    create_depth_list(1, names(.) %>% stringr::str_split("/") %>% purrr::map_int(length) %>% unique())

  filter_traj_list(output, subset = subset, method = method)
}

filter_traj_list <- function(traj, subset = NULL, method = NULL) {
  assertthat::assert_that(all(subset %in% names(traj)))
  
  if (is.null(subset) & is.null(method)) {
    message("Nothing to filter, returning the same object")
    return(traj)
  } else if (!is.null(subset) & is.null(method)) {
    return(traj[subset])
  } else if (is.null(subset) & !is.null(method)) {
    return(purrr::map(traj, ~ magrittr::extract(., method) %>% purrr::compact()) %>%
      purrr::compact())
  } else if (!is.null(subset) & !is.null(method)) {
    return(traj[subset] %>%
      purrr::map(~ magrittr::extract(., method) %>% purrr::compact()) %>%
      purrr::compact())
  }
}

recursive_names <- function(x, parent_name = NULL) {
  if (rlang::is_bare_list(x)) {
    if (!is.null(parent_name)) {
      x_names <- paste(parent_name, names(x), sep = "_-_")
    } else {
      x_names <- names(x)
    }

    purrr::map2(
      x,
      x_names,
      recursive_names
    ) %>%
      purrr::reduce(c)
  } else {
    return(parent_name)
  }
}

squashed_name <- function(x) {
  purrr::set_names(rlang::squash(x), recursive_names(x))
}

prepare_eval <- function(trajectories) {
  # Function to prepare the trajectory for evaluation
  # a list with sub list of size 2 to make 2 by 2 comparisons of the trajectory
  traj_squashed <- squashed_name(trajectories) # remember the methods of trajectory inference
  purrr::imap(traj_squashed, function(x, y, l) {
    purrr::imap(l, ~ list(.x, x) %>% purrr::set_names(c(.y, y))) %>% unname()
  }, l = traj_squashed) %>%
    purrr::flatten() %>%
    purrr::map(~ purrr::discard(., names(.) %>% unique() %>% length() == 1 %>% rep(2))) %>% # Remove entry with 2 times the same trajectory
    purrr::compact() %>%
    purrr::keep(~ all(.[[1]]$cell_ids %in% .[[2]]$cell_ids) & all(.[[2]]$cell_ids %in% .[[1]]$cell_ids)) # Keep only the trajectory where cell_ids match (avoid error when computing metrics)
}
