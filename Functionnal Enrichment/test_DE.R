lock_gestion <- function(path, test, DE_type) {
  assertthat::assert_that(assertthat::is.dir(path), assertthat::is.writeable(path))
  assertthat::assert_that(is.character(test))
  assertthat::assert_that(is_in(DE_type, c("1vs1", "1vsAll", "SvsS", "conditions")))

  contrast_lock <- fs::path_dir(path) %>%
    fs::path("contrast.txt.lock") %>%
    as.character()
  lock_files <- fs::path(path, test, stringr::str_c(test, "_", DE_type, "_markers.txt.lock")) %>%
    fs::file_exists() %>%
    tibble::enframe() %>%
    tibble::add_row(name = contrast_lock, value = fs::file_exists(contrast_lock))

  dplyr::filter(lock_files, value == FALSE) %>%
    dplyr::pull(name) %>%
    fs::file_create()

  dplyr::filter(lock_files, value == TRUE) %>%
    dplyr::pull(name) %>%
    fs::file_delete()
}

run_seurat_test <- function(params, DE_type, output_path, sobj, batch, cluster_conditions, verbose = FALSE, pb = NULL) {
  new_sobj <- calc_new.cluster(sobj, params$first_group, params$second_group, levels(Seurat::Idents(sobj))) # renommage des clusters
  if (!params$test %in% c("LR", "negbinom", "poisson")) { # batch doesn't work with MAST test and categorical batch
    batch <- NULL
  }

  Seurat_test(DE_type, output_path, params$test, new_sobj, params$first_group_name, params$second_group_name, params$min.pct, batch, verbose, cluster_conditions)
  if (!is.null(pb)) pb$tick()

  values <- c(DE_type, params$test, params$first_group_name, params$second_group_name, params$min.pct, ifelse(is.null(batch), NA, batch)) %>%
    stringr::str_replace_na() %>%
    stringr::str_c(collapse = "\t")
  lock <- filelock::lock(path = fs::path_dir(output_path) %>% fs::path("contrast.txt.lock"), exclusive = TRUE, timeout = Inf)
  readr::write_lines(values, path = fs::path_dir(output_path) %>% fs::path("contrast.txt"), append = TRUE)
  filelock::unlock(lock)
}

run_limma_edge <- function(params, clusters, counts = NULL, DE_type, design_str, batch, cluster_conditions, output_path, pb = NULL) {
  # For 1 vs ALL rename second group to All
  if (DE_type == "1vsAll") {
    params$second_group_name <- "All"
  }

  # counts <- sobj[["RNA"]]@counts
  cells.1 <- counts[, clusters$cluster %in% params$first_group]
  cells.2 <- counts[, clusters$cluster %in% params$second_group]
  cells <- c(colnames(cells.1), colnames(cells.2))

  # Calcul % d'expression pour chaque gène de chaque groupe, vérif si ça passe le min.pct et retourne les gènes restants
  features <- calc_features_pct(list(cells.1, cells.2), params$min.pct)
  if (length(x = features) == 0) {
    message("No features pass min.pct threshold")
  }
  # Objet DGEList + facteur de normalisation (profondeur)
  counts <- edgeR::DGEList(counts[rownames(features), cells])
  counts <- edgeR::calcNormFactors(counts)

  # récupération des effets batch
  if (!is.null(batch_list)) {
    for (batch in batch_list) {
      assign(x = batch, value = sobj@meta.data[cells, batch]) # récupération des effets batch
    }
  }

  # design
  if (stringr::str_detect(design_str, "C_")) {
    C_ <- clusters$cluster # nécessaire pour pourvoir modifier le nom dans le design ("clusters[,1]", en tant que titre de colonne, ne peut pas être modifié)
    design <- model.matrix(as.formula(design_str))
    colnames(design) <- c(paste0("C_", levels(C_)), as.vector(batch))
  } else {
    # new_clusters <- factor(clusters[c(colnames(cells.1), colnames(cells.2)),], exclude=NULL)
    new_clusters <- clusters %>%
      dplyr::filter(cells %in% c(colnames(cells.1), colnames(cells.2))) %>%
      dplyr::arrange(match(cluster, c(params$first_group, params$second_group))) %>%
      dplyr::pull(cluster) %>%
      forcats::fct_drop()

    design <- model.matrix(as.formula(design_str))
    colnames(design) <- sub("new_clusters", "C_", colnames(design))
    colnames(design) <- c(colnames(design)[grep("C_", colnames(design))], as.vector(batch))
  }

  # Contrast
  contr.matrix <- calc_contr.matrix(params$first_group, params$second_group, design)
  
  calc_DE_genes <- switch(params$test,
    Limma_Voom = Limma_voom_test,
    Limma_Trend = Limma_trend_test,
    EdgeR_LRT = EdgeR_LRT_test,
    EdgeR_QL = EdgeR_QL_test
  )
  calc_DE_genes(DE_type, output_path, counts, clusters, design, contr.matrix, params$first_group_name, params$second_group_name, params$min.pct, cluster_conditions)

  if (!is.null(pb)) pb$tick()

  values <- c(DE_type, params$test, params$first_group_name, params$second_group_name, params$min.pct, ifelse(is.null(batch), NA, batch)) %>%
    stringr::str_replace_na() %>%
    stringr::str_c(collapse = "\t")
  lock <- filelock::lock(path = fs::path_dir(output_path) %>% fs::path("contrast.txt.lock"), exclusive = TRUE, timeout = Inf)
  readr::write_lines(values, path = fs::path_dir(output_path) %>% fs::path("contrast.txt"), append = TRUE)
  filelock::unlock(lock)
}

run_limma_edge_parallel <- function(params, clusters, counts, DE_type, design_str, cluster_conditions, batch = NULL, output_path, pb = NULL) {
  n <- length(params)
  globals <- stringr::str_subset(ls(envir = globalenv()), "sobj", negate = TRUE)
  future_opt <- furrr::future_options(globals = globals, packages = NULL, seed = FALSE, lazy = FALSE, scheduling = 1) # set which object should be exported as a global (not the sobj)
  furrr::future_walk(params, run_limma_edge, DE_type = DE_type, clusters = clusters, design_str = design_str, output_path = output_path, counts = counts, batch = NULL, pb = NULL, cluster_conditions = cluster_conditions, .options = future_opt)
  if (!is.null(pb)) pb$tick(n)
}

run_DE <- function(sobj = NULL, clusters = NULL, DE_type = "1vsAll", method = c("Limma_Voom", "Limma_Trend", "EdgeR_LRT", "EdgeR_QL"), testToUse = c("wilcox", "bimod", "t", "poisson", "negbinom", "LR", "MAST"), first_group = NULL, second_group = NULL, cluster_conditions = NULL, min.pct = 0.1, path = NULL, batch = NULL, .parallel = FALSE, verbose = TRUE) {
  # convert NA to NULL
  na_args <- rlang::current_env() %>%
    as.list() %>%
    purrr::map_lgl(~ any(is.na(.))) %>%
    purrr::discard(isFALSE) %>%
    names()

  if (!rlang::is_empty(na_args)) {
    purrr::map(na_args, ~ assign(., NULL, envir = rlang::env_parent()))
  }

  assertthat::assert_that(is_seurat(sobj))
  assertthat::assert_that(is_in(DE_type, c("1vs1", "1vsAll", "SvsS", "conditions")))

  path <- fs::path_tidy(path)
  assertthat::assert_that(assertthat::is.dir(path), assertthat::is.writeable(path))
  assertthat::assert_that(is.numeric(min.pct), min.pct >= 0, min.pct <= 1)

  if (verbose) message("Running differential expression : ", DE_type)

  # Setting Default assay to RNA
  Seurat::DefaultAssay(sobj) <- "RNA"
  # check if data are normalized in RNA assay
  is.norm <- sobj@commands %>%
    names() %>%
    stringr::str_detect("NormalizeData.RNA") %>%
    any()
  if (!is.norm) {
    sobj <- Seurat::NormalizeData(sobj, verbose = FALSE)
  }

  if (DE_type == "conditions") {
    sobj <- subset(sobj, idents = cluster_conditions)
  }

  if (is.null(clusters) && (DE_type == "conditions")) {
    Seurat::Idents(sobj) <- "conditions"
  } else if (rlang::is_character(clusters, 1)) {
    assertthat::assert_that(clusters %in% colnames(sobj@meta.data), msg = paste(clusters, "is not in the metadata of the Seurat object"))
    Seurat::Idents(sobj) <- clusters
  }
  clusters <- Seurat::Idents(sobj) %>% tibble::enframe(name = "cells", value = "cluster")
  # check if cluster is a factor
  assertthat::assert_that(is.factor(clusters$cluster))

  clusters_list <- levels(clusters$cluster)

  if ((is.null(first_group)) && (!is.null(second_group))) {
    message("\tfirst_group is empty. Switching second_group and first_group.")
    first_group <- second_group
    second_group <- NULL
  }
  if ((DE_type == "SvsS" | DE_type == "1vsAll") & length(first_group) == 1 & length(second_group) == 1) {
    message("\tOnly one cluster in first_group and second_group, so switching DE_type to '1vs1'.")
    DE_type <- "1vs1"
  }
  if (DE_type == "1vsAll" & !is.null(second_group)) {
    warning("The second_group argument must be empty for comparisons_type to '1vsAll'. \n Set DE_type to '1vs1' or 'SvsS' to use it.")
  }
  if ((DE_type == "SvsS") & (!is.null(second_group)) & (length(intersect(as.character(first_group), as.character(second_group))) != 0)) {
    stop("first_group and second_group can not overlap.")
  }
  if ((DE_type == "SvsS") && (is.null(first_group))) {
    stop("In DE_type 'SvsS', first_group must not be empty.")
  }
  if ((DE_type == "1vs1" | DE_type == "1vsAll") & is.null(first_group)) {
    message("\tWill use all the known cluster for the first group")
    first_group <- as.character(clusters_list)
  }
  if (DE_type == "1vs1" & is.null(second_group)) {
    message("\tWill use all the known cluster for the second group")
    second_group <- as.character(clusters_list)
  }
  if (DE_type == "1vsAll") {
    second_group <- clusters_list %>%
      as.character() %>%
      list() %>%
      rep(length(first_group)) %>%
      purrr::map2(as.character(first_group), dplyr::setdiff)
  }
  if ((DE_type == "1vs1" | DE_type == "1vsAll") & !is.null(first_group)) {
    first_group <- as.character(first_group)
  }
  if (DE_type == "1vs1" & !is.null(second_group)) {
    second_group <- as.character(second_group)
  }
  if (DE_type == "SvsS" & is.null(second_group)) {
    second_group <- dplyr::setdiff(as.character(clusters_list), as.character(first_group))
  }
  if (DE_type == "SvsS") {
    first_group <- list(first_group)
    second_group <- list(second_group)
  }

  #### formation design_str ####
  if (DE_type == "1vsAll") {
    design_str <- "~0 + C_"
  } else if ((DE_type == "1vs1") | (DE_type == "SvsS") | (DE_type == "conditions")) {
    design_str <- "~0 + new_clusters"
  }
  if (!is.null(batch)) {
    batch <- factor(unlist(strsplit(batch, ",")))
    for (batches in batch) {
      design_str <- paste0(design_str, " + ", batches) # formation formule du design
    }
  }

  path.old <- path
  path <- fs::path(path, "DIFFERENTIAL_EXPRESSION", DE_type)
  # if(DE_type == "conditions") path <- fs::path(path, stringr::str_c("cluster_",stringr::str_c(cluster_conditions, collapse = "-")))

  # Create folders for Seurat function
  # Walk is equivalent to map but onlyu capture side effects
  if (verbose) message("\tCreating folders for outputs")
  create_folders_DE(path, test = c(testToUse, method))
  lock_gestion(path, c(testToUse, method), DE_type)
  on.exit(lock_gestion(path, c(testToUse, method), DE_type)) # when exiting (error or not) function run lock gestion to delete .lock file

  if (fs::path(path.old, "DIFFERENTIAL_EXPRESSION", "contrast.txt") %>% fs::file_exists()) {
    contrast_done <- readr::read_tsv(fs::path(path.old, "DIFFERENTIAL_EXPRESSION", "contrast.txt"),
      col_names = TRUE,
      col_types = readr::cols(
        DE_type = readr::col_character(),
        test = readr::col_character(),
        tested_cluster = readr::col_character(),
        control_cluster = readr::col_character(),
        min.pct = readr::col_double(),
        batch = readr::col_character()
      )
    ) %>%
      dplyr::filter(DE_type == !!DE_type) %>%
      dplyr::select(test, tested_cluster, control_cluster, min.pct)
  } else {
    headers <- stringr::str_c("DE_type", "test", "tested_cluster", "control_cluster", "min.pct", "batch", sep = "\t")
    readr::write_lines(headers, path = fs::path(path.old, "DIFFERENTIAL_EXPRESSION", "contrast.txt"), append = FALSE)
    contrast_done <- tibble::tibble(test = character(), tested_cluster = character(), control_cluster = character(), min.pct = numeric(), .rows = 0L)
  }


  #### Create grid of parameters for seurat ####
  if (length(testToUse) > 0) {
    seurat_params <- purrr::cross(list("test" = testToUse, "first_group" = first_group, "second_group" = second_group, "min.pct" = min.pct),
      .filter = function(test, first_group, second_group, min.pct) any(second_group == first_group)
    ) %>%
      purrr::map(~ purrr::list_modify(.,
        second_group_name = stringr::str_c(.$second_group, collapse = "_"),
        first_group_name = stringr::str_c(.$first_group, collapse = "_")
      )) %>%
      tibble::enframe() %>%
      tidyr::unnest_wider(col = "value") %>%
      dplyr::anti_join(contrast_done %>% dplyr::filter(test %in% c("wilcox", "bimod", "t", "poisson", "negbinom", "LR", "MAST", "roc", "DEseq2")), by = c("test", "first_group_name" = "tested_cluster", "second_group_name" = "control_cluster", "min.pct")) %>%
      tidyr::nest(value = -name) %>%
      dplyr::mutate(value = purrr::map(value, as.list)) %>%
      tibble::deframe() %>%
      unname() %>%
      purrr::map(purrr::flatten)
    if (verbose) {
      message("\tSeurat Tests")
      pb_seurat <- progress::progress_bar$new(format = "\t\tRunning : |:bar| :percent (:elapsedfull)", total = length(seurat_params), complete = "█", incomplete = " ", current = " ", clear = FALSE)
      pb_seurat$tick(0)
    } else {
      pb_seurat <- NULL
    }
    if (length(seurat_params)) purrr::walk(seurat_params, run_seurat_test, DE_type = DE_type, output_path = path, sobj = sobj, batch = batch, cluster_conditions = cluster_conditions, verbose = FALSE, pb = pb_seurat)
  }

  if (length(method) > 0) {
    limma_edgeR_params <- purrr::cross(list("test" = method, "first_group" = first_group, "second_group" = second_group, "min.pct" = min.pct),
      .filter = function(test, first_group, second_group, min.pct) any(second_group == first_group) # remove comp with itself
    ) %>%
      purrr::map(~ purrr::list_modify(.,
        second_group_name = stringr::str_c(.$second_group, collapse = "_"),
        first_group_name = stringr::str_c(.$first_group, collapse = "_")
      )) %>%
      tibble::enframe() %>%
      tidyr::unnest_wider(col = "value") %>%
      dplyr::anti_join(contrast_done %>% dplyr::filter(any(test %in% c("EdgeR_LRT", "EdgeR_QL", "Limma_Voom", "Limma_Trend"))), by = c("test", "first_group_name" = "tested_cluster", "second_group_name" = "control_cluster", "min.pct")) %>%
      tidyr::nest(value = -name) %>%
      dplyr::mutate(value = purrr::map(value, as.list)) %>%
      tibble::deframe() %>%
      unname() %>%
      purrr::map(purrr::flatten) %>%
      purrr::set_names(purrr::map_chr(., `[[`, "test") %>% make.names(unique = TRUE)) %>%
      .[stringr::str_sort(names(.))]

    counts <- sobj[["RNA"]]@counts
    rm(sobj)
    gc(verbose = FALSE)
    if (verbose) {
      message("\tLimma and/or EdgeR")
      pb_limma <- progress::progress_bar$new(format = "\t\tRunning : |:bar| :percent (:elapsedfull)", total = length(limma_edgeR_params), complete = "█", incomplete = " ", current = " ", clear = FALSE)
      pb_limma$tick(0)
    } else {
      pb_limma <- NULL
    }

    if (.parallel) {
      n <- length(limma_edgeR_params) / future::availableCores()["system"]
      limma_edgeR_params_parallel <- unname(split(limma_edgeR_params, cut(seq_along(limma_edgeR_params), n, labels = FALSE))) # cut in x intervals of size n
      purrr::walk(limma_edgeR_params_parallel, run_limma_edge_parallel, DE_type = DE_type, clusters = clusters, design_str = design_str, output_path = path, counts = counts, batch = batch, pb = pb_limma, cluster_conditions = cluster_conditions)
    } else {
      if (length(limma_edgeR_params)) purrr::walk(limma_edgeR_params, run_limma_edge, DE_type = DE_type, clusters = clusters, design_str = design_str, output_path = path, counts = counts, batch = batch, pb = pb_limma, cluster_conditions = cluster_conditions)
    }
  }
}

run_multiple_DE <- function(parameters = NULL, sobj, path, min.pct = 0.1, .parallel = FALSE, verbose = TRUE) {
  assertthat::assert_that(is_seurat(sobj))
  path <- fs::path_tidy(path)
  assertthat::assert_that(assertthat::is.dir(path), assertthat::is.writeable(path))
  assertthat::assert_that(is.numeric(min.pct), min.pct >= 0, min.pct <= 1)
  clusters <- "seurat_clusters"
  if (is.null(parameters)) {
    source("Functionnal Enrichment/shiny/app.R", local = TRUE)
    parameters <- runApp(multiple_DE_app, launch.browser = TRUE) %>%
      dplyr::mutate(
        first_group = purrr::map(first_group, ~ stringr::str_split(., ",") %>% unlist()),
        second_group = purrr::map(second_group, ~ stringr::str_split(., ",") %>% unlist()),
        method = purrr::map(method, ~ stringr::str_split(., ",") %>% unlist()),
        testToUse = purrr::map(testToUse, ~ stringr::str_split(., ",") %>% unlist()),
        cluster_conditions = purrr::map(cluster_conditions, ~ stringr::str_split(., ",") %>% unlist()),
        batch = purrr::map(batch, ~ stringr::str_split(., ",") %>% unlist())
      )
  }
  # save the table of parameters
  fs::path(path, "DIFFERENTIAL_EXPRESSION") %>% fs::dir_create()
  dplyr::mutate(
    parameters,
    first_group = purrr::map(first_group, ~ stringr::str_c(., collapse = ",")),
    second_group = purrr::map(second_group, ~ stringr::str_c(., collapse = ",")),
    method = purrr::map(method, ~ stringr::str_c(., collapse = ",")),
    testToUse = purrr::map(testToUse, ~ stringr::str_c(., collapse = ",")),
    cluster_conditions = purrr::map(cluster_conditions, ~ stringr::str_c(., collapse = ",")),
    batch = purrr::map(batch, ~ stringr::str_c(., collapse = ","))
  ) %>%
    tidyr::unnest(cols = c(first_group, second_group, method, testToUse, cluster_conditions, batch)) %>%
    readr::write_tsv(path = fs::path(path, "DIFFERENTIAL_EXPRESSION", "parameters", ext = "txt"))

  purrr::pwalk(parameters, run_DE, sobj = sobj, min.pct = min.pct, path = path, .parallel = .parallel, verbose = verbose)
}
