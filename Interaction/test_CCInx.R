run_CCInx <- function(path, test, clusters, species) {
  assertthat::assert_that(assertthat::is.dir(path))
  assertthat::assert_that(assertthat::is.readable(path))
  assertthat::assert_that(is.character(clusters))
  fs::path(path, "INTERACTION", "CCInx") %>% fs::dir_create()
  species <- switch(species,
    mouse = "mmusculus",
    human = "hsapiens"
  )
  DE_type <- "1vsAll"

  # check if DE_folder exist
  dir_exist <- fs::path(path, "DIFFERENTIAL_EXPRESSION", DE_type) %>% fs::dir_exists()

  if (!dir_exist) {
    stringr::str_c("Differential expression was not run. Please run differential expression with DE_type set to ", DE_type) %>%
      rlang::abort()
  }

  available_test <- fs::path(path, "DIFFERENTIAL_EXPRESSION", DE_type) %>%
    fs::dir_ls() %>%
    fs::path_file()

  if (!is.null(test)) {
    assertthat::assert_that(is.character(test))
  } else {
    test <- available_test
  }

  if (any(!test %in% available_test)) {
    wrong_test <- test[!test %in% available_test]
    test <- setdiff(test, wrong_test)
    warning("The following test are not available, they will be ignored :\n", stringr::str_c(" - ", wrong_test, collapse = "\n"))
  }

  # read the different test
  DE_data <- purrr::map_df(test, ~ fs::path(path, "DIFFERENTIAL_EXPRESSION", DE_type, .x, paste(.x, DE_type, "markers", sep = "_"), ext = "txt") %>%
    readr::read_table2(col_types = test_cols_spec(.x)) %>%
    dplyr::select(gene, logFC, adj.P.Val, tested_cluster)) # not filtering allow information in plot on p.value

  # select unique rows based on genes and cluster
  DE_data %<>% dplyr::distinct(gene, tested_cluster, .keep_all = TRUE)

  # check clusters
  # if clusters is missing
  available_clust <- DE_data[["tested_cluster"]] %>%
    unique() %>%
    as.character()

  wrong_clust <- !all(clusters %in% available_clust)

  if (wrong_clust) {
    stringr::str_c(
      "One or more of the specified cluster does not exist in differential analysis results.\n",
      "Available clusters are : ",
      stringr::str_c(available_clust, collapse = ",")
    ) %>%
      rlang::abort()
  }

  # make it a named list by cluster
  GeneStatList <- available_clust %>%
    purrr::set_names(., nm = .) %>%
    purrr::map(~ dplyr::filter(DE_data, tested_cluster == .) %>% tibble::column_to_rownames("gene")) %>%
    purrr::keep(names(.) %in% clusters) # only the requested cluster by the user

  interaction_network <- CCInx::BuildCCInx(GeneStatList, GeneMagnitude = "logFC", GeneStatistic = "adj.P.Val", Species = species)
  return(interaction_network)
}

test_cols_spec <- function(test) {
  if (any(c("wilcox", "bimod", "t", "poisson", "negbinom", "LR", "MAST", "roc", "DEseq2") == test)) {
    readr::cols(
      p_val = readr::col_number(),
      logFC = readr::col_number(),
      pct.1 = readr::col_number(),
      pct.2 = readr::col_number(),
      adj.P.Val = readr::col_number(),
      tested_cluster = readr::col_factor(),
      control_cluster = readr::col_factor(),
      gene = readr::col_character(),
      min.pct = readr::col_number(),
      cluster_conditions = readr::col_character()
    )
  } else if (stringr::str_detect(test, "EdgeR")) {
    readr::cols(
      logFC = readr::col_number(),
      logCPM = readr::col_number(),
      F = readr::col_number(),
      adj.P.Val = readr::col_number(),
      FDR = readr::col_number(),
      tested_cluster = readr::col_factor(),
      control_cluster = readr::col_factor(),
      gene = readr::col_character(),
      min.pct = readr::col_number(),
      cluster_conditions = readr::col_character()
    )
  } else if (stringr::str_detect(test, "Limma")) {
    readr::cols(
      logFC = readr::col_number(),
      AveExpr = readr::col_number(),
      t = readr::col_number(),
      P.Value = readr::col_number(),
      adj.P.Val = readr::col_number(),
      B = readr::col_number(),
      tested_cluster = readr::col_factor(),
      control_cluster = readr::col_factor(),
      gene = readr::col_character(),
      min.pct = readr::col_number(),
      cluster_conditions = readr::col_character()
    )
  }
}

# path <- "/home/aurel/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/Project/NF1-MPNST/NF1-MPNST_20190212/"
#
# DE_type <- "1vsAll"
# test <- c("t", "bimod") # as many as you want

# Result visualization
# library(shiny)
# CCInx::ViewCCInx(test_ccinx)
