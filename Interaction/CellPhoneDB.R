run_cellphonedb <- function(sobj, gene_type, species = "human", threshold = 0.1, iterations = 1000, pval = 0.05, path, filename, n_threads = future::availableCores()["system"] / 2, verbose = TRUE) {
  assertthat::assert_that(is_seurat(sobj))
  assertthat::assert_that(is_in(gene_type, c("ensembl", "gene_name", "hgnc_symbol")))
  assertthat::assert_that(is.numeric(threshold), threshold >= 0, threshold <= 1)
  assertthat::assert_that(is.numeric(pval), pval >= 0, pval <= 1)
  assertthat::assert_that(is.numeric(iterations))
  assertthat::assert_that(assertthat::is.dir(path), assertthat::is.writeable(path))
  assertthat::assert_that(is.character(filename), length(filename) == 1)
  assertthat::assert_that(is.numeric(n_threads))
  assertthat::assert_that(is.logical(verbose))

  path <- fs::path(path, "INTERACTION", "cellphoneDB") %>% fs::dir_create()
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
  # if mouse need to convert the Gene name to human
  if (species == "mouse") {
    data(geneinfo_human, package = "nichenetr")

    gene_name_converted <- nichenetr::convert_mouse_to_human_symbols(rownames(sobj@assays$RNA@data))
    counts <- sobj@assays$RNA@data %>%
      tibble::as_tibble(rownames = "Gene") %>%
      dplyr::filter(Gene %in% names(gene_name_converted)) %>%
      dplyr::mutate(Gene = gene_name_converted[Gene]) -> test
  } else {
    counts <- sobj@assays$RNA@data %>% tibble::as_tibble(rownames = "Gene")
  }

  count_file <- fs::path(path, paste0(filename, "_counts.txt"))
  meta_file <- fs::path(path, paste0(filename, "_meta.txt"))

  # Now write the matrix to the disk
  counts %>% readr::write_tsv(path = count_file)
  # Write cell type information
  Seurat::Idents(sobj) %>%
    tibble::enframe("cell", "cell_type") %>%
    readr::write_tsv(path = meta_file)

  # Search for cellphonedb executable (conda env named cellphoneDB)
  cmd <- reticulate::conda_list() %>%
    dplyr::filter(name == "cellphoneDB") %>%
    dplyr::pull(python) %>%
    fs::path_dir() %>%
    fs::path("cellphonedb")

  # Create args list
  args <- c(
    "method", "statistical_analysis",
    meta_file, count_file,
    paste0("--counts-data=", gene_type),
    paste0("--output-path=", path),
    paste0("--threshold=", threshold),
    paste0("--iterations=", iterations),
    paste0("--threads=", n_threads),
    ifelse(verbose, "--verbose", "--quiet")
  )

  # Run analysis
  phonedb_out <- processx::run(
    command = cmd,
    args = args,
    echo = verbose,
    echo_cmd = verbose,
    spinner = !verbose,
    error_on_status = FALSE,
    cleanup_tree = TRUE
  )

  # plot
  # dot plot
  dot_plot(
    means_path = fs::path(path, "means.txt"),
    pvalues_path = fs::path(path, "pvalues.txt")
  ) %>% purrr::iwalk(function(plot, num) {
    pbuild <- ggplot2::ggplot_build(plot)

    y <- pbuild$layout$panel_params[[1]]$y.range %>% diff()
    x <- pbuild$layout$panel_params[[1]]$x.range %>% diff()

    ggplot2::ggsave(
      filename = paste0(filename, "_dotplots_", num, ".png"),
      plot = plot,
      width = x * 18.5, height = y * 7.4, units = "mm", limitsize = FALSE,
      path = path
    )
  })

  # heatmap plot
  # cellphonedb plot heatmap_plot test_phonedb_meta.txt
  args <- c("plot", "heatmap_plot", meta_file, paste0("--output-path=", path), paste0("--pvalue=", pval), verbose)
  processx::run(
    command = cmd,
    args = args,
    echo = FALSE,
    echo_cmd = FALSE,
    spinner = TRUE,
    error_on_status = FALSE,
    cleanup_tree = TRUE
  )
}

dot_plot <- function(
                     means_path = "./means.txt",
                     pvalues_path = "./pvalues.txt",
                     means_separator = "\t",
                     pvalues_separator = "\t") {
  all_pval <- read.table(pvalues_path, header = T, stringsAsFactors = F, sep = means_separator, comment.char = "", check.names = F)
  all_means <- read.table(means_path, header = T, stringsAsFactors = F, sep = pvalues_separator, comment.char = "", check.names = F)

  intr_pairs <- all_pval$interacting_pair
  all_pval <- all_pval[, -c(1:11)]
  all_means <- all_means[, -c(1:11)]

  selected_rows <- intr_pairs
  selected_columns <- colnames(all_pval)

  sel_pval <- all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means <- all_means[match(selected_rows, intr_pairs), selected_columns]

  df_names <- expand.grid(selected_rows, selected_columns)
  pval <- unlist(sel_pval)
  pval[pval == 0] <- 0.0009
  plot.data <- cbind(df_names, pval)
  pr <- unlist(as.data.frame(sel_means))
  pr[pr == 0] <- 1
  plot.data <- cbind(plot.data, log2(pr))
  colnames(plot.data) <- c("pair", "clusters", "pvalue", "mean")
  plot.data %<>% dplyr::mutate(clust1 = stringr::str_split(as.character(clusters), "\\|") %>% purrr::map(1) %>% unlist())
  plot.data <- dplyr::filter(plot.data, pvalue <= 0.05)

  my_palette <- grDevices::colorRampPalette(c("black", "blue", "yellow", "red"), alpha = TRUE)(n = 399)

  purrr::map(plot.data$clust1 %>% unique(), ~ plot.data %>%
    dplyr::filter(clust1 == .x) %>%
    ggplot2::ggplot(ggplot2::aes(x = clusters, y = pair)) +
    ggplot2::geom_point(ggplot2::aes(size = -log10(pvalue), color = mean)) +
    ggplot2::scale_color_gradientn("Log2 mean (Molecule 1, Molecule 2)", colors = my_palette) +
    ggplot2::theme_bw())
}
