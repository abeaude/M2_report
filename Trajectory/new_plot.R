# Can work on asimple wrap trajectory with dynverse
# or with a complex list like the one returned by run ti
add_gene_importances <- function(trajectories, qval_thres = 0.05, n_perm = 5, ntree = 10000, ntree_perm = ntree / 10, n_threads = 8, adjust.methods = "BH",pb = NULL) {
  assertthat::assert_that(qval_thres >= 0, qval_thres < 1)
  assertthat::assert_that(is.numeric(n_perm), is.numeric(ntree), is.numeric(ntree_perm), is.numeric(n_threads))
  assertthat::assert_that(is_in(adjust.methods,p.adjust.methods))
  if (list_depth(trajectories) >= 1) {
    purrr::map(trajectories, add_gene_importances, qval_thres = qval_thres, n_perm = n_perm, ntree = ntree, ntree_perm = ntree_perm, n_threads = n_threads, adjust.methods = adjust.methods, pb = pb)
  } else if (list_depth(trajectories) == 0 & dynwrap::is_data_wrapper(trajectories)) {
    # only if linear trajectory 
    gr <- igraph::graph_from_data_frame(
      d = trajectories$milestone_network %>% dplyr::rename(weight = length),
      directed = any(trajectories$milestone_network$directed),
      vertices = trajectories$milestone_ids
    )
    linear_degree <- (igraph::vcount(gr) - 2) * 2 + 2
    total_degree <- igraph::degree(gr) %>% sum()
    
    if(total_degree == linear_degree){ # only linear trajectory otherwise do it by branch ??
      gene_imp <- SCORPIUS::gene_importances(as.matrix(trajectories$expression), trajectories$pseudotime, num_permutations = n_perm, ntree = ntree, ntree_perm = ntree_perm, num_threads = n_threads, verbose = FALSE)
      gene_imp$qvalue <- p.adjust(gene_imp$pvalue, adjust.methods, length(gene_imp$pvalue))
      gene_sel <- gene_imp$gene[gene_imp$qvalue < qval_thres]
      trajectories[["gene_importances"]] <- gene_sel
    } # else : by lineages ? 
    
    if(!is.null(pb)){
      pb$tick(1)
    }
    return(trajectories)
  } else {
    rlang::abort("Change ME")
  }
}

heatmap_trajectory <- function(expression, time, label = NULL, genes_highlight = NULL) {

  # checks input
  if (is(expression, "sparseMatrix")) { # some bugs with sparse matrix
    expression <- as.matrix(expression)
  }
  assertthat::assert_that(is.matrix(expression))
  assertthat::assert_that(is.numeric(time))
  assertthat::assert_that(rlang::is_named(time), msg = "time is not a named vector")
  assertthat::assert_that(identical(rownames(expression), names(time)))
  if (!is.null(label)) {
    assertthat::assert_that(is.character(label))
  }
  if (!is.null(genes_highlight)) {
    assertthat::assert_that(is.character(genes_highlight))
  }
  # order cell according to time
  cell_order <- time %>%
    sort() %>%
    names()

  # Identify different modules to order the genes along the trajectory
  modules <- SCORPIUS::extract_modules(dynutils::scale_quantile(expression), time, verbose = FALSE)
  ordered_gene <- dplyr::pull(modules, feature)

  # Reshape the expression matrix to a tibble usable by
  # ggplot2. gene and cells are transformed to factor
  # to preserve their order when plotted
  expression <- dynutils::scale_quantile(expression) %>%
    tibble::as_tibble(rownames = "cell") %>%
    tidyr::pivot_longer(cols = -cell, names_to = "gene", values_to = "expression") %>%
    dplyr::mutate(cell = factor(cell, cell_order), gene = factor(gene, ordered_gene))

  # Principal plot : heatmap
  plot <- ggplot2::ggplot(expression, ggplot2::aes(x = cell, y = gene, fill = expression)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(name = "Scaled Expression", colours = c("blue", "yellow", "red"), na.value = "white") +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin()
    )

  # Create a color bar to represent the pseudotime evolution
  # of the cells
  time_color <- scales::cscale(sort(time), palette = scales::seq_gradient_pal("darkblue", "red")) %>% t()

  bar <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotation_raster(time_color, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(plot.margin = ggplot2::margin())

  # Label the 2 extremities of the time bar
  txt1 <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = label[1], fontface = "bold") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), expand = c(0, 0),
      breaks = NULL, labels = NULL, name = NULL
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), expand = c(0, 0),
      breaks = NULL, labels = NULL, name = NULL
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin()
    )

  txt2 <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = label[2], fontface = "bold") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), expand = c(0, 0),
      breaks = NULL, labels = NULL, name = NULL
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), expand = c(0, 0),
      breaks = NULL, labels = NULL, name = NULL
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin()
    )

  # On y-axis only display the name of the genes supplied by the user
  if (!is.null(genes_highlight)) {
    genes <- levels(expression$gene)
    n_genes <- expression$gene %>% nlevels()
    genes[which(!levels(expression$gene) %in% genes_highlight)] <- ""

    axis <- ggplot2::ggplot(
      data.frame(
        y = 1:n_genes,
        gene = genes
      ),
      ggplot2::aes(x = 0, y = y, label = gene)
    ) +
      ggrepel::geom_text_repel(
        min.segment.length = grid::unit(0, "pt"),
        color = "grey30", ## ggplot2 theme_grey() axis text
        size = 0.8 * 11 / ggplot2::.pt ## ggplot2 theme_grey() axis text
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1), expand = c(0, 0),
        breaks = NULL, labels = NULL, name = NULL
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0.5, n_genes + .5), expand = c(0, 0),
        breaks = NULL, labels = NULL, name = NULL
      ) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        plot.margin = ggplot2::margin()
      )
  } else { # otherwise create an empty plot
    axis <- patchwork::plot_spacer()
  }

  # Wrap all the plot together
  patchwork::wrap_plots(txt1, 
    bar,
    txt2, 
    patchwork::plot_spacer(),
    plot,
    axis,
    ncol = 3,
    nrow = 2,
    guides = "collect",
    heights = c(0.5, 5),
    widths = c(.5, 5, 0.5)
  )
}
