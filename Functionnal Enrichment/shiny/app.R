library(shiny)
library(shinyjs, warn.conflicts = FALSE)
library(shinyWidgets, warn.conflicts = FALSE)

ui <- fluidPage(
  useShinyjs(),
  shiny::titlePanel("Create your set of parameters to run differential expression"),
  sidebarLayout(
    sidebarPanel(
      tags$p("This application let you define a set of different parameters to run differential expression on your seurat object.
             Three types of analysis are available : "),
      tags$ol(
        tags$li(tags$b("1vs1"), " : 1 cluster versus another cluster"),
        tags$li(tags$b("1vsAll"), " : 1 cluster vs all the other clusters"),
        tags$li(tags$b("SvsS"), " : A group of cluster versus another group of cluster"),
        tags$li(tags$b("conditions"), " : Compare a group of cells beetwen 2 conditions")
      ),
      tags$p("When selecting 1vsAll, you can only specify the first group. When selecting SvsS first group and second group cannot overlap"),
      tags$h3("Example:"),
      tags$p("For example if you select cluster 1,2,3 in the first group and 2,3 in the second group and 1vs1 as method of comparisons, the algorithm will run the 4 comparisons : 1vs2, 1vs3, 2vs3 and 3vs2  (homotypic comparisons are removed : 2vs2 and 3vs3), for each choosen methods/test."),
      tags$h3("Instructions:"),
      tags$ol(
        tags$li("Select the type of differential expression you want to perform"),
        tags$li("Select the cluster you want to include in the first group"),
        tags$li("Select the cluster you want to include in the second group. You can include cluster selected in the first group (depending of the choosen DE_type."),
        tags$li("Choose the Seurat test and/or method (EdgeR/Limma) you want to use"),
        tags$li("Click 'Add' to add the chosen parameters"),
        tags$li("When finished Click 'Done'"),
        tags$li("You can now close the tab")
      ),
      actionButton("add", "Add"),
      actionButton("rm", "Remove"),
      actionButton("clear", "Clear"),
      actionButton("myBtn", "DONE")
    ),

    mainPanel(
      tags$h3("Select your parameters"),
      selectInput("DE_type", label = "DE type", choices = c("1vs1", "1vsAll", "SvsS", "conditions"), multiple = FALSE),
      selectInput("Grouping", label = "Grouping (clustering or conditions)", choices = colnames(sobj@meta.data), multiple = FALSE, selected = clusters),
      tabsetPanel(
        id = "pan_switch",
        type = "hidden",
        selected = "empty_pan",
        tabPanel("cluster_cond_panel", checkboxGroupInput("cluster_cond", label = "Cluster conditions", choices = levels(sobj@meta.data[["seurat_clusters"]]), inline = TRUE)),
        tabPanel("empty_pan", tags$p(""))
      ),
      checkboxGroupInput("first_group", label = "First Group", choices = levels(sobj@meta.data[[clusters]]), inline = TRUE),
      checkboxGroupInput("second_group", label = "Second Group", choices = levels(sobj@meta.data[[clusters]]), inline = TRUE),
      selectInput("batch", label = "Batch effect", choices = colnames(sobj@meta.data), multiple = TRUE),
      fluidRow(
        column(width = 4, selectInput("method", label = "Method", choices = c("Limma_Voom", "Limma_Trend", "EdgeR_QL", "EdgeR_LRT"), multiple = TRUE)),
        column(width = 4, selectInput("Seurat_test", label = "Seurat test", choices = c("wilcox", "bimod", "t", "poisson", "negbinom", "LR", "MAST"), multiple = TRUE))
      ),
      tags$h3("Table previsualization"),
      tableOutput("data")
    )
  )
)

server <- function(input, output, session) {
  values <- reactiveValues()
  values$data_print <- NULL
  values$warn <- TRUE
  values$cluster <- NULL

  observeEvent(input$DE_type,
    {
      if (input$DE_type == "1vsAll") {
        choice_scnd_grp <- character(0)
        choice_first_grp <- isolate(values$cluster)
        updateSelectInput(session, "Grouping", selected = "seurat_clusters")
        updateTabsetPanel(session, "pan_switch", selected = "empty_pan")
      } else if (input$DE_type == "conditions") {
        updateTabsetPanel(session, "pan_switch", selected = "cluster_cond_panel")
        updateSelectInput(session, "Grouping", selected = "conditions")
        choice_scnd_grp <- NULL
        choice_first_grp <- NULL
      } else {
        choice_scnd_grp <- isolate(values$cluster)
        choice_first_grp <- isolate(values$cluster)
        updateSelectInput(session, "Grouping", selected = "seurat_clusters")
        updateTabsetPanel(session, "pan_switch", selected = "empty_pan")
      }
      updateCheckboxGroupInput(session, "second_group",
        label = "Second group",
        choices = choice_scnd_grp,
        selected = NULL,
        inline = TRUE
      )
      updateCheckboxGroupInput(session, "first_group",
        label = "First group",
        choices = choice_first_grp,
        selected = NULL,
        inline = TRUE
      )
      to_reset <- c("first_group", "second_group", "method", "Seurat_test", "cluster_cond", "batch")
      purrr::map(to_reset, reset)
    },
    ignoreInit = TRUE
  )

  observeEvent(input$Grouping, {
    isolate(values$cluster <- levels(sobj@meta.data[[input$Grouping]]))
    if (input$DE_type == "1vsAll") {
      choice_first_group <- isolate(values$cluster)
      choice_scnd_group <- character(0)
    } else {
      choice_first_group <- isolate(values$cluster)
      choice_scnd_group <- isolate(values$cluster)
    }
    updateCheckboxGroupInput(session, "first_group",
      label = "First group",
      choices = choice_first_group,
      selected = NULL,
      inline = TRUE
    )
    updateCheckboxGroupInput(session, "second_group",
      label = "Second group",
      choices = choice_scnd_group,
      selected = NULL,
      inline = TRUE
    )

    to_reset <- c("first_group", "second_group", "method", "Seurat_test", "cluster_cond", "batch")
    purrr::map(to_reset, reset)
  })

  observeEvent(input$first_group, {
    if (input$DE_type == "SvsS") {
      choice_scnd_group <- setdiff(isolate(values$cluster), input$first_group)
      updateCheckboxGroupInput(session, "second_group",
        label = "Second group",
        choices = choice_scnd_group,
        selected = input$second_group,
        inline = TRUE
      )
    }
  })

  observeEvent(input$second_group, {
    if (input$DE_type == "SvsS") {
      choice_first_group <- setdiff(isolate(values$cluster), input$second_group)
      updateCheckboxGroupInput(session, "first_group",
        label = "First group",
        choices = choice_first_group,
        selected = input$first_group,
        inline = TRUE
      )
    }
  })

  observeEvent(input$add, {
    null_to_na <- function(x) {
      if (is.null(x)) {
        return(NA)
      } else {
        return(x)
      }
    }
    second_group <- null_to_na(input$second_group)
    method <- null_to_na(input$method)
    Seurat_test <- null_to_na(input$Seurat_test)
    cluster_conditions <- null_to_na(input$cluster_cond)
    batch <- null_to_na(input$batch)

    data_p <- tibble::tibble(
      DE_type = input$DE_type,
      first_group = input$first_group %>% stringr::str_c(collapse = ","),
      second_group = second_group %>% stringr::str_c(collapse = ","),
      method = method %>% stringr::str_c(collapse = ","),
      testToUse = Seurat_test %>% stringr::str_c(collapse = ","),
      cluster_conditions = cluster_conditions %>% stringr::str_c(collapse = ","),
      batch = batch %>% stringr::str_c(collapse = ","),
      clusters = input$Grouping
    )

    isolate(values$data_print <- rbind(values$data_print, data_p))
    output$data <- renderTable(values$data_print)
    to_reset <- c("first_group", "second_group", "method", "Seurat_test", "cluster_cond", "batch")
    purrr::map(to_reset, reset)
  })

  observeEvent(input$method, {
    if ("EdgeR_LRT" %in% input$method & values$warn) {
      sendSweetAlert(
        session,
        title = "Warning",
        text = "You choose the method EdgeR_LRT, this method take more than 1 hour to run",
        type = "warning",
        btn_labels = "Ok",
        closeOnClickOutside = TRUE,
        showCloseButton = FALSE,
        width = NULL
      )
      isolate(values$warn <- FALSE)
    }
  })

  observeEvent(input$myBtn, {
    stopApp(values$data_print)
  })
  observeEvent(input$clear, {
    isolate(values$data_print <- NULL)
  })

  observeEvent(input$rm, {
    if (!is.null(values$data_print)) {
      if (nrow(values$data_print) > 0) {
        isolate(values$data_print %<>% dplyr::slice(-nrow(.)))
      } else {
        isolate(values$data_print <- NULL)
      }
      # soutput$data <- renderTable(values$data_print)
    }
  })
}

multiple_DE_app <- shinyApp(ui, server)
