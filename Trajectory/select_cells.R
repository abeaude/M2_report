select_cells <- function(sobj, clusters, dimred){
  library(shiny)
  library(ggiraph)
  
  ui <- fluidPage(
    shiny::titlePanel("Select cells and get their names"),
    sidebarLayout(
      sidebarPanel(
        tags$p("This application let you select cells on the UMAP and get their names. Choose as many cells as you want adn when your finish it the DONE button "),
        tags$p("The plot might take a few seconds to get displyed"),
        actionButton("clear", "Clear"),
        actionButton("myBtn", "DONE")
      ),
      
      mainPanel(
        tags$h3("Select your parameters"),
        girafeOutput("Dimplot", width = "100%", height = "100%" ),
        tags$h3("Table previsualization"),
        verbatimTextOutput("test")
      )
    )
  )
  
  
  server <- function(input, output, session) {
    output$Dimplot <- renderGirafe({
      clust_name <- as.name(clusters)
      # dimred_name <- names(sobj@reductions) %>% stringr::str_subset("umap$")
      dimkey <- sobj@reductions[[dimred]]@key
      dimkey_1 <- as.name(paste0(dimkey,1))
      dimkey_2 <- as.name(paste0(dimkey,2))
      gg_scatter <- sobj@reductions[[dimred]]@cell.embeddings %>% 
        tibble::as_tibble(rownames = "cells") %>%
        dplyr::inner_join(tibble::rownames_to_column(sobj@meta.data[clusters],var = "cells"), by = 'cells') %>%
        ggplot2::ggplot(ggplot2::aes(x = {{dimkey_1}} , y = {{dimkey_2}}, color = {{clust_name}}, tooltip = cells, data_id = cells)) +
        geom_point_interactive(size = 4)
      girafe(ggobj = gg_scatter, 
             options = list(opts_selection(type = "multiple")), width_svg = 15, height_svg = 10 )
      
    })
    
    output$test <- renderText({paste(input$Dimplot_selected, collapse = "\n")})
    
    observeEvent(input$myBtn, {
      stopApp(input$Dimplot_selected)
    })
    observeEvent(input$clear, {
      session$sendCustomMessage(type = 'Dimplot_set', message = character(0))
    })
  }
  
  cells_app <- shinyApp(ui, server)
  cells <- runApp(cells_app)
  
  return(cells)
}


