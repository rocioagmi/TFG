mostrarTabla <- function(df){
  
  require(DT)
  require(shiny)
  #require(data.table)
  
  if (is.null(df) || nrow(df) == 0) {
    stop("El dataframe está vacío. No hay resultados para visualizar.")
  }
  
  #setDT(df)
  
  ui <- fluidPage(
    titlePanel("Resultados de Búsqueda - Microbiota Intestinal en Esclerosis Múltiple"),
    
    tags$h4("Explora y selecciona los estudios que deseas descargar"),
    tags$p("Usa los filtros de la tabla para acotar resultados.
            Si no seleccionas ninguna fila, se usarán todos los estudios visibles."),
    
    fluidRow(column(12, DTOutput("tabla", height = "600px"))),
    
    br(),
    
    fluidRow(column(8, verbatimTextOutput("info_seleccion")),
             column(4, uiOutput("btn_dinamico"))
    )
  )

  server <- function(input, output, session) {
    
    output$tabla <- renderDT({
      datatable(
        df,
        rownames = FALSE,
        selection = "multiple",
        filter = "top",
        extensions = c("Buttons", "Scroller", "FixedHeader"),
        options = list(
          pageLength = 50,
          lengthMenu = c(25, 50, 100, 250, 500),
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel"),
          scrollY = 480,
          scrollX = TRUE,
          deferRender = TRUE,
          scroller = TRUE,
          fixedHeader = TRUE,
          stateSave = FALSE))
      }, server = FALSE)
    
    output$info_seleccion <- renderText({
      n_sel <- length(input$tabla_rows_selected)
      n_visible <- length(input$tabla_rows_all)
      n_total <- nrow(df)
      
      if (n_sel == 0) {
        sprintf("Visibles: %d de %d | Sin selección → se usarán los %d visibles",
                n_visible, n_total, n_visible)
      } else {
        sprintf("Visibles: %d de %d | Seleccionados: %d",
                n_visible, n_total, n_sel)
      }
    })
    
    output$btn_dinamico <- renderUI({
      n_sel <- length(input$tabla_rows_selected)
      n_visible <- length(input$tabla_rows_all)
      
      label <- if (n_sel == 0) {
        sprintf("Continuar con %d estudio(s) visible(s)", n_visible)
      } else {
        sprintf("Continuar con %d estudio(s) seleccionado(s)", n_sel)
      }
      
      actionButton("btn_continuar", label, class = "btn-lg btn-primary", width = "100%")
    })
    
    observeEvent(input$btn_continuar, {
      
      filas_a_usar <- if (length(input$tabla_rows_selected) == 0) {
        input$tabla_rows_all  
      } else {
        input$tabla_rows_selected
      }
      
      datos_seleccionados <- df[filas_a_usar, ]
      stopApp(returnValue = datos_seleccionados)
    })
  }
 
  seleccion <- shiny::runApp(list(ui = ui, server = server), launch.browser = TRUE)
  return(seleccion)
}
