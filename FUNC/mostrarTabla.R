mostrarTabla <- function(df){
  
  require(DT)
  require(shiny)
  require(data.table)
  
  if (is.null(df) || nrow(df) == 0) {
    stop("El dataframe está vacío. No hay resultados para visualizar.")
  }
  
  setDT(df)
  
  ui <- fluidPage(
    titlePanel("Resultados de Búsqueda - Microbiota Intestinal en Esclerosis Múltiple"),
    
    tags$h4("Explora y selecciona las filas que deseas descargar"),
    tags$p("Si no seleccionas ninguna fila y continúas con descarga, se usarán TODAS las filas."),
    
    fluidRow(column(12, DTOutput("tabla", height = "650px"))),
    
    br(),
    
    fluidRow(column(8, verbatimTextOutput("info_seleccion")),
             column(4, actionButton("btn_continuar", "Descargar filas seleccionadas", 
                                    class = "btn-lg btn-primary", width = "100%"))))


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
          lengthMenu = c(25, 50, 100, 250, 500, 1000),
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel"),
          scrollY = 520,
          scrollX = TRUE,
          deferRender = TRUE,
          scroller = TRUE,
          fixedHeader = TRUE,
          stateSave = TRUE))
      }, server = TRUE)
    
    output$info_seleccion <- renderText({
      n <- length(input$tabla_rows_selected)
      n_total <- nrow(df)
      
      if (n == 0) {
        paste0("Ninguna fila seleccionada, se descargarán ", n_total, " filas.")
      } else {
        paste0(n, " fila(s) seleccionada(s) de ", n_total, " totales")
      }
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
