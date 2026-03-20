library(readr)
library(httr)
library(dplyr)

descargarMuestras <- function(df, idEstudios){
  
  # --------------------------------------------------------------------------------------
  # CAMBIAR ESTA FUNCIÓN:
  # - Nombrar a las secuencias por alias
  # --------------------------------------------------------------------------------------
  
  idAcceso <- trimws(unlist(strsplit(idEstudios, ",")))
  idAcceso <- idAcceso[nchar(idAcceso) > 0]
  
  df_procesado <- df %>%
    filter(study_accession %in% idAcceso) %>%
    mutate(
      A = na_if(submitted_ftp, ""),
      B = na_if(fastq_ftp, ""),
      enlaces =  coalesce(A, B)) %>%
    arrange(enlaces)
  
  enlaces_validos <- df_procesado$enlaces
  vacios <- sum(is.na(enlaces_validos))
  
  if (vacios > 0) {
    cat(sprintf("No se pueden descargar algunas muestras."))
    enlaces_validos <- enlaces_validos[!is.na(enlaces_validos)]
  }
  
  if (length(enlaces_validos) == 0) {
    cat("No hay enlaces válidos para descargar.\n")
    return(invisible(NULL))
  }
  
  fallidos <- c()
  
  enlaces <- strsplit(enlaces_validos, ";")
  
  for (enl in enlaces){
    for (i in enl){
      i <- trimws(i)
      if (nchar(i) == 0) next
      
      destino <- strsplit(i, "/")[[1]]
      directorioEnl <- paste0("INPUT/DATA/",destino[length(destino)])
      tryCatch({
        download.file(i,directorioEnl, mode = "wb")
      }, error = function(e) {
        cat(sprintf(" Error en %s: %s\n", directorioEnl, e$message))
        fallidos <<- c(fallidos, i)
      })
    }
  }
  
  cat(sprintf("\n--- DESCARGA COMPLETADA ---\n"))
  cat(sprintf("Archivos descargados: %d\n", length(unlist(enlaces)) - length(fallidos)))
  if (length(fallidos) > 0) {
    cat(sprintf("Archivos fallidos (%d):\n", length(fallidos)))
    for (f in fallidos) cat(sprintf("  - %s\n", f))
  }
}
