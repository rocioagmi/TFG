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
    filter(study_accession %in% idAccesso) %>%
    mutate(
      A = na_if(submitted_ftp, ""),
      B = na_if(fastq_ftp, ""),
      enlaces <-  coalesce(A, B)) %>%
    arrange(enlaces)
  
  enlaces_validos <- df_procesado$enlaces
  vacios <- sum(is.na(enlaces_validos))
  
  if (vacios > 0) {
    sprintf("No se pueden descargar algunas muestras.")
    enlaces_validos <- enlaces_validos[!is.na(enlaces_validos)]
  }
  
  enlaces <- strsplit(enlaces_validos, ";")
  
  for (enl in enlaces){
    for (i in enl){
      i <- trimws(i)
      destino <- strsplit(i, "/")[[1]]
      directorioEnl <- paste0("INPUT/DATA/",destino[length(destino)])
      download.file(i,directorioEnl)
    }
  }
}
