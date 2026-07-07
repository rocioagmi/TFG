
procesarSecuencias <- function(union, estudio) {
  
  require(dada2)
  
  directorio_rds <- file.path("OUTPUT", estudio, "RDS")
  dir.create(directorio_rds, recursive = TRUE, showWarnings = FALSE)
  
  # CONSTRUIR LA TABLA DE SECUENCIAS
  message("Construyendo tabla de secuencias")
  seqtab <- makeSequenceTable(union)
  
  message(sprintf("Dimensiones: %d muestras x %d ASV's", nrow(seqtab), ncol(seqtab)))
  message("Distribución de longitudes de amplicón:")
  print(table(nchar(getSequences(seqtab))))
  
  saveRDS(seqtab, file.path(directorio_rds, "seqtab.Rds"))

  # ELIMINAR QUIMERAS
  message("Eliminando quimeras")
  tabSinQuim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                   multithread = FALSE, verbose = TRUE)
  pct_limpias <- round(sum(tabSinQuim) / sum(seqtab) * 100, 2)
  
  message(sprintf("Lecturas conservadas tras eliminar quimeras: %s%%", pct_limpias))
  message(sprintf("ASV's: %d  - %d", ncol(seqtab), ncol(tabSinQuim)))
  
  saveRDS(tabSinQuim, file.path(directorio_rds, "tabSinQuim.Rds"))
  message("Tabla de secuencias sin quimeras guardada en %s", file.path(directorio_rds, "tabSinQuim.Rds"))
  
  return(tabSinQuim)
}