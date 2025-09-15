asignarTaxonomia <- function(tabSinQuim){
  directorioBBDD <- dlgInput(message = "Ruta completa al directorio con la Base de Datos: ")$res
  taxa <- assignTaxonomy(tabSinQuim, directorioBBDD, multithread = TRUE)
  
  texto <- dlgInput(message = "Introduzca el nombre del archivo (sin extensión)")$res
  nombre <- paste0("OUTPUT/RDS/", texto, ".Rds")
  saveRDS(taxa, nombre)  
}