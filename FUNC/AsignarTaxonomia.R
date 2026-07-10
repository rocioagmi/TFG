
asignarTaxonomia <- function(tabSinQuim, estudio){
  
  require(dada2)
  
  dir_bbdd_genero <-"BBDD/silva_nr99_v138.2_toGenus_trainset.fa"
  dir_bbdd_especie <-"BBDD/silva_v138.2_assignSpecies.fa"
  
  if (!file.exists(dir_bbdd_genero)) {
    stop("No se encuentra la base de datos de referencia (género) en: ", dir_bbdd_genero)
  }
  if (!file.exists(dir_bbdd_especie)) {
    stop("No se encuentra la base de datos de referencia (especie) en: ", dir_bbdd_especie)
  }
  
  taxa <- assignTaxonomy(tabSinQuim, dir_bbdd_genero, multithread = FALSE)
  taxa <- assignTaxonomy(tabSinQuim, dir_bbdd_especie, multithread = FALSE)
  
  directorio_salida <- file.path("OUTPUT", estudio, "TAXONOMIA")
  dir.create(directorio_salida, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(taxa, file.path(directorio_salida, "taxonomia.Rds")) 
  
  return(taxa)
}