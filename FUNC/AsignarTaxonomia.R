asignarTaxonomia <- function(tabSinQuim){
  
  taxa <- assignTaxonomy(tabSinQuim, directorioBBDD, multithread = TRUE)
  saveRDS(taxa, paste0("OUTPUT/RDS", "/taxHit.Rds"))  
}