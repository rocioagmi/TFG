library(dada2)

eliminarQuimeras <- function(seqtab){
  tabSinQuim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
  dim(tabSinQuim)
  sum(tabSinQuim)/sum(seqtab)
  
  saveRDS(tabSinQuim, paste0("OUTPUT/RDS", "/tabSinQuim.Rds"))
}