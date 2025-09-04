construirTablaSecuencias <- function(union){
  seqtab <- makeSequenceTable(union)
  dim(seqtab)
  table(nchar(getSequences(seqtab)))
  
  saveRDS(seqtab, paste0("OUTPUT/RDS","/seqtab.Rds"))
}