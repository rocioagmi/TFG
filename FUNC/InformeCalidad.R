informeCalidadSR <- function(directorioMuestras){
  qaSummary <- qa(directorioMuestras, type = "fastq")
  browseURL(report(qaSummary)) 
}

informeCalidadR <- function(directorioMuestras){
  batch_size <- 50
  batches <- split(directorioMuestras, ceiling(seq_along(directorioMuestras) / batch_size))
  qa_list <- list()
  
  for (i in seq_along(batches)) {
    cat("Procesando lote", i, "de", length(batches), "\n")
    qa_list[[i]] <- rqcQA(batches[[i]], sample = TRUE, n = 1e6, workers = 1)
  }
  
  qa <- do.call(c, qa_list)
  
  rqcReport(qa, outdir = "INPUT/DATA/REPORT", openBrowser = TRUE)
}