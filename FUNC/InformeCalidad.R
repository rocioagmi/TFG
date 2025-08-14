informeCalidadSR <- function(directorioMuestras){
  qaSummary <- qa(directorioMuestras, type = "fastq")
  browseURL(report(qaSummary)) 
}

informeCalidadR <- function(directorioMuestras){
  param <- SnowParam(workers = 4, type = "SOCK")
  register(param)
  
  batch_size <- 50
  batches <- split(directorioMuestras, ceiling(seq_along(directorioMuestras) / batch_size))
  
  qa_list <- list()
  
  for (i in seq_along(batches)) {
    cat("Procesando lote", i, "de", length(batches), "\n")
    qa_list[[i]] <- rqcQA(batches[[i]], sample = TRUE, n = 1e6)
  }
  
  qa <- do.call(c, qa_list)
  
  rqcReport(qa, outdir = "INPUT/DATA/REPORT")
  browseURL(file.path("INPUT/DATA/REPORT", "rqc_report.html"))
}