
filtrarMuestras <- function(muestras, parametros, estudio_id){
  
  require(dada2)
  
  R1 <- muestras$R1
  R2 <- muestras$R2
  
  if(length(R1) == 0 || length(R2) == 0) {
    stop("No se encontraron archivos de entrada para el estudio: ", estudio_id)
  }
  
  campos_requeridos <- c("truncLen", "trimLeft", "maxEE", "truncQ")
  campos_faltantes <- setdiff(campos_requeridos, names(parametros))
  
  if (length(campos_faltantes) > 0) {
    stop("El objeto 'parámetros' no contiene los campos requeridos: ",
         paste(campos_faltantes, collapse = ", "))
  }
  
  dir_filtrado <- file.path("OUTPUT", estudio_id,"FILTRADO")
  dir.create(dir_filtrado, recursive = TRUE, showWarnings = FALSE)
  
  salida_R1 <- file.path(dir_filtrado, basename(R1))
  salida_R2 <- file.path(dir_filtrado, basename(R2))
  
  cat(sprintf("\n Iniciando el filtrado %s: %d muestras\n",
              estudio_id, length(R1)))
  
  filtrado <- filterAndTrim(fwd = R1, filt = salida_R1, 
                            rev = R2, filt.rev = salida_R2, 
                            truncLen = parametros$truncLen,
                            trimLeft = parametros$trimLeft, 
                            maxEE = parametros$maxEE, 
                            truncQ = parametros$truncQ, 
                            maxN = 0, rm.phix = TRUE, 
                            compress = TRUE, multithread = FALSE)
  
  
  n_total <- nrow(filtrado)
  n_reads <- sum(filtrado[, 2] > 0)
  pct_retencion <- round(100 * sum(filtrado[, 2]) / max(sum(filtrado[, 1]), 1), 1)
  
  cat(sprintf("Filtrado completado %s: %d/%d muestras conservaron reads - %.1f%% reads retenidos.\n", 
              estudio_id, n_reads, ntotal, pct_retencion))
  
  muestras_vacias <- rownames(filtrado)[filtrado[, 2] == 0]
  if (length(muestras_vacias) > 0) {
    warning(sprintf("%s - Las siguientes muestras quedaron vacías tras el filtrado:\n",
            estudio_id, paste(muestras_vacias, collapse = "\n")))
  }
  
  ruta_params <- file.path("OUTPUT", estudio_id, "parametros_filtrado.rds")
  saveRDS(parametros, ruta_params)
  
  cat(sprintf("Parámetros de filtrado %s almacenados en: %s\n",
              estudio_id, ruta_params))
  
  return(list(R1 = nombres_R1[filtrado[,2] > 0], 
              R2 = nombres_R2[filtrado[,2] > 0],
              stats = filtrado))
}
