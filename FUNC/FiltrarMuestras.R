
filtrarMuestras <- function(muestras, parametros, estudio_id){
  
  require(dada2)
  
  R1 <- muestras$R1
  R2 <- muestras$R2
  
  if(length(R1) == 0 || length(R2) == 0) {
    stop("No se encontraron archivos de entrada para el estudio: ", estudio_id)
  }
  
  if (length(R1) != length(R2)) {
    stop("El número de archivos R1 y R2 no coincide para el estudio: ", estudio_id)
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
              estudio_id, n_reads, n_total, pct_retencion))
  
  muestras_vacias <- rownames(filtrado)[filtrado[, 2] == 0]
  if (length(muestras_vacias) > 0) {
    warning(sprintf("%s - Las siguientes muestras quedaron vacías tras el filtrado: %s\n",
            estudio_id, paste(muestras_vacias, collapse = "\n")))
  }
  
  ruta_params <- file.path("OUTPUT", estudio_id, "parametros_filtrado.rds")
  saveRDS(parametros, ruta_params)
  
  cat(sprintf("Parámetros de filtrado %s almacenados en: %s\n",
              estudio_id, ruta_params))
  
  return(list(R1 = salida_R1[filtrado[,2] > 0], 
              R2 = salida_R2[filtrado[,2] > 0],
              stats = filtrado))
}


proponerParametros <- function(datos_qa, umbral_calidad) {
  
  require(dplyr)
  
  campos_necesarios <- c("calidad_ciclo", "estadisticas")
  campos_faltantes  <- setdiff(campos_necesarios, names(datos_qa))
  if (length(campos_faltantes) > 0) {
    stop("datos_qa no contiene los campos requeridos: ",
         paste(campos_faltantes, collapse = ", "))
  }
  
  # Mediana de calidades por hebra y ciclo
  resumen_calidad <- datos_qa$calidad_ciclo %>%
    group_by(hebra, ciclo) %>%
    summarise(q10_md = median(q10, na.rm = TRUE),
              q25_md = median(q25, na.rm = TRUE),
              mediana_md = median(mediana, na.rm = TRUE),
              ee_mediana_md = median(ee_mediana, na.rm = TRUE),
              .groups = "drop"
    )
  
  # Longitud max por hebra
  max_R1 <- max(datos_qa$estadisticas$longitud_max[
    datos_qa$estadisticas$hebra == "R1"], na.rm = TRUE)
  max_R2 <- max(datos_qa$estadisticas$longitud_max[
    datos_qa$estadisticas$hebra == "R2"], na.rm = TRUE)
  
  # Funcion interna devuelve el último ciclo donde la curva >= umbral
  ultimo_ciclo <- function(ciclos, valores, umbral) {
    ciclos_buenos <- ciclos[valores >= umbral]
    if (length(ciclos_buenos) == 0) {
      warning("La curva de calidad no supera el umbral en ningún ciclo.")
      return(1)
    }
    max(ciclos_buenos)
  }
  
  # Cálculo truncLen & truncLen_rango
  for_R1 <- resumen_calidad[resumen_calidad$hebra == "R1", ]
  for_R2 <- resumen_calidad[resumen_calidad$hebra == "R2", ]
  
  trunc_R1 <- ultimo_ciclo(for_R1$ciclo, for_R1$q25_md, umbral_calidad)
  trunc_R2 <- ultimo_ciclo(for_R2$ciclo, for_R2$q25_md, umbral_calidad)
  
  rango_R1 <- c(
    min = ultimo_ciclo(for_R1$ciclo, for_R1$q10_md, umbral_calidad),
    max = ultimo_ciclo(for_R1$ciclo, for_R1$mediana_md, umbral_calidad)
  )
  rango_R2 <- c(
    min = ultimo_ciclo(for_R2$ciclo, for_R2$q10_md, umbral_calidad),
    max = ultimo_ciclo(for_R2$ciclo, for_R2$mediana_md, umbral_calidad)
  )
  
  # EE estimado en el pto de corte propuesto
  ee_R1 <- for_R1$ee_mediana_md[for_R1$ciclo == trunc_R1]
  ee_R2 <- for_R2$ee_mediana_md[for_R2$ciclo == trunc_R2]
  
  if (length(ee_R1) == 0) ee_R1 <- NA
  if (length(ee_R2) == 0) ee_R2 <- NA
  
  # Parámetros filtrado
  parametros_filtrado <- list(truncLen = c(R1 = trunc_R1, R2 = trunc_R2),
                              truncLen_rango = list(R1 = rango_R1, R2 = rango_R2),
                              trimLeft = c(R1 = 0,  R2 = 0),
                              maxEE = c(R1 = 2,   R2 = 2),
                              truncQ = 2,
                              max_R1 = max_R1,
                              max_R2 = max_R2,
                              ee_estimado_truncLen = c(R1 = round(ee_R1, 2), 
                                                       R2 = round(ee_R2, 2)))
  
  return(parametros_filtrado)
}


confirmarParametros <- function(parametros, interactivo = TRUE) {
  
  require(svDialogs)
  
  if (!interactivo) {
    return(parametros)
  }
  
  parsear_par <- function(texto) {
    partes <- strsplit(trimws(texto), "[,.\\s]+")[[1]]
    if (length(partes) != 2) return(NULL)
    valores <- suppressWarnings(as.integer(partes))
    if (any(is.na(valores))) return(NULL)
    if (any(valores < 0)) return(NULL)
    return(valores)
  }
  
  # truncLen y el EE estimado en ese punto
  ee_R1 <- parametros$ee_estimado_truncLen["R1"]
  ee_R2 <- parametros$ee_estimado_truncLen["R2"]
  
  aviso_ee <- ""
  if (!is.na(ee_R2) && ee_R2 > parametros$maxEE["R2"]) {
    aviso_ee <- sprintf(
      "\nAVISO: EE mediano de R2 en el punto de corte: %.2f (supera maxEE = %.0f).\n  Considera acortar R2 para mejorar la retención de reads.",
      ee_R2, parametros$maxEE["R2"])
  }
  
  mensaje_truncLen <- sprintf(
    paste0(
      "Estudio: longitud máxima de reads — R1: %d pb | R2: %d pb\n\n",
      "Valor propuesto:  R1 = %d | R2 = %d\n",
      "Rango recomendado: R1 [%d – %d] | R2 [%d – %d]\n",
      "EE mediano estimado en el punto de corte: R1 = %.2f | R2 = %.2f",
      "%s\n\n",
      "Introduce dos valores separados por coma (R1,R2):"
    ),
    parametros$max_R1, parametros$max_R2,
    parametros$truncLen["R1"], parametros$truncLen["R2"],
    parametros$truncLen_rango$R1["min"], parametros$truncLen_rango$R1["max"],
    parametros$truncLen_rango$R2["min"], parametros$truncLen_rango$R2["max"],
    ee_R1, ee_R2,
    aviso_ee
  )
  
  repeat {
    respuesta <- dlgInput(message = mensaje_truncLen, 
                          default = paste(parametros$truncLen["R1"], 
                                          parametros$truncLen["R2"], sep = ","))$res
    
    if (is.null(respuesta) || respuesta == "") {
      stop("Parámetro truncLen no confirmado. Filtrado cancelado por el usuario.")
    }
    
    valores <- parsear_par(respuesta)
    
    if (is.null(valores)) {
      dlg_message(sprintf("Entrada no válida: '%s'.\nIntroduce exactamente dos números enteros separados por coma (ej: 230,160).",
                  respuesta), type = "ok")
      next
    }
    
    if (valores[1] > parametros$max_R1 || valores[2] > parametros$max_R2) {
      dlg_message(sprintf(
        "Los valores introducidos (%d, %d) superan la longitud máxima de los reads (R1: %d, R2: %d).\nIntroduce valores dentro del rango.",
        valores[1], valores[2], parametros$max_R1, parametros$max_R2), type = "ok")
      next
    }
    
    break
  }
  truncLen_confirmado <- valores
  
  # trimLeft
  mensaje_trimLeft <- paste0(
    "Eliminación de nucleótidos en el extremo 5' de cada read (trimLeft).\n\n",
    "Valor propuesto: R1 = 0 | R2 = 0\n\n",
    "Si conoces la longitud de los primers utilizados en este estudio,\n",
    "introdúcela aquí para eliminarlos antes del filtrado.\n",
    "Si no la conoces, acepta el valor 0,0.\n\n",
    "Introduce dos valores separados por coma (R1,R2):"
  )
  
  repeat {
    respuesta <- dlgInput(
      message = mensaje_trimLeft,
      default = paste(parametros$trimLeft["R1"], parametros$trimLeft["R2"], sep = ",")
    )$res
    
    if (is.null(respuesta) || respuesta == "") {
      stop("Parámetro trimLeft no confirmado. Filtrado cancelado por el usuario.")
    }
    
    valores <- parsear_par(respuesta)
    
    if (is.null(valores)) {
      dlg_message(sprintf(
        "Entrada no válida: '%s'.\nIntroduce exactamente dos números enteros separados por coma (ej: 19,20).",
        respuesta), type = "ok")
      next
    }
    
    break
  }
  trimLeft_confirmado <- valores
  
  parametros$truncLen <- c(R1 = truncLen_confirmado[1], R2 = truncLen_confirmado[2])
  parametros$trimLeft <- c(R1 = trimLeft_confirmado[1], R2 = trimLeft_confirmado[2])
  
  return(parametros)
}