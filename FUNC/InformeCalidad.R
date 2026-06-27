informeCalidad <- function(estudio_id, muestras, umbral_calidad = 20, umbral_ee = 2, abrir_navegador = interactive()){
  
  require(rmarkdown)
  
  cat(sprintf("[%s] Procesando %d archivos.\n", estudio_id, length(muestras$R1) + length(muestras$R2)))
  
  qa_R1 <- extraerCalidad(muestras$R1, hebra = "R1")
  qa_R2 <- extraerCalidad(muestras$R2, hebra = "R2")
  
  datosqa <- list(
    calidad_ciclo = rbind(qa_R1$calidad_ciclo, qa_R2$calidad_ciclo),
    estadisticas = rbind(qa_R1$estadisticas, qa_R2$estadisticas),
    gc_dist = rbind(qa_R1$gc_dist, qa_R2$gc_dist),
    composicion_ciclo = rbind(qa_R1$composicion_ciclo, qa_R2$composicion_ciclo),
    longitud_dist = rbind(qa_R1$longitud_dist, qa_R2$longitud_dist),
    n_muestras = length(muestras$R1),
    n_archivos = qa_R1$n_archivos + qa_R2$n_archivos
  )
  
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  directorio_salida <- file.path("OUTPUT", estudio_id, "REPORT",
                                 paste0("InformeCalidad_", timestamp))
  dir.create(directorio_salida, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(datosqa, file.path(directorio_salida, "datos_qa.rds"))
  
  cat("Creando informe Rmarkdown\n")
  archivoRmd <- generarRmd(directorio_salida, timestamp, umbral_calidad, umbral_ee)
  rmarkdown::render(archivoRmd, quiet = TRUE, output_file = "informe.html")
  
  ruta <- file.path(directorio_salida, "informe.html")
  if (abrir_navegador) browseURL(ruta)
  
  cat("Informe disponible en:", ruta, "\n")
  return(directorio_salida)
}


extraerCalidad <- function(archivos, hebra, subsample = 100000, ciclos_composicion = 30) {
  
  require(ShortRead)
  
  calidad_ciclo <- data.frame()
  estadisticas_generales <- data.frame()
  contenido_gc <- data.frame()
  composicion_ciclo <- data.frame()
  longitud_dist <- data.frame()
  
  if (length(archivos) == 0) {
    return(list(calidad_ciclo = calidad_ciclo,
                estadisticas = estadisticas_generales,
                gc_dist = contenido_gc,
                composicion_ciclo = composicion_ciclo,
                longitud_dist = longitud_dist,
                n_archivos = 0))
  }
  
  for (i in seq_along(archivos)) {
    
    archivo <- archivos[i]
    nombre_muestra <- gsub("\\.fastq(\\.gz)?$", "", basename(archivo))
    cat(sprintf(" [%d/%d] (%s) Procesando: %s\n", i, length(archivos), hebra, nombre_muestra))
    
    fq_stream <- FastqStreamer(archivo, n = subsample)
    
    fq <- tryCatch(
      yield(fq_stream),
      error = function(e) {
        warning(sprintf("Error leyendo %s: %s", archivo, e$message))
        return(NULL)
      }
    )
    
    close(fq_stream)
    
    if (is.null(fq) || length(fq) == 0) next
    
    # CALIDAD POR CICLO
    matrix_qa <- as(quality(fq), "matrix")
    resumen_pos <- apply(matrix_qa, 2, function(col) {
      q <- quantile(col, probs = c(0.10, 0.25, 0.50, 0.75), na.rm = TRUE)
      c(media = mean(col, na.rm = TRUE),
        p10 = q[[1]], q25 = q[[2]], mediana = q[[3]], q75 = q[[4]])
    })
    resumen_pos <- t(resumen_pos)
    
    # ERROR ESPERADO (EE)
    error_matrix <- 10^(-matrix_qa / 10)
    ee_acumulado <- t(apply(error_matrix, 1, cumsum))
    resumen_ee <- apply(ee_acumulado, 2, function(col) {
      q <- quantile(col, probs = c(0.50, 0.75, 0.90), na.rm = TRUE)
      c(ee_mediana = q[[1]], ee_q75 = q[[2]], ee_p90 = q[[3]])
    })
    resumen_ee <- t(resumen_ee)
    
    calidad_ciclo <- rbind(calidad_ciclo, data.frame(
      muestra = nombre_muestra,
      hebra = hebra,
      ciclo = seq_len(nrow(resumen_pos)),
      calidad_media = resumen_pos[, "media"],
      q25 = resumen_pos[, "q25"],
      mediana = resumen_pos[, "mediana"],
      ee_mediana = resumen_ee[, "ee_mediana"],
      ee_q75 = resumen_ee[, "ee_q75"],
      ee_p90 = resumen_ee[, "ee_p90"],
      stringsAsFactors = FALSE
    ))
    
    # ESTADÍSTICAS GENERALES
    longitudes <- width(fq)
    tabla_long <- as.data.frame(table(longitudes), stringsAsFactors = FALSE)
    colnames(tabla_long) <- c("longitud", "n_reads")
    tabla_long$longitud <- as.integer(as.character(tabla_long$longitud))
    tabla_long$muestra <- nombre_muestra
    tabla_long$hebra <- hebra
    longitud_dist <- rbind(longitud_dist, tabla_long)
    secuencias <- sread(fq)
    gc <- letterFrequency(secuencias, "GC", as.prob = TRUE)
    n_freq <- letterFrequency(secuencias, "N")
    
    estadisticas_generales <- rbind(estadisticas_generales, data.frame(
      muestra = nombre_muestra,
      hebra = hebra,
      reads_analizados = length(fq),
      calidad_media = mean(matrix_qa, na.rm = TRUE),
      longitud_min = min(longitudes),
      longitud_media = mean(longitudes),
      longitud_max = max(longitudes),
      gc_content = mean(gc) * 100,
      pct_reads_N = 100 * sum(n_freq > 0) / length(fq),
      stringsAsFactors = FALSE
    ))
    
    # DISTRIBUCIÓN GC
    contenido_gc <- rbind(contenido_gc, data.frame(
      muestra = nombre_muestra,
      hebra = hebra,
      porcentaje_gc = as.numeric(gc) * 100,
      stringsAsFactors = FALSE
    ))
    
    # COMPOSICIÓN POR BASE
    limite_comp <- pmin(ciclos_composicion, width(fq))
    subsecuencias <- subseq(secuencias, start = 1L, end = limite_comp)
    mat_comp <- consensusMatrix(subsecuencias)
    mat_comp <- mat_comp[c("A", "C", "G", "T"), , drop = FALSE]
    
    composicion_ciclo <- rbind(composicion_ciclo, data.frame(
      muestra = nombre_muestra,
      hebra = hebra,
      ciclo = rep(seq_len(ncol(mat_comp)), each = 4),
      base = rep(rownames(mat_comp), times = ncol(mat_comp)),
      count = as.vector(mat_comp),
      stringsAsFactors = FALSE
    ))
  }
  
  return(list(calidad_ciclo = calidad_ciclo,
              estadisticas = estadisticas_generales,
              gc_dist = contenido_gc,
              composicion_ciclo = composicion_ciclo,
              longitud_dist = longitud_dist,
              n_archivos = nrow(estadisticas_generales)))
}


generarRmd <- function(directorio_salida, timestamp, umbral_calidad, umbral_ee) {
  
  archivoRmd <- file.path(directorio_salida, "informe.Rmd")
  
  Sys.setlocale("LC_TIME", "Spanish")
  fecha_informe <- format(Sys.time(), "%d de %B de %Y, %H:%M")
  
  contenido <- c(
    '---',
    'title: "Control de Calidad"',
    paste0('date: "', fecha_informe, '"'),
    'output:', 
    '  html_document:',
    '    toc: true',
    '    toc_float: true',
    '    theme: flatly',
    '    code_folding: hide',
    '---',
    '',
    '```{r setup, include = FALSE}',
    'knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,
                          fig.width = 12, fig.height = 7)',
    'library(ggplot2)',
    'library(dplyr)',
    'library(knitr)',
    'library(DT)',
    '',
    'datos <- readRDS("datos_qa.rds")',
    'calidad_ciclo <- datos$calidad_ciclo',
    'estadisticas <- datos$estadisticas',
    'gc_dist <- datos$gc_dist',
    'composicion_ciclo <- datos$composicion_ciclo',
    'longitud_dist <- datos$longitud_dist',
    'n_muestras <- datos$n_muestras',
    'n_archivos <- datos$n_archivos',
    'alpha_normal <- if (n_muestras <= 10) 0.6 else if (n_muestras <= 50) 0.3 else 0.15',
    paste0('umbral <- ', umbral_calidad),
    paste0('umbral_ee <- ', umbral_ee),
    '```',
    '',
    '## Resumen General',
    '',
    'La tabla recoge las estadísticas básicas de cada archivo, ordenados por calidad media ascendente.  ',
    'Los valores se calculan sobre un subconjunto de hasta 100000 reads por archivo.  ',
    'R1 y R2 se evalúan de forma independiente, por lo que cada muestra ocupa dos filas.',
    '',
    '```{r tabla_general}',
    'datatable(', 
    '  estadisticas %>%',
    '    arrange(calidad_media) %>%',
    '    mutate(',
    '      reads_analizados = format(reads_analizados, big.mark = ","),',
    '      calidad_media = round(calidad_media, 1),',
    '      longitud_min = as.integer(longitud_min),',
    '      longitud_media = round(longitud_media, 0),',
    '      longitud_max = as.integer(longitud_max),',
    '      gc_content = round(gc_content, 1),',
    '      pct_reads_N = round(pct_reads_N, 2)',
    '    ),',
    '    colnames = c("Archivo", "Hebra", "Reads Analizados", "Calidad Media",',
    '                  "Longitud Mín", "Longitud Media", "Longitud Máx", "GC (%)", "% con N"),',
    '    options = list(scrollX = TRUE, scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),',
    '    rownames = FALSE',
    '    )',
    '```',
    '',
    '---',
    '',
    '## 1. Calidad por posición y hebra',
    '',
    'Puntuación de calidad Phred por posición en el read.  ',
    'Cada muestra se representa con dos líneas: la mediana en azul y el primer cuartil Q25 en naranja.  ',
    'Las líneas horizontales discontinuas marcan Q30 y Q20 como referencias en control de calidad de los datos.',
    '',
    '```{r calidad_distribucion, fig.height = 6}',
    '',
    'ggplot(calidad_ciclo, aes(x = ciclo, group = muestra)) +',
    '  geom_line(aes(y = mediana), color = "steelblue", alpha = alpha_normal, linewidth = 0.7) +',
    '  geom_line(aes(y = q25), color = "orange", alpha = alpha_normal, linewidth = 0.7) +',
    '  geom_hline(yintercept = 30, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +',
    '  geom_hline(yintercept = umbral, linetype = "dashed", color = "red", linewidth = 0.8) +',
    '  facet_wrap(~ hebra) +',
    '  scale_y_continuous(limits = c(0, 42), breaks = seq(0, 40, 5)) +',
    '  labs(x = "Posición en el read (pb)", y = "Phred Score") +',
    '  theme_minimal(base_size = 13)',
    '```',
    '',
    '---',
    '',
    '## 2. Distribución de longitudes',
    '',
    'Número de reads por longitud de secuencia (pb).',
    '',
    '```{r longitud_distribucion, fig.height = 4, results = "asis"}',
    'longitudes_unicas <- unique(longitud_dist$longitud)',
    'if (length(longitudes_unicas) == 1) {',
    '  cat(sprintf(',
    '    "Todas las secuencias presentan longitud uniforme de %d pb.", longitudes_unicas))',
    '} else {',
    '  print(',
    '    ggplot(longitud_dist, aes(x = longitud, y = n_reads, group = muestra)) +',
    '      geom_line(alpha = alpha_normal, color = "steelblue", linewidth = 0.7) +',
    '      geom_point(alpha = alpha_normal, color = "steelblue", size = 1.5) +',
    '      facet_wrap(~ hebra) +',
    '      labs(x = "Longitud (pb)", y = "Número de reads") +',
    '      theme_minimal(base_size = 13)',
    '  )',
    '}',
    '```',
    '',
    '---',
    '',
    '## 3. Error Esperado acumulado (EE)',
    '',
    'Error esperado (EE) acumulado a lo largo de la posición en el read.  ',
    'Para cada muestra se trazan tres curvas, la mediana en azul, el percentil 75 en naranja y el percentil 90 en morado.  ',
    'La línea roja discontinua señala el valor de maxEE configurado como umbral de referencia.',
    '',
    '```{r ee_acumulado, fig.height = 6}',
    'ggplot(calidad_ciclo, aes(x = ciclo, group = muestra)) +',
    '  geom_line(aes(y = ee_mediana), color = "steelblue", alpha = alpha_normal, linewidth = 0.7) +',
    '  geom_line(aes(y = ee_q75), color = "orange", alpha = alpha_normal, linewidth = 0.7) +',
    '  geom_line(aes(y = ee_p90), color = "purple4", alpha = alpha_normal, linewidth = 0.7) +',
    '  geom_hline(yintercept = umbral_ee, linetype = "dashed", color = "red", linewidth = 0.8) +',
    '  facet_wrap(~ hebra) +',
    '  labs(x = "Posición en el read (pb)", y = "Error esperado acumulado") +',
    '  theme_minimal(base_size = 13)',
    '```',
    '',
    '---',
    '',
    '## 4. Composición de bases por posición y hebra',
    '',
    'Porcentaje de cada base (A, C, G, T) por posición en los primeros 30 ciclos, agrupando todas las muestras.',
    '',
    '```{r composicion_bases, fig.height = 5}',
    'composicion_resumen <- composicion_ciclo %>%',
    '  group_by(hebra, ciclo, base) %>%',
    '  summarise(total = sum(count), .groups = "drop") %>%',
    '  group_by(hebra, ciclo) %>%',
    '  mutate(porcentaje = 100 * total / sum(total)) %>%',
    '  ungroup()',
    '',
    'ggplot(composicion_resumen, aes(x = ciclo, y = porcentaje, color = base)) +',
    '  geom_line(linewidth = 1) +',
    '  scale_color_manual(values = c(A = "darkgreen", C = "blue", G = "black", T = "red")) +',
    '  facet_wrap(~ hebra) +',
    '  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +',
    '  scale_x_continuous(breaks = seq(0, 30, 5)) +',
    '  labs(x = "Posición en el read (pb)", y = "Composición (%)", color = "Base") +',
    '  theme_minimal(base_size = 13)',
    '```',
    '',
    '---',
    '',
    '## 5. Archivos por debajo del umbral de calidad',
    '',
    'Comprobación del umbral de calidad media por archivo.  ', 
    'Si se detectan archivos por debajo del umbral mínimo de calidad, se muestra la tabla de estadísticas y el perfil de calidad por posición de los 10 primeros.',
    '',
    '```{r muestras_bajas, results = "asis"}',
    'muestras_bajas <- estadisticas %>%',
    '  filter(calidad_media < umbral) %>%',
    '  arrange(calidad_media)',
    '',
    'if (nrow(muestras_bajas) > 0) {',
    '  cat("\\nSe identificaron", nrow(muestras_bajas), "archivo(s) con calidad media inferior a Q", umbral, ":\\n\\n")',
    '  cat("*R1 y R2 se evalúan por separado; una misma muestra biológica puede aparecer hasta dos veces si ambas hebras quedan por debajo del umbral.*\\n\\n")',
    '  print(kable(',
    '    muestras_bajas %>%',
    '      mutate(reads_analizados    = format(reads_analizados, big.mark = ","),',
    '             calidad_media  = round(calidad_media, 1),',
    '             longitud_min   = as.integer(longitud_min),',
    '             longitud_media = round(longitud_media, 0),',
    '             longitud_max   = as.integer(longitud_max),',
    '             gc_content     = round(gc_content, 1),',
    '             pct_reads_N    = round(pct_reads_N, 2)),',
    '    col.names = c("Archivo", "Hebra", "Reads Analizados", "Cal. Media",',
    '                  "Long. Mín", "Long. Media", "Long. Máx", "GC (%)", "% con N"),',
    '    caption = paste0("Archivos con calidad media < Q", umbral)',
    '  ))',  
    '} else {',
    '  cat("\\nTodos los archivos presentan calidad media ≥ Q", umbral, ".\\n")',
    '}',
    '```',
    '',
    '```{r graficos_problematicas, fig.height = 4, results = "asis"}',
    'if (nrow(muestras_bajas) > 0) {',
    '  cat("\\n### Perfil por archivo\\n\\n")',
    '  muestras_mostrar <- head(muestras_bajas$muestra, 10)',
    '  for (m in muestras_mostrar) {',
    '    cat("\\n#### Archivo:", m, "\\n\\n")',    
    '    datos_muestra <- filter(calidad_ciclo, muestra == m)',
    '    p_individual <- ggplot(datos_muestra, aes(x = ciclo, y = calidad_media)) +',
    '      geom_line(color = "darkred", linewidth = 1.2) +',
    '      geom_ribbon(aes(ymin = 0, ymax = calidad_media), fill = "darkred", alpha = 0.2) +',
    '      geom_hline(yintercept = 30, linetype = "dashed", color = "darkgreen") +',
    '      geom_hline(yintercept = umbral, linetype = "dashed", color = "red") +',
    '      scale_y_continuous(limits = c(0, 42)) +',     
    '      labs(x = "Posición (pb)", y = "Phred Score") +',
    '      theme_minimal(base_size = 12)',
    '   print(p_individual)',
    '   cat("\\n")',
    '  }',
    '  if (nrow(muestras_bajas) > 10)',
    '    cat("\\n*Se muestran los 10 archivos con menor calidad media.*\\n")',
    '}',
    '```',
    '',
    '---',
    '',
    '## 6. Distribución de calidad media',
    '',
    'Distribución de la calidad media Phred por archivo, diferenciando R1 y R2.',
    '',
    '```{r distribucion_calidad, fig.height = 5}',
    'ggplot(estadisticas, aes(x = calidad_media, fill = hebra)) +',
    '  geom_histogram(bins = 30, color = "white", alpha = 0.6, position = "identity") +',
    '  scale_fill_manual(values = c(R1 = "steelblue", R2 = "coral")) +',
    '  geom_vline(xintercept = umbral, linetype = "dashed", color = "red", linewidth = 1) +',
    '  labs(x = "Calidad Media (Phred Score)", y = "Número de Archivos", fill = "Hebra") +',
    '  theme_minimal(base_size = 13)',
    '```',
    '',
    '---',
    '',
    '## 7. Contenido GC',
    '',
    'Distribución de densidad del contenido GC por lectura, agrupando R1 y R2.  ',
    'La línea discontinua marca el 50% como referencia visual.',
    '',
    '```{r gc_content, fig.height = 5}',
    'ggplot(gc_dist, aes(x = porcentaje_gc)) +',
    '  geom_density(fill = "coral", alpha = 0.6, color = "coral4", linewidth = 1) +',
    '  geom_vline(xintercept = 50, linetype = "dotted", color = "gray30", linewidth = 1) +',
    '  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +',
    '  labs(x = "Contenido GC (%)", y = "Densidad") +',
    '  theme_minimal(base_size = 13)',
    '```',
    '',
    '---',
    '',
    '## 8. Conclusiones',
    '',
    '```{r conclusiones, results="asis"}',
    'reads_analizados <- sum(estadisticas$reads_analizados, na.rm = TRUE)',
    'calidad_global <- mean(estadisticas$calidad_media, na.rm = TRUE)',
    'n_problematicas <- sum(estadisticas$calidad_media < umbral, na.rm = TRUE)',
    'porc_problematicas <- round(100 * n_problematicas / n_archivos, 1)',
    '',
    'cat("Total de reads procesados:", format(reads_analizados, big.mark = ","), "\\n\\n")',
    'cat("Archivos evaluados:", n_archivos, "\\n\\n")',
    'cat("Archivos con calidad media < Q", umbral, ":", n_problematicas, "de",',
    '    n_archivos, paste0("(", porc_problematicas, "%)"), "\\n\\n")',
    '',
    'if (n_problematicas > 0) {',
    '  cat(\'<div class="alert alert-info">Los archivos listados en el apartado 5',
    '      presentan calidad media por debajo del umbral configurado.  ',
    '      Se recomienda revisarlos antes de proceder con el filtrado.</div>\\n\')',
    '}',
    '```',
    '',
    '---',
    '',
    paste0('*Informe generado el ', fecha_informe, '*')
  )
  
  writeLines(contenido, archivoRmd)
  return(archivoRmd)
}