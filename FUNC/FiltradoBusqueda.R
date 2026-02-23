library(svDialogs)
library(stringr)

filtrarBusqueda <- function(df) {
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df_filtrado <- df
  
  cat("INICIANDO FILTRADO\n")
  cat(sprintf("Total de muestras inicial: %d\n\n", nrow(df_filtrado)))
  
  campos_texto <- c("scientific_name", "experiment_title", "study_title",
                    "sample_title", "run_alias", "sample_alias")
  campos_disponibles <- campos_texto[campos_texto %in% colnames(df_filtrado)]
  
  campos_con_datos <- campos_disponibles[sapply(campos_disponibles, function(col) {
    any(!is.na(df_filtrado[[col]]) & df_filtrado[[col]] != "")
  })]
  
  if (length(campos_con_datos) == 0) {
    cat("No hay campos de texto con datos para filtrar.\n\n")
    return(df_filtrado)
  }
  
  metodo <- dlg_list(
    choices = c(
      "Usar términos sugeridos del diccionario",
      "Introducir mis propios términos"),
    title = "¿Cómo deseas filtrar?",
    multiple = FALSE)$res
  
  if (length(metodo) == 0) {
    cat("No se seleccionó método. Se devuelven las muestras sin filtrar.\n")
    return(df_filtrado)
  }
  
  # TÉRMINOS SUGERIDOS DEL DICCIONARIO - TERMINAR DE MIRAR ESTO
  # filtro automático que garantiza muestras humanas, intestinales o fecales de microbioma y/o 16S
  
  if (metodo == "Usar términos sugeridos del diccionario") {
    
    filtros_obligatorios <- list(
      "Humano" = "\\b(human|homo\\s*sapiens)\\b",
      "Intestinal/Fecal" = "\\b(fecal|stool|feces|faecal|gut|intestin(al|e)|colon|rectal|bowel)\\b")
      "Microbioma/Metagenómica" = "\\b(microbiom(e|a)|metagenomic(s|a)?|16S(\\s*(rRNA|RNA|ribosomal))?)\\b"
    
    aplica_patron <- function(df_local, patron, campos) {
      matches <- rep(FALSE, nrow(df_local))
      for (col in campos) {
        valores <- df_local[[col]]
        valores[is.na(valores)] <- ""
        matches <- matches | grepl(patron, valores, ignore.case = TRUE, perl = TRUE)
      }
      matches
    }
    
    for (nombre_grupo in names(filtros_obligatorios)) {
      patron <- filtros_obligatorios[[nombre_grupo]]
      mascara <- aplica_patron(df_filtrado, patron, campos_con_datos)
      
      n_antes  <- nrow(df_filtrado)
      df_filtrado <- df_filtrado[mascara, , drop = FALSE]
      n_despues <- nrow(df_filtrado)
      
      cat(sprintf("  [%s] %d -> %d muestras (%.1f%% superaron)\n",
                  nombre_grupo, n_antes, n_despues,
                  ifelse(n_antes > 0, (n_despues / n_antes) * 100, 0)))
      
      if (nrow(df_filtrado) == 0) {
        cat("\nNo quedan muestras después de aplicar el filtrado.\n")
        break
      }
    }
    
    
    # INTRODUCIR MIS PROPIOS TÉRMINOS
    # El usuario escribe uno o varios términos separados por comas y se buscan en TODOS los campos con datos.
    
  } else if (metodo == "Introducir mis propios términos") {
    
    terminos <- dlgInput(
      message = "Introduce términos separados por comas:", default = "")$res
    
    if (length(terminos) == 0 || nchar(trimws(terminos)) == 0) {
      cat("No se introdujeron términos. Se devuelven las muestras sin filtrar.\n")
      return(df_filtrado)
    }
    
    terminos_usuario <- trimws(unlist(strsplit(terminos, ",")))
    terminos_usuario <- terminos_usuario[nchar(terminos_usuario) > 0]
    
    patron <- paste0("\\b(", paste(terminos_usuario, collapse = "|"), ")\\b")
    
    cat(sprintf("\nBuscando términos: %s\n", paste(terminos_usuario, collapse = ", ")))
    cat(sprintf("Campos donde se busca: %s\n", paste(campos_con_datos, collapse = ", ")))
    
    df_filtrado <- df_filtrado %>%
      filter(if_any(all_of(campos_con_datos),
                    ~grepl(patron, ., ignore.case = TRUE, perl = TRUE)))
    
    cat(sprintf("Muestras tras filtro: %d\n", nrow(df_filtrado)))
  } 
  
  # ===========================================================
  # RESUMEN FINAL
  # ===========================================================
  cat(sprintf("\n ---- FILTRADO COMPLETADO ----\n"))
  cat(sprintf("Muestras iniciales : %d\n", nrow(df)))
  cat(sprintf("Muestras finales   : %d\n", nrow(df_filtrado)))
  cat(sprintf("Muestras eliminadas: %d (%.1f%%)\n",
              nrow(df) - nrow(df_filtrado),
              ifelse(nrow(df) > 0, ((nrow(df) - nrow(df_filtrado)) / nrow(df)) * 100, 0)))
  
  return(df_filtrado)
}
