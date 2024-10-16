
# AUTOMATIZACIÓN DE LA DESCARGA DE LAS MUESTRAS DEL ENA

library(readr)
library(dplyr)
library(Biostrings)
library(svDialogs)

#   DESCARGA DEL TSV 

#   DESCARGA DE SECUENCIAS
directorio <- dlgInput(message = "Introduzca el directorio al tsv: ")$res

tabla <- read_delim(directorio, col_names = TRUE, delim ="\t")

tabla_ordenada <- tabla[order(tabla$submitted_ftp), ]

enlaces <- strsplit(tabla_ordenada$submitted_ftp, ";")

for (enl in enlaces){
  for (i in enl){
    destino <- strsplit(i, "/")[[1]]
    download.file(i,destino[length(destino)])
  }
}