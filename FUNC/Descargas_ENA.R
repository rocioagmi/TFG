library(readr)
library(httr)

descargas_ENA <- function(nAcceso){

  # DESCARGA DEL TSV
  enlaceTSV <- paste("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=","&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_bytes,fastq_ftp,submitted_bytes,submitted_ftp,sra_bytes,sra_ftp,bam_ftp,bam_bytes&format=tsv&download=true&limit=0", sep = nAcceso)
  directorio <- mkdir("INPUT/DATA")
  directorioTSV <- paste("INPUT/DATA/filereport_read_run_","_tsv.tsv",sep = nAcceso)
  tsv <- download.file(enlaceTSV, directorioTSV)

  # DESCARGA DE LAS MUESTRAS
  tabla <- read_delim(directorioTSV, col_names = TRUE, delim ="\t")
  tabla_ordenada <- tabla[order(tabla$submitted_ftp), ]
  enlaces <- strsplit(tabla_ordenada$submitted_ftp, ";")

  for (enl in enlaces){
    for (i in enl){
      destino <- strsplit(i, "/")[[1]]
      directorioEnl <- paste("INPUT/DATA/",destino[length(destino)])
      download.file(i,directorioEnl)
    }
  }
}