graficosCalidad <- function(R1,R2){
  timestamp <- format(Sys.time(), "%H_%M_%S")
  
  grafR1 <- plotQualityProfile(R1[1:10])
  ggsave(sprintf("OUTPUT/FIGURES/graficoCalidad_R1_%s.png", timestamp), grafR1, width = 10, height = 7)
  
  grafR2 <- plotQualityProfile(R2[1:10])
  ggsave(sprintf("OUTPUT/FIGURES/graficoCalidad_R2_%s.png", timestamp), grafR2, width = 10, height = 7)
}