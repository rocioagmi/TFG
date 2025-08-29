graficosCalidad <- function(R1,R2){
  plotQualityProfile(R1[1:10])
  ggsave("OUTPUT/FIGURES/graficoCalidad_")
  plotQualityProfile(R2[1:10])
}