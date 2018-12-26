
analysisDirs <- c(
  "SCM.analysis", "visualizingPoissonAlgo", "visualizingSolutionSpace",
  "syntheticMultipletsFromCounts"
)

for(i in 1:length(analysisDirs)) {
  setwd(analysisDirs[i])
  source('./scripts/script.R')
  rmarkdown::render('./analysis/analysis.Rmd')
}
