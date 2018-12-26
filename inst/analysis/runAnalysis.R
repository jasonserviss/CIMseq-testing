baseDirs <- './inst/analysis'
analysisDirs <- c(
  "SCM.analysis", "visualizingPoissonAlgo", "visualizingSolutionSpace",
  "syntheticMultipletsFromCounts"
)
directories <- file.path(baseDirs, analysisDirs)
for(i in 1:length(directories)) {
  setwd(directories[i])
  source('./scripts/script.R')
  rmarkdown::render('./analysis/analysis.Rmd')
}
