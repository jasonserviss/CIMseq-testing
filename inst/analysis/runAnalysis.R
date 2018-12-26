baseDirs <- '~/Github/CIMseq.testing/inst/analysis'
analysisDirs <- c(
  "SCM.analysis", "visualizingPoissonAlgo", "visualizingSolutionSpace",
  "syntheticMultipletsFromCounts"
)
fullPaths <- file.path(baseDirs, analysisDirs)
for(i in 1:length(analysisDirs)) {
  setwd(analysisDirs[i])
  source('./scripts/script.R')
  rmarkdown::render('./analysis/analysis.Rmd')
}
