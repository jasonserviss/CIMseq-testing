library(here)
baseDir <- 'analysis'
analysisDirs <- c(
  "SCM.analysis", "visualizingPoissonAlgo", "visualizingSolutionSpace",
  "syntheticMultipletsFromCounts"
)
directories <- system.file(file.path(baseDir, analysisDirs), package = 'CIMseq.testing')

for(i in 1:length(directories)) {
  source(file.path(directories[i], 'scripts/script.R'))
  rmarkdown::render(file.path(directories[i], 'analysis/analysis.Rmd'))
}
