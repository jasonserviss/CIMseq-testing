
baseDir <- '~/Github/CIMseq.testing/inst/analysis'
analysisDirs <- c(
  "SCM.analysis", "visualizingPoissonAlgo", "visualizingSolutionSpace",
  "syntheticMultipletsFromCounts"
)
directories <- file.path(baseDir, analysisDirs)

for(i in 1:length(directories)) {
  setwd(directories[i])
  source(file.path(directories[i], 'scripts/script.R'))
  rmarkdown::render(file.path(directories[i], 'analysis/analysis.Rmd'))
  setwd('~/Github/CIMseq.testing')
}
