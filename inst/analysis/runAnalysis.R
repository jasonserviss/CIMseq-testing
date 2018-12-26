
baseDir <- '~/Github/CIMseq.testing/inst/analysis'
analysisDirs <- c(
  "SCM.analysis", "visualizingPoissonAlgo", "visualizingSolutionSpace",
  "syntheticMultipletsFromCounts"
)
directories <- file.path(baseDir, analysisDirs)

for(i in 1:length(directories)) {
  print(paste0("Processing ", basename(directories[i])))
  setwd(directories[i])
  print(paste0("Working dir ", getwd()))
  print(list.files(directories[i]))
  if(file.exists('./scripts/script.R')) source(file.path(directories[i], 'scripts/script.R'))
  rmarkdown::render(file.path(directories[i], 'analysis/analysis.Rmd'))
  setwd('~/Github/CIMseq.testing')
}
