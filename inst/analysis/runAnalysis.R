
directories <- file.path('~/Github/CIMseq.testing/inst/analysis', c(
  "SCM.analysis", "visualizingPoissonAlgo", "visualizingSolutionSpace",
  "syntheticMultipletsFromCounts"
))

for(i in 1:length(directories)) {
  print(paste0("Processing ", basename(directories[i])))
  setwd(directories[i])
  print(paste0("Working dir ", getwd()))
  print(list.files(directories[i], recursive = TRUE))
  if(file.exists('./scripts/script.R')) {
    print(paste0("Running ", file.path(directories[i], 'scripts/script.R')))
    source(file.path(directories[i], 'scripts/script.R'))
  }
  rmarkdown::render(file.path(directories[i], 'analysis/analysis.Rmd'))
  setwd('~/Github/CIMseq.testing')
}
