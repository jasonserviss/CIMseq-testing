
runAnalysis <- function(rootPath) {
  directories <- file.path(rootPath, 'inst/analysis', c(
    "SCM.analysis", "visualizingPoissonAlgo", "visualizingSolutionSpace",
    "syntheticMultipletsFromCounts", "MGA.analysis_SI", "poissonCostStability",
    "singletVsMultipletCellFreq"
  ))
  
  for(i in 1:length(directories)) {
    print(paste0("Processing ", basename(directories[i])))
    setwd(directories[i])
    print(paste0("Working dir ", getwd()))
    print(list.files(directories[i], recursive = TRUE))
    scripts <- list.files('scripts', full.names = TRUE)
    for(y in 1:length(scripts)){
      print(paste0("Running ", scripts[y]))
      source(scripts[y])
    }
    rmarkdown::render(file.path(directories[i], 'analysis/analysis.Rmd'))
  }
}

