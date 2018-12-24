#SCM.analysis
setwd('~/Github/CIMseq.testing/analysis/SCM.analysis')
source('./scripts/runCountsSorted2.R')
rmarkdown::render('./analysis/analysis.Rmd')

#visualizingPoissonAlgo
setwd('~/Github/CIMseq.testing/analysis/visualizingPoissonAlgo')
source('./scripts/script.R')
rmarkdown::render('./analysis/analysis.Rmd')

#visualizingSolutionSpace
setwd('~/Github/CIMseq.testing/analysis/visualizingSolutionSpace')
source('./scripts/script.R')
rmarkdown::render('./analysis/version2.Rmd')
