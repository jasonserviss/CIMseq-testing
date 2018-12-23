#SCM.analysis
setwd('~/Github/CIMseq.testing/analysis/SCM.analysis')
source('./scripts/runCountsSorted2.R')
rmarkdown::render('./analysis/analysis.Rmd')

