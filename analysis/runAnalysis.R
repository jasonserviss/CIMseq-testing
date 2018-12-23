#SCM.analysis
setwd('/home/Github/CIMseq.testing/analysis/SCM.analysis')
source('./scripts/runCountsSorted2.R')
rmarkdown::render('./analysis/analysis.Rmd')

