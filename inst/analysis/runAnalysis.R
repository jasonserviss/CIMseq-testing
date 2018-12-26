analysisDirs <- list.dirs(
  '~/Github/CIMseq.testing/inst/analysis', recursive = FALSE
)

for(i in 1:length(analysisDirs)) {
  setwd(analysisDirs[i])
  source('./scripts/script.R')
  rmarkdown::render('./analysis/analysis.Rmd')
}
