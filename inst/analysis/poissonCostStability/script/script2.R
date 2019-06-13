library(CIMseq)
library(future.apply)
library(tidyverse)
load('../MGA.analysis_SI/data/CIMseqData.rda')
options(future.globals.maxSize = Inf)
currPath <- getwd()

set.seed(344905)
idx <- sample(1:ncol(getData(cObjMul, "counts")), 40, FALSE)
cObjMul.small <-  CIMseqMultiplets(
  getData(cObjMul, "counts")[, idx], 
  getData(cObjMul, "counts.ercc")[, idx], 
  getData(cObjMul, "features")
)

nSyntheticMultiplets <- c(10, 50, 100, 500, 1000, 5000, 10000)
res <- map(nSyntheticMultiplets, function(n) {
  init <- swarmInit(cObjSng, 2)
  future::plan(future::multiprocess, workers = min(8, length(idx)))
  sObj <- CIMseqSwarm(
    cObjSng, cObjMul.small, maxiter = 100, swarmsize = ncol(init), 
    nSyntheticMultiplets = n, swarmInit = init
  )
})

costs <- map(res, function(r) getData(r, "costs"))
names(costs) <- nSyntheticMultiplets
save(costs, file = file.path(currPath, "data/costs2.rda"))

