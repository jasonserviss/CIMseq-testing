library(CIMseq)
library(future.apply)
library(tidyverse)
load('../MGA.analysis_all/data/CIMseqData.rda')
currPath <- getwd()

idx <- sample(1:ncol(getData(cObjMul, "counts")), 21, FALSE)
cObjMul.small <-  CIMseqMultiplets(
  getData(cObjMul, "counts")[, idx], 
  getData(cObjMul, "counts.ercc")[, idx], 
  getData(cObjMul, "features")
)

nSyntheticMultiplets <- c(10, 100, 400, 800, 1200, 1600, 2000)
res <- map(nSyntheticMultiplets, function(n) {
  future::plan(multiprocess)
  sObj <- CIMseqSwarm(
    cObjSng, cObjMul.small, maxiter = 10, swarmsize = 50, nSyntheticMultiplets = n
  )
})

costs <- map(res, function(r) getData(r, "costs"))
save(costs, file = file.path(currPath, "data/costs2.rda"))

