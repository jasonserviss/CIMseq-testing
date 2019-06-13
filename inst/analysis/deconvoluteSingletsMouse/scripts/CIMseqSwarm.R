#ARGS
args <- commandArgs(TRUE)

#PACKAGES
packages <- c("tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
target <- '../MGA.analysis_enge/data/CIMseqData.rda'
if(file.exists(target)) {
  load(target)
}

if(!is.na(args[1])) {
  library(CIMseq, lib.loc = "/home/jason/R/x86_64-redhat-linux-gnu-library/3.5")
  library(CIMseq.data, lib.loc = "/home/jason/R/x86_64-redhat-linux-gnu-library/3.5")
  library(CIMseq.testing, lib.loc = "/home/jason/R/x86_64-redhat-linux-gnu-library/3.5")
  out <- if(is.na(args[2])) {out <- 'tmp'} else {out <- as.character(args[2])}

  cObjMul.2 <- CIMseqMultiplets(
    getData(cObjSng, "counts"),
    getData(cObjSng, "counts.ercc"),
    getData(cObjMul, "features")
  )

  #gives 2070 swarm members
  init <- cbind(
    swarmInit(cObjSng, 2, null.weight = 1, seed = 35466),
    swarmInit(cObjSng, 3, null.weight = 1, seed = 35466)
  )

  
  runSwarmUppmax(
    cObjSng, cObjMul.2, swarmInit = init, maxiter = 100,
    nSyntheticMultiplets = 400, eps.stagnate = 1, maxit.stagnate = 5,
    currPath = out, args = args
  )
  
} else {
  library(CIMseq)
  library(CIMseq.data)
  library(CIMseq.testing)
  
  cObjMul.2 <- CIMseqMultiplets(
    getData(cObjSng, "counts"),
    getData(cObjSng, "counts.ercc"),
    getData(cObjMul, "features")
  )

  #gives 2070 swarm members
  init <- cbind(
    swarmInit(cObjSng, 2, null.weight = 1, seed = 35466),
    swarmInit(cObjSng, 3, null.weight = 1, seed = 35466)
  )

  runSwarmMultiprocess(
    cObjSng, cObjMul.2, swarmInit = init, maxiter = 100,
    nSyntheticMultiplets = 400, eps.stagnate = 1, maxit.stagnate = 5,
    currPath = getwd()
  )
}

writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_CIMseqSwarm.txt"))
print("finished")
