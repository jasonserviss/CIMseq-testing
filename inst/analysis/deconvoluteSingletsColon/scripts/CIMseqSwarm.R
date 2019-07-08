#ARGS
args <- commandArgs(TRUE)

#PACKAGES
packages <- c("tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
target <- '../MGA.analysis_SI20/data/CIMseqData.rda'
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

  #gives 620 swarm members
  baseseed <- 72839
  init <- map(1:2, function(i) {
    cbind(
      swarmInit(cObjSng, 2, null.weight = 1, seed = baseseed + i),
      swarmInit(cObjSng, 3, null.weight = 1, seed = baseseed + i)
    )
  }) %>%
    do.call(cbind, .)
  
  runSwarmUppmax(
    cObjSng, cObjMul.2, swarmInit = init, maxiter = 100,
    nSyntheticMultiplets = 2000, eps.stagnate = 1, maxit.stagnate = 5,
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

  #gives 620 swarm members
  baseseed <- 72839
  init <- map(1:2, function(i) {
    cbind(
      swarmInit(cObjSng, 2, null.weight = 1, seed = baseseed + i),
      swarmInit(cObjSng, 3, null.weight = 1, seed = baseseed + i)
    )
  }) %>%
    do.call(cbind, .)

  runSwarmMultiprocess(
    cObjSng, cObjMul.2, swarmInit = init, maxiter = 100,
    nSyntheticMultiplets = 2000, eps.stagnate = 1, maxit.stagnate = 5,
    currPath = getwd()
  )
}

writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_CIMseqSwarm.txt"))
print("finished")
