#ARGS
args <- commandArgs(TRUE)

#PACKAGES
packages <- c("tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
if(file.exists(file.path('../MGA.analysis_engeSimple', 'data/CIMseqData.rda'))) {
  load(file.path('../MGA.analysis_engeSimple', 'data/CIMseqData.rda'))
}

if(!is.na(args[1])) {
  library(CIMseq, lib.loc = "/home/jason/R/x86_64-redhat-linux-gnu-library/3.5")
  library(CIMseq.data, lib.loc = "/home/jason/R/x86_64-redhat-linux-gnu-library/3.5")
  library(CIMseq.testing, lib.loc = "/home/jason/R/x86_64-redhat-linux-gnu-library/3.5")
  out <- if(is.na(args[2])) {out <- 'tmp'} else {out <- as.character(args[2])}

  options(future.globals.maxSize = Inf)  
  runSwarmUppmax(
    cObjSng, cObjMul, swarmInit = NULL, maxiter = 100,
    nSyntheticMultiplets = 5000, eps.stagnate = 1, maxit.stagnate = 5,
    currPath = out, args = args
  )
  
} else {
  library(CIMseq)
  library(CIMseq.data)
  library(CIMseq.testing)
  
  runSwarmMultiprocess(
    cObjSng, cObjMul, swarmInit = NULL, maxiter = 100,
    nSyntheticMultiplets = 400, eps.stagnate = 1, maxit.stagnate = 5,
    currPath = getwd()
  )
}

writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_CIMseqSwarm.txt"))
print("finished")
