#PACKAGES
packages <- c("CIMseq", "CIMseq.data", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

##DATA
if(file.exists('../MGA.analysis_SI/data/CIMseqData.rda')) {
  load('../MGA.analysis_SI/data/CIMseqData.rda')
}

idx <- sample(1:ncol(getData(cObjMul, "counts")), 8, FALSE)
cObjMul.small <-  CIMseqMultiplets(
  getData(cObjMul, "counts")[, idx], 
  getData(cObjMul, "counts.ercc")[, idx], 
  getData(cObjMul, "features")
)

##TEST
nSwarmMembers <- c(100, 500, 1000, 5000)
res <- map(nSwarmMembers, function(n) {
  plan(multiprocess)
  CIMseqSwarm(
    cObjSng, cObjMul, maxiter = 100, swarmsize = n, nSyntheticMultiplets = 400
  )
})

save(res, file = file.path(currPath, "data/results.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))