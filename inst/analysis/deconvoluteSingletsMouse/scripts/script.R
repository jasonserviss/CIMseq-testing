#PACKAGES
packages <- c("CIMseq", "CIMseq.data", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
target <- '../MGA.analysis_enge/data/CIMseqData.rda'
if(file.exists(target)) {
  load(target)
} else {
  stop("Could not load CIMseq data")
}
rm(cObjMul)
cObjMul <- cObjSng

##deconvolution
print(paste0("Starting deconvolution at ", Sys.time()))
future::plan(multiprocess)
sObj <- CIMseqSwarm(
  cObjSng, cObjMul, maxiter = 10, swarmsize = 150, nSyntheticMultiplets = 400
)
#maxiter = 100, swarmsize = 500
print(paste0("Finished deconvolution at ", Sys.time()))

save(sObj, file = file.path(currPath, "data/sObj.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spSwarm.txt"))
print("finished")
