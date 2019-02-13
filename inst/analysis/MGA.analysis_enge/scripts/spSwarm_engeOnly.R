#PACKAGES
packages <- c("CIMseq", "CIMseq.data", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
if(file.exists(file.path(currPath, 'data/CIMseqData.rda'))) {
  load(file.path(currPath, 'data/CIMseqData.rda'))
}

# Testing
# cObjMul <- CIMseqMultiplets(
#   getData(cObjMul, "counts")[, 1:2],
#   getData(cObjMul, "counts.ercc")[, 1:2],
#   getData(cObjMul, "features")
# )

#future::plan(multiprocess)
options(future.wait.interval = 10000.0)
options(future.wait.timeout = 1e9)
future::plan(
  future.batchtools::batchtools_slurm,
  template = "/crex/proj/snic2018-8-151/private/batchtools.slurm.tmpl",
  resources = list(
    account = "snic2018-8-151", partition = "core", ntasks = 1L,
    time = "24:00:00", jobname = "MGA_enge",
    modules = "R_packages/3.5.2", R = "R/3.5.2", log.file = file.path(currPath, "logs/slurm.txt")
   ),
  workers = 100
)

print(paste0("Starting deconvolution at ", Sys.time()))
sObj <- CIMseqSwarm(
  cObjSng, cObjMul, maxiter = 100, swarmsize = 500, nSyntheticMultiplets = 400
)
print(paste0("Finished deconvolution at ", Sys.time()))

save(sObj, file = file.path(currPath, "data/sObj.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spSwarm.txt"))
print("finished")
