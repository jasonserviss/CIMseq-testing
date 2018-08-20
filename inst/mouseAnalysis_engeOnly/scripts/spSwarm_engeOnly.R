#PACKAGES
packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

s <- str_detect(colnames(countsMgfp), "^s")
e <- colnames(countsMgfp) %in% filter(countsMgfpMeta, is.na(GFP) & !filtered)$sample
boolSng <- s & e
boolMul <- !s & e

#setup spCounts
cObjSng <- spCounts(countsMgfp[, boolSng], countsMgfpERCC[, boolSng])
cObjMul <- spCounts(countsMgfp[, boolMul], countsMgfpERCC[, boolMul])

print("spCounts done")

#spUnsupervised
if(file.exists(file.path(currPath, 'data/uObj_engeOnly.rda'))) {
  load(file.path(currPath, 'data/uObj_engeOnly.rda'))
}

print("spUnsupervised done")

##spSwarm
selectIdx <- spTopVar(cObjSng, 2000)

#future::plan(multiprocess)
options(future.wait.interval = 10.0)
future::plan(
  future.batchtools::batchtools_slurm,
  template = "/crex/proj/snic2018-8-151/private/batchtools.slurm.tmpl",
  resources = list(
    account = "snic2018-8-151", partition = "core", ntasks = 1L,
    time = "24:00:00", jobname = "testingPoissonSorted",
    modules = "R_packages/3.5.0", R = "R/3.5.0", log.file = file.path(currPath, "logs/slurm.txt")
   ),
  workers = 100
)

print(paste0("Starting deconvolution at ", Sys.time()))
sObj <- spSwarm(
  cObjSng, cObjMul, uObj, maxiter = 100, swarmsize = 500,
  nSyntheticMultiplets = 400, selectInd = selectIdx
)
print(paste0("Finished deconvolution at ", Sys.time()))

save(sObj, file = file.path(currPath, "data/sObj_engeOnly.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.spSwarm_engeOnly.txt"))
print("finished")
