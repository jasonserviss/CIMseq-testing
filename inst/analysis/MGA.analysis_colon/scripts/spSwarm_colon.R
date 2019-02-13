#PACKAGES
packages <- c("CIMseq", "CIMseq.data", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
if(file.exists(file.path(currPath, 'data/CIMseqData.rda'))) {
  load(file.path(currPath, 'data/CIMseqData.rda'))
}

#TESTING
set.seed(9348)
sections <- MGAS.Meta %>% filter(!filtered & Section != "Undefined") %>% group_by(Section) %>% sample_n(20) %>% pull(sample)
sharedGenes <- intersect(rownames(getData(cObjMul, "counts")), rownames(MGAS.Counts))
multiplets <- cbind(getData(cObjMul, "counts")[sharedGenes, ], MGAS.Counts[sharedGenes, sections])
singlets <- getData(cObjSng, "counts")[sharedGenes, ]
features <- which(rownames(multiplets) %in% rownames(getData(cObjMul, "counts"))[getData(cObjMul, "features")])

cObjMul	 <- CIMseqMultiplets(
  multiplets,
  cbind(getData(cObjMul, "counts.ercc"), MGAS.CountsERCC[, sections]),
  features
)
cObjSng <- CIMseqSinglets(
  singlets,
  getData(cObjSng, "counts.ercc"),
  getData(cObjSng, "dim.red"),
  getData(cObjSng, "classification")
)
#cObjMul <- CIMseqMultiplets(
#  getData(cObjMul, "counts")[, 1:2],
#  getData(cObjMul, "counts.ercc")[, 1:2],
#  getData(cObjMul, "features")
#)

#RUN
#future::plan(multiprocess)
options(future.wait.interval = 10000.0)
options(future.wait.timeout = 1e9)
future::plan(
  future.batchtools::batchtools_slurm,
  template = "/crex/proj/snic2018-8-151/private/batchtools.slurm.tmpl",
  resources = list(
    account = "snic2018-8-151", partition = "core", ntasks = 1L,
    time = "28:00:00", jobname = "CIMseq",
    modules = "R_packages/3.5.1", R = "R/3.5.1", log.file = file.path(currPath, "logs/slurm.txt")
  ),
  workers = 100
)

#run deconvolution
print(paste0("Starting deconvolution at ", Sys.time()))
sObj <- CIMseqSwarm(
  cObjSng, cObjMul, maxiter = 100, swarmsize = 500, nSyntheticMultiplets = 400
)
print(paste0("Finished deconvolution at ", Sys.time()))

save(sObj, file = file.path(currPath, "data/sObj.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spSwarm.SI.txt"))
print("finished")
