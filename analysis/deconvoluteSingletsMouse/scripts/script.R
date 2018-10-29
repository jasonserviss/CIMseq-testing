#PACKAGES
packages <- c("CIMseq", "sp.scRNAseqData", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

s <- str_detect(colnames(MGA.Counts), "^s")
commonGenes <- intersect(rownames(MGA.Counts), rownames(RSI.Counts))

sng <- cbind(MGA.Counts[commonGenes, s], RSI.Counts[commonGenes, ])
mul <- MGA.Counts[commonGenes, !s]

erccSng <- cbind(MGA.CountsERCC[, s], RSI.CountsERCC)
erccMul <- cbind(MGA.CountsERCC[, !s])

#setup spCounts
cObjSng <- CIMseqSinglets(sng, erccSng)
set.seed(29843)
cObjMul <- spCounts(sng[, sample(1:ncol(sng), 100, replace = FALSE)], erccMul)

print("spCounts done")

#spUnsupervised
if(file.exists(file.path('~/Github/sp.scRNAseqTesting/inst/testingPoissonMouse', 'data/uObj.rda'))) {
  load(file.path('~/Github/sp.scRNAseqTesting/inst/testingPoissonMouse', "data/uObj.rda"))
} else {
  stop("Could not load spUnsupervised data")
}

print("spUnsupervised done")

##spSwarm
selectIdx <- spTopVar(cObjSng, 2000)
print(paste0("Starting deconvolution at ", Sys.time()))

future::plan(multiprocess)
sObj <- spSwarm(
  cObjSng, cObjMul, uObj, maxiter = 100, swarmsize = 500,
  nSyntheticMultiplets = 400, selectInd = selectIdx
)
print(paste0("Finished deconvolution at ", Sys.time()))

save(sObj, file = file.path(currPath, "data/sObj.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spSwarm.txt"))
print("finished")
