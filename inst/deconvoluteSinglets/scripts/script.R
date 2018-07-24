#PACKAGES
packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

s <- str_detect(colnames(countsMgfp), "^s")
commonGenes <- intersect(rownames(countsMgfp), rownames(countsRegev))

sng <- cbind(countsMgfp[commonGenes, s], countsRegev[commonGenes, ])
mul <- countsMgfp[commonGenes, !s]

erccSng <- cbind(
  countsMgfpERCC[, s], 
  matrix(NA, nrow = nrow(countsMgfpERCC), ncol = ncol(countsRegev))
)
erccMul <- cbind(countsMgfpERCC[, !s])

#setup spCounts
cObjSng <- spCounts(sng, erccSng)
#cObjMul <- spCounts(mul, erccMul)

print("spCounts done")

#spUnsupervised
if(file.exists(file.path(currPath, 'data/uObj.rda'))) {
  load(file.path(currPath, "data/uObj.rda"))
}

print("spUnsupervised done")

##spSwarm
future::plan(multiprocess)
selectIdx <- spTopVar(cObjSng, 2000)

print(paste0("Starting deconvolution at ", Sys.time()))
sObj <- spSwarm(
  cObjSng, cObjSng, uObj, maxiter = 100, swarmsize = 500,
  nSyntheticMultiplets = 400, selectInd = selectIdx
)
print(paste0("Finished deconvolution at ", Sys.time()))

save(sObj, file = file.path(currPath, "data/sObj.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spSwarm.txt"))
print("finished")
