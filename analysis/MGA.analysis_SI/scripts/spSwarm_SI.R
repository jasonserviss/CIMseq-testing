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

boolMulC <- colnames(mul) %in% filter(countsMgfpMeta, tissue == "colon")$sample
boolSngC <- colnames(sng) %in% filter(countsMgfpMeta, tissue == "colon")$sample

#setup spCounts
cObjSngSi <- spCounts(sng[, !boolSngC], erccSng[, !boolSngC])
cObjMulSi <- spCounts(mul[, !boolMulC], erccMul[, !boolMulC])

print("spCounts done")

#spUnsupervised
if(file.exists(file.path(currPath, 'data/uObjSi.rda'))) {
  load(file.path(currPath, 'data/uObjSi.rda'))
}

print("spUnsupervised done")

##spSwarm
future::plan(multiprocess)
selectIdx <- spTopVar(cObjSngSi, 2000)

print(paste0("Starting deconvolution at ", Sys.time()))
sObjSi <- spSwarm(
  cObjSngSi, cObjMulSi, uObjSi, maxiter = 100, swarmsize = 500,
  nSyntheticMultiplets = 400, selectInd = selectIdx
)
print(paste0("Finished deconvolution at ", Sys.time()))

save(sObjSi, file = file.path(currPath, "data/sObjSi.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.spSwarm_Si.txt"))
print("finished")
