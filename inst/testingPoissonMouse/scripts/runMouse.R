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
cObjSng <- spCounts(sng, erccSng)
cObjMul <- spCounts(mul[, 1:2], erccMul[, 1:2])

print("spCounts done")

#spUnsupervised
if(file.exists(file.path(currPath, 'data/uObj.rda'))) {
  load(file.path(currPath, "data/uObj.rda"))
} else {
  uObj <- spUnsupervised(cObjSng, max_iter = 7000, initial_dims = ncol(sng), seed = 87689)
  save(uObj, file = file.path(currPath, "data/uObj.rda"))
}

print("spUnsupervised done")

##spSwarm
future::plan(multiprocess)

print(paste0("Starting deconvolution at ", Sys.time()))
sObj <- spSwarm(cObjSng, cObjMul, uObj, maxiter = 100, swarmsize = 500, nSyntheticMultiplets = 400)
print(paste0("Finished deconvolution at ", Sys.time()))

save(sObj, file = file.path(currPath, "data/sObj.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))
print("finished")
