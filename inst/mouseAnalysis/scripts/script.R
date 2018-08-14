packages <- c("sp.scRNAseq", "sp.scRNAseqData", "seqTools", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#setup spCounts
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
cObjMul <- spCounts(mul, erccMul)

cObjSngSi <- spCounts(sng[, !boolSngC], erccSng[, !boolSngC])
cObjSngC <- spCounts(sng[, boolSngC], erccSng[, boolSngC])

cObjMulSi <- spCounts(mul[, !boolMulC], erccMul[, !boolMulC])
cObjMulC <- spCounts(mul[, boolMulC], erccMul[, boolMulC])

#only colon
print(paste0("Starting colon cells analysis at ", Sys.time()))
uObjC <- spUnsupervised(cObjSngC, max_iter = 7000, initial_dims = sum(boolSngC), max = 1000, seed = 876901)
print(paste0("Done colon cells analysis at ", Sys.time()))

#only small intestine
print(paste0("Starting SI cells analysis at ", Sys.time()))
uObjSi <- spUnsupervised(cObjSngSi, max_iter = 7000, initial_dims = sum(!boolSngC), max = 1000, seed = 23)
print(paste0("Done SI cells analysis at ", Sys.time()))

#all cells
print(paste0("Starting all cells analysis at ", Sys.time()))
uObj <- spUnsupervised(cObjSng, max_iter = 7000, initial_dims = ncol(sng), seed = 876901)
print(paste0("Done all cells analysis at ", Sys.time()))

#save
save(uObjC, uObjSi, uObj, file = file.path(currPath, "data/uObjs.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))
