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

#cObjSngSi <- spCounts(sng[, !boolSngSi], erccSng[, !boolSngSi])
#cObjSngC <- spCounts(sng[, boolSngC], erccSng[, boolSngC])

#cObjMulSi <- spCounts(mul[, !boolMulSi], erccMul[, !boolMulSi])
#cObjMulC <- spCounts(mul[, boolMulC], erccMul[, boolMulC])

#only colon
#print(paste0("Starting colon cells analysis at ", Sys.time()))
#uObjC <- spUnsupervised(cObjSngC, max_iter = 7000, initial_dims = sum(boolSngC), max = 1000, seed = 876901)
#print(paste0("Done colon cells analysis at ", Sys.time()))

#only small intestine
#print(paste0("Starting SI cells analysis at ", Sys.time()))
#uObjSi <- spUnsupervised(cObjSngSi, max_iter = 7000, initial_dims = sum(!boolSngC), max = 1000, seed = 23)
#print(paste0("Done SI cells analysis at ", Sys.time()))

#all cells
print(paste0("Starting all cells analysis at ", Sys.time()))
uObj <- spUnsupervised(cObjSng, max_iter = 5000, initial_dims = ncol(sng), seed = 7234834)
print(paste0("Done all cells analysis at ", Sys.time()))

#rename classes
#classes <- tibble(
#  oldClass = getData(uObj, "classification"),
#  newClass = case_when(
#    oldClass %in% c("Q1", "O1", "N1", "R1", "B1", "I1", "T1", "L1") ~ "SI.Stem",
#    oldClass %in% c("F1", "H1") ~ "SI.Enterocyte",
#    oldClass %in% c("G1", "P1") ~ "SI.Goblet",
#    oldClass %in% c("S1") ~ "SI.Paneth",
#    oldClass %in% c("K1") ~ "SI.Tufft",
#    oldClass %in% c("A1", "E1") ~ "C.Stem",
#    oldClass %in% c("J1") ~ "C.Colonocyte",
#    oldClass %in% c("C1") ~ "C.Goblet",
#    oldClass %in% c("M1") ~ "Endocrine",
#    oldClass %in% c("D1") ~ "Blood",
#    TRUE ~ "error"
#  )
#)

#classification(uObj) <- classes$newClass
#tsneMeans(uObj) <- tsneGroupMeans(getData(uObj, "tsne"), getData(uObj, "classification"))

#save
#save(uObjC, file = file.path(currPath, "data/uObjC.rda"))
#save(uObjSi, file = file.path(currPath, "data/uObjSi.rda"))
save(uObj, file = file.path(currPath, "data/uObj.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))
