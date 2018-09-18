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
cObjSngC <- spCounts(sng[, boolSngC], erccSng[, boolSngC])
cObjMulC <- spCounts(mul[, boolMulC], erccMul[, boolMulC])

#only colon
print(paste0("Starting colon cells analysis at ", Sys.time()))
uObjC <- spUnsupervised(cObjSngC, max_iter = 7000, initial_dims = sum(boolSngC), max = 1000, seed = 32932)
print(paste0("Done colon cells analysis at ", Sys.time()))

#rename classes
#plotUnsupervisedMarkers(uObjC, cObjSngC, c("Lgr5", "Muc2", "Ptprc", "Chga", "Slc40a1", "Dclk1"), pal = RColorBrewer::brewer.pal(8, "Set1"))
#
#classes <- tibble(
#  oldClass = getData(uObjC, "classification"),
#  newClass = case_when(
#    oldClass %in% c("I1", "A1", "J1", "C1", "D1", "H1") ~ "C.Stem",
#    oldClass %in% c("E1", "F1") ~ "C.Colonocyte",
#    oldClass %in% c("K1", "B1", "G1") ~ "C.Goblet",
#    oldClass %in% c("N1") ~ "Endocrine",
#    oldClass %in% c("M1") ~ "Blood",
#    oldClass %in% c("L1") ~ "C.Tufft",
#    TRUE ~ "error"
#  )
#)
#
#classification(uObjC) <- classes$newClass
#tsneMeans(uObjC) <- tsneGroupMeans(getData(uObjC, "tsne"), getData(uObjC, "classification"))

#save
save(uObjC, file = file.path(currPath, "data/uObjC.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))
