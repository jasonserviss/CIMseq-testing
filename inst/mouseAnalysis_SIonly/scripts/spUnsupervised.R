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

boolMulSi <- colnames(mul) %in% filter(countsMgfpMeta, tissue == "SI")$sample
boolSngSi <- colnames(sng) %in% filter(countsMgfpMeta, tissue == "SI")$sample

#setup spCounts
cObjSngSi <- spCounts(sng[, boolSngSi], erccSng[, boolSngSi])
cObjMulSi <- spCounts(mul[, boolMulSi], erccMul[, boolMulSi])

#only small intestine
print(paste0("Starting SI cells analysis at ", Sys.time()))
uObjSi <- spUnsupervised(cObjSngSi, max_iter = 7000, initial_dims = sum(boolSngSi), max = 1000, seed = 348592)
print(paste0("Done SI cells analysis at ", Sys.time()))

#rename classes
plotUnsupervisedMarkers(
  uObjSi, cObjSngSi,
  c("Lgr5", "Muc2", "Ptprc", "Chga", "Alpi", "Lyz1", "Dclk1"),
  pal = RColorBrewer::brewer.pal(7, "Set2")
)

classes <- tibble(
  oldClass = getData(uObjSi, "classification"),
  newClass = case_when(
    oldClass %in% c("P1", "K1", "N1", "B1", "L1") ~ "SI.Stem",
    oldClass %in% c("Q1", "O1") ~ "SI.Enterocyte",
    oldClass %in% c("M1") ~ "SI.Goblet",
    oldClass %in% c("R1") ~ "SI.Paneth",
    oldClass %in% c("H1") ~ "SI.Tufft",
    oldClass %in% c("M1") ~ "Endocrine",
    oldClass %in% c("D1") ~ "Blood",
    TRUE ~ "error"
  )
)

classification(uObjSi) <- classes$newClass
tsneMeans(uObjSi) <- tsneGroupMeans(getData(uObjSi, "tsne"), getData(uObjSi, "classification"))

#save
#save(uObjSi, file = file.path(currPath, "data/uObjSi.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))
