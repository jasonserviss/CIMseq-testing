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
print("spCounts done")

#run spUnsupervised
uObj <- spUnsupervised(cObjSng, max_iter = 5000, initial_dims = ncol(sng), seed = 976723)

#plotUnsupervisedMarkers(
#  uObj, cObjSng,
#  c("Lgr5", "Muc2", "Ptprc", "Chga", "Alpi", "Lyz1", "Dclk1", "Slc40a1"),
#  pal = RColorBrewer::brewer.pal(8, "Set1")
#)

classes <- tibble(
  sample = rownames(getData(uObj, "tsne")),
  oldClass = getData(uObj, "classification")
) %>%
mutate(newClass = case_when(
  oldClass %in% c("P1", "O1", "U1", "M1", "Q1", "J1", "A1") ~ "SI.Stem",
  oldClass %in% c("T1", "F1") ~ "SI.Enterocyte",
  oldClass %in% c("H1") ~ "SI.Goblet",
  oldClass %in% c("L1") ~ "SI.Tufft",
  oldClass %in% c("V1") ~ "SI.Paneth",
  oldClass %in% c("B1", "D1", "S1") ~ "C.Stem",
  oldClass %in% c("E1") ~ "C.Colonocyte",
  oldClass %in% c("K1", "G1", "I1") ~ "C.Goblet",
  oldClass %in% c("C1") ~ "Blood",
  oldClass %in% c("R1", "N1") ~ "Endocrine",
  TRUE ~ "error"
))

classification(uObj) <- classes$newClass[match(rownames(getData(uObj, "tsne")), classes$sample)]
groupMeans(uObj) <- averageGroupExpression(cObjSng, getData(uObj, "classification"), FALSE)
tsneMeans(uObj) <- tsneGroupMeans(getData(uObj, "tsne"), getData(uObj, "classification"))

save(uObj, file = file.path(currPath, "data/uObj.rda"))
print("spUnsupervised done")
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))
print("finished")
