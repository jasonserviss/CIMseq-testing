packages <- c("sp.scRNAseq", "sp.scRNAseqData", "seqTools", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#setup spCounts
s <- str_detect(colnames(countsMgfp), "^s")
e <- colnames(countsMgfp) %in% filter(countsMgfpMeta, is.na(GFP) & !filtered)$sample
boolSng <- s & e
boolMul <- !s & e

#setup spCounts
cObjSng <- spCounts(countsMgfp[, boolSng], countsMgfpERCC[, boolSng])
cObjMul <- spCounts(countsMgfp[, boolMul], countsMgfpERCC[, boolMul])

#all cells
print(paste0("Starting all cells analysis at ", Sys.time()))
uObj <- spUnsupervised(cObjSng, max_iter = 5000, initial_dims = ncol(sng), seed = 7234834, max = 1000)
print(paste0("Done all cells analysis at ", Sys.time()))

#rename classes
#paneth <- plotUnsupervisedMarkers(uObj, cObjSng, "Lyz1") %>%
#  plotData() %>%
#  filter(Classification == "L1" & `t-SNE dim 1` < -100) %>%
#  pull(Sample)
#
#classes <- tibble(
#  sample = rownames(getData(uObj, "tsne")),
#  oldClass = getData(uObj, "classification"),
#  newClass = case_when(
#    oldClass %in% c("O1", "N1", "J1", "K1") ~ "SI.Stem",
#    oldClass %in% c("M1") ~ "SI.Enterocyte",
#    oldClass %in% c("L1") & sample %in% paneth ~ "SI.Paneth",
#    oldClass %in% c("L1") ~ "SI.Goblet",
#    oldClass %in% c("E1") ~ "SI.Tufft",
#    oldClass %in% c("H1", "A1", "I1", "C1") ~ "C.Stem",
#    oldClass %in% c("D1") ~ "C.Colonocyte",
#    oldClass %in% c("B1") ~ "C.Goblet",
#    oldClass %in% c("G1") ~ "Endocrine",
#    oldClass %in% c("F1") ~ "Blood",
#    TRUE ~ "error"
#  )
#)
#
#classification(uObj) <- classes$newClass
#tsneMeans(uObj) <- tsneGroupMeans(getData(uObj, "tsne"), getData(uObj, "classification"))

#save
#save(uObjC, file = file.path(currPath, "data/uObjC.rda"))
#save(uObjSi, file = file.path(currPath, "data/uObjSi.rda"))
save(uObj, file = file.path(currPath, "data/uObj_engeOnly.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised_engeOnly.txt"))
