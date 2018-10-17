packages <- c("sp.scRNAseq", "sp.scRNAseqData", "seqTools", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#setup spCounts
keep.plates.SI <- c("NJA01202", "NJA01301", "NJA01302")
keep.plates.colon <- c("NJA01303", "NJA01401", "NJA01205", "NJA00609")
keep <- c(keep.plates.SI, keep.plates.colon)

s <- str_detect(colnames(MGA.Counts), "^s")
samples <- filter(MGA.Meta, !filtered & plate %in% keep & sub_tissue == "small_intestine")$sample
e <- colnames(MGA.Counts) %in% samples
boolSng <- s & e
boolMul <- !s & e

#setup spCounts
cObjSng <- spCounts(MGA.Counts[, boolSng], MGA.CountsERCC[, boolSng])

#only small intestine
print(paste0("Starting SI cells analysis at ", Sys.time()))
uObjSi <- spUnsupervised(
  cObjSng, max_iter = 3000, initial_dims = sum(boolSng), seed = 7976898,
  max = 1500, perplexity = 30, pcVarPercent = 0.35, kNN = 50, distCut = 0.75,
  classCut = 4
)
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
