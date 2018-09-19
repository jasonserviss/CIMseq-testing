packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

#check package version
algoV <- sessionInfo()$otherPkgs$sp.scRNAseq$Version
last3 <- paste(
  strsplit(sessionInfo()$otherPkgs$sp.scRNAseq$Version, "\\.")[[1]][2:4],
  collapse = ""
)
if(!as.numeric(last3) >= 102) {
  stop("sp.scRNAseq package version too low. Must be >= 0.1.0.2")
}

currPath <- getwd()

#setup spCounts
keep.plates.SI <- c("NJA01202", "NJA01301", "NJA01302")
keep.plates.colon <- c("NJA01303", "NJA01401", "NJA01205", "NJA00609")

s <- str_detect(colnames(countsMgfp), "^s")
e <- colnames(countsMgfp) %in% filter(countsMgfpMeta, is.na(GFP) & !filtered & plate %in% c(keep.plates.SI, keep.plates.colon))$sample
boolSng <- s & e
boolMul <- !s & e

#setup spCounts
cObjSng <- spCounts(countsMgfp[, boolSng], countsMgfpERCC[, boolSng])

#all cells
print(paste0("Starting all cells analysis at ", Sys.time()))
uObj <- spUnsupervised(
  cObjSng, max_iter = 3000, initial_dims = sum(boolSng), seed = 7976898,
  max = 1500, perplexity = 30, pcVarPercent = 0.35, kNN = 50, distCut = 0.75,
  classCut = 4
)
print(paste0("Done all cells analysis at ", Sys.time()))

#save
save(uObj, file = file.path(currPath, "data/uObj.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))
