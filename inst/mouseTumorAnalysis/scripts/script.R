packages <- c("sp.scRNAseq", "sp.scRNAseqData", "seqTools", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#only tumor cells
st <- str_detect(colnames(countsMgfpTumor), "^s")
cObjSng <- spCounts(countsMgfpTumor[, st], countsMgfpTumorERCC[, st])

print(paste0("Starting only tumor analysis at ", Sys.time()))
uObj_tumor <- spUnsupervised(cObjSng, max_iter = 7000, initial_dims = sum(st), seed = 324)
print(paste0("Done only tumor analysis at ", Sys.time()))

#all cells
sn <- str_detect(colnames(countsMgfp), "^s")
commonGenes <- intersect(rownames(countsMgfp), rownames(countsRegev))
commonGenes <- intersect(commonGenes, rownames(countsMgfpTumor))

sng <- cbind(countsMgfp[commonGenes, sn], countsMgfpTumor[commonGenes, st], countsRegev[commonGenes, ])
mul <- cbind(countsMgfp[commonGenes, !sn], countsMgfpTumor[commonGenes, !st])

erccSng <- cbind(
  countsMgfpERCC[, sn],
  countsMgfpTumorERCC[, st],
  matrix(NA, nrow = nrow(countsMgfpERCC), ncol = ncol(countsRegev))
)
erccMul <- cbind(countsMgfpERCC[, !sn], countsMgfpTumorERCC[, !st])

cObjSng <- spCounts(sng, erccSng)

print(paste0("Starting all cells analysis at ", Sys.time()))
uObj_all <- spUnsupervised(cObjSng, max_iter = 7000, initial_dims = ncol(sng), seed = 324)
print(paste0("Done all cells analysis at ", Sys.time()))

#save
save(uObj_tumor, uObj_all, file = file.path(currPath, "data/uObjs.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))
