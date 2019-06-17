packages <- c("CIMseq", "CIMseq.data", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
path <- '../MGA.analysis_enge/data/CIMseqData.rda'
if(file.exists(path)) {
  load(path)
}
rm(path)

colon.samples <- filter(MGA.Meta, sub_tissue == "colon")$sample
boolSng <- colnames(getData(cObjSng, "counts")) %in% colon.samples
boolMul <- colnames(getData(cObjMul, "counts"))  %in% colon.samples
cObjSng <- CIMseqSinglets(
  getData(cObjSng, "counts")[, boolSng],
  getData(cObjSng, "counts.ercc")[, boolSng],
  getData(cObjSng, "dim.red")[boolSng, ],
  getData(cObjSng, "classification")[boolSng]
)

cObjMul <- CIMseqMultiplets(
  getData(cObjMul, "counts")[, boolMul],
  getData(cObjMul, "counts.ercc")[, boolMul],
  getData(cObjMul, "features")
)

save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_CIMseqData.txt"))