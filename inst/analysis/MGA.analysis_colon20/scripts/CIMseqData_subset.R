packages <- c("CIMseq", "CIMseq.data", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
path1 <- '../MGA.analysis_enge20/data/seurat.rda'
path2 <- '../MGA.analysis_enge20/data/CIMseqData.rda'

if(file.exists(path1)) {
  load(path1)
}
if(file.exists(path2)) {
  load(path2)
}
rm(list=c("path1", "path2"))

#subset small intestine samples
keep.plates.colon <- c(
  "NJA01203", "NJA01205","NJA01303", "NJA01401", "NJA01503", "NJA01504",
  "NJA01801", "NJA01803", "NJD00101", "NJD00102", "NJD00103", "NJD00104",
  "NJA01901", "NJA02001"
)
sng.samples <- filter(MGA.Meta, !filtered & unique_key %in% keep.plates.colon & cellNumber == "Singlet")$sample
stopifnot(all(sng.samples %in% colnames(mca@raw.data)))
mca.c <- SubsetData(mca, sng.samples)

#detect features in small intestine subset
markers <- FindAllMarkers(
  object = mca.c, only.pos = TRUE, min.diff.pct = 0.2, logfc.threshold = log(2),
  test.use = "roc"
)
select <- which(rownames(getData(cObjSng, "counts")) %in% unique(markers$gene))

#setup CIMseqData objects
cObjSng <- CIMseqSinglets(
  getData(cObjSng, "counts")[, sng.samples],
  getData(cObjSng, "counts.ercc")[, sng.samples],
  getData(cObjSng, "dim.red")[sng.samples, ],
  getData(cObjSng, "classification")[sng.samples]
)

mul.samples <- filter(MGA.Meta, !filtered & unique_key %in% keep.plates.colon & cellNumber == "Multiplet")$sample
cObjMul <- CIMseqMultiplets(
  getData(cObjMul, "counts")[, mul.samples],
  getData(cObjMul, "counts.ercc")[, mul.samples],
  select
)

#save
save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))
save(mca.c, markers, file = file.path(currPath, "data/seurat.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_CIMseqData.txt"))