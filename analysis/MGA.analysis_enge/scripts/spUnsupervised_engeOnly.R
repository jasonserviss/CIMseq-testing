packages <- c("CIMseq", "sp.scRNAseqData", "tidyverse", "Seurat")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

#check package version
algoV <- sessionInfo()$otherPkgs$CIMseq$Version
last3 <- paste(strsplit(algoV, "\\.")[[1]][2:4], collapse = "")
if(!as.numeric(last3) >= 100) {
  stop("sp.scRNAseq package version too low. Must be >= 0.2.0.0")
}

currPath <- getwd()

#load individual analysis from each tissue
load('../MGA.analysis_colon/data/CIMseqData.rda')
cObjSng.colon <- cObjSng
cObjMul.colon <- cObjMul
rm(cObjSng); rm(cObjMul)
load('../MGA.analysis_SI/data/CIMseqData.rda')
cObjSng.SI <- cObjSng
cObjMul.SI <- cObjMul
rm(cObjSng); rm(cObjMul)

#combine
singlets <- cbind(
  getData(cObjSng.SI, "counts"),
  getData(cObjSng.colon, "counts")
)
singletERCC <- cbind(
  getData(cObjSng.SI, "counts.ercc"),
  getData(cObjSng.colon, "counts.ercc")
)
multiplets <- cbind(
  getData(cObjMul.SI, "counts"),
  getData(cObjMul.colon, "counts")
)
multipletERCC <- cbind(
  getData(cObjMul.SI, "counts.ercc"),
  getData(cObjMul.colon, "counts.ercc")
)
classes <- c(
  paste(getData(cObjSng.SI, "classification"), "SI", sep = "."),
  paste(getData(cObjSng.colon, "classification"), "C", sep = ".")
)
names(classes) <- colnames(singlets)

#perform dimensionality reduction with all data
mca <- CreateSeuratObject(raw.data = singlets)
mca <- NormalizeData(
  object = mca, normalization.method = "LogNormalize", scale.factor = 1e6
)
mca <- FindVariableGenes(
  object = mca, mean.function = ExpMean, dispersion.function = LogVMR,
  do.plot = FALSE, x.low.cutoff = 1
)
mca <- ScaleData(
  object = mca, genes.use = mca@var.genes, display.progress = FALSE, do.par = TRUE,
  num.cores = 4
)
mca <- RunPCA(
  object = mca, pc.genes = mca@var.genes, pcs.compute = 100, do.print = FALSE,
  seed.use = 983209
)
mca <- JackStraw(object = mca, num.replicate = 100, display.progress = TRUE, num.pc = 50)
mca <- JackStrawPlot(object = mca, PCs = 1:50)
PCp <- mca@dr$pca@jackstraw@overall.p.values
pcs <- PCp[PCp[, 2] < 10^-9, 1]
print(paste0("Using ", max(pcs), " principal components."))

mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.5,
  n_neighbors = 10, seed.use = 79356
)

# DimPlot(
#   object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
#   vector.friendly = FALSE, pt.size = 1
# ) + scale_colour_manual(values = col40())
#
#
# FeaturePlot(
#   mca,
#   c("Lgr5", "Ptprc", "Chga", "Dclk1", "Alpi", "Slc26a3", "Atoh1", "Lyz1", "Mki67", "Hoxb13"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
#   vector.friendly = FALSE
# )

markers <- FindAllMarkers(
  object = mca, only.pos = TRUE, min.diff.pct = 0.2, logfc.threshold = log(2),
  test.use = "roc"
)
print(paste0("Done all cells analysis at ", Sys.time()))

#extract data
var.genes <- unique(markers$gene)
select <- which(rownames(singlets) %in% var.genes)
dim.red <- mca@dr$umap@cell.embeddings
colnames(dim.red) <- NULL

#setup CIMseq objects
cObjSng <- CIMseqSinglets(singlets, singletERCC, dim.red, classes)
cObjMul <- CIMseqMultiplets(multiplets, multipletERCC, select)

#save
save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))
