
#PACKAGES
packages <- c(
  "CIMseq", "sp.scRNAseqData", "tidyverse", "future", "future.apply", "Seurat"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#check package version
algoV <- sessionInfo()$otherPkgs$CIMseq$Version
last3 <- paste(strsplit(algoV, "\\.")[[1]][2:4], collapse = "")
if(!as.numeric(last3) >= 100) {
  stop("sp.scRNAseq package version too low. Must be >= 0.2.0.0")
}

currPath <- getwd()

##spCounts
s <- str_detect(colnames(SCM.Counts), "^s")
sng <- SCM.Counts[, s]
sngERCC <- SCM.CountsERCC[, s]
mul <- SCM.Counts[, !s]
mulERCC <- SCM.CountsERCC[, !s]
meta <- filter(SCM.Meta, sample %in% colnames(sng))

#Dimensionality reduction and classification
print(paste0("Starting dim.red and classification analysis at ", Sys.time()))
mca <- CreateSeuratObject(
  raw.data = sng, meta.data = meta, project = "Enge_only"
)
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
  object = mca, pc.genes = mca@var.genes, do.print = FALSE
)

#PCElbowPlot(object = mca, num.pc = 20) + scale_x_continuous(breaks = seq(0, 20, 1))

set.seed(7239832)

mca <- FindClusters(
  object = mca, reduction.type = "pca", dims.use = 1:6, resolution = 1,
  save.SNN = TRUE, n.start = 100, nn.eps = 0.5, print.output = FALSE,
  force.recalc = TRUE
)
mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = 1:6, min_dist = 0.5,
  n_neighbors = 15
)

#DimPlot(
#  object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
#  vector.friendly = FALSE, pt.size = 1
#) + scale_colour_manual(values = col40())

#FeaturePlot(
#  mca, c("CD74", "ANXA3", "ACTG2"),
#  reduction.use = "umap", dark.theme = FALSE, pt.size = 1,
#  vector.friendly = FALSE
#)
print(paste0("Done all dim.red and classification analysis at ", Sys.time()))

#extract data
if(identical(colnames(sng), rownames(FetchData(mca, "ident")))) {
  classes <- as.character(FetchData(mca, "ident")[[1]])
} else {
  stop("Naming mismatch")
}
var.genes <- mca@var.genes
select <- which(rownames(sng) %in% var.genes)
dim.red <- mca@dr$umap@cell.embeddings
colnames(dim.red) <- NULL

#setup CIMseqData objects
cObjSng <- CIMseqSinglets(sng, sngERCC, dim.red, classes)
cObjMul <- CIMseqMultiplets(mul, mulERCC, select)

save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))

##spSwarm
future::plan(multiprocess)

print(paste0("Starting deconvolution at ", Sys.time()))
sObj <- CIMseqSwarm(
  cObjSng, cObjMul, maxiter = 10, swarmsize = 150, nSyntheticMultiplets = 400
)
print(paste0("Finished deconvolution at ", Sys.time()))

save(sObj, file = file.path(currPath, "data/sObj.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))
print("finished")
