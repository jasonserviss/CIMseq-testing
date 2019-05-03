packages <- c("CIMseq", "CIMseq.data", "Seurat", "stringr", "dplyr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#setup CIMseqSinglets
singlets <- TMD.Counts

#Dimensionality reduction and classification
print(paste0("Starting all cells analysis at ", Sys.time()))
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
# mca <- JackStraw(object = mca, num.replicate = 100, display.progress = TRUE, num.pc = 50, do.par = TRUE, num.cores = 8)
# mca <- JackStrawPlot(object = mca, PCs = 1:50)
# PCp <- mca@dr$pca@jackstraw@overall.p.values
# pcs <- PCp[PCp[, 2] < 10^-6, 1]
pcs <- c(
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
  18, 19, 20, 21, 22, 23, 24, 25
)
print(paste0("Using ", max(pcs), " principal components."))

mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.3,
  n_neighbors = 20, seed.use = 56424
)

mca <- FindClusters(
  object = mca, reduction.type = "pca", dims.use = pcs, resolution = 1,
  n.start = 100, n.iter = 1000, nn.eps = 0, k.param = 20, prune.SNN = 1/15,
  algorithm = 1, save.SNN = TRUE, print.output = FALSE, plot.SNN = FALSE,
  force.recalc = TRUE, random.seed = 93820
)

DimPlot(
  object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
  vector.friendly = FALSE, pt.size = 1
) + scale_colour_manual(values = col40())

FeaturePlot(
  mca,
  c("Lgr5", "Ptprc", "Chga", "Dclk1", "Slc26a3", "Atoh1", "Lyz1", "Alpi", "Hoxb13", "Mki67"),
  reduction.use = "umap", dark.theme = FALSE, pt.size = 0.1,
  vector.friendly = FALSE
)

#find differentially expressed genes
# markers <- FindAllMarkers(
#   object = mca, only.pos = TRUE, min.diff.pct = 0.2, logfc.threshold = log(2),
#   test.use = "roc"
# )

save(mca, file = "data/seuratObj.rda", compress = "bzip2")
