packages <- c("CIMseq", "sp.scRNAseqData", "Seurat", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

#check package version
algoV <- sessionInfo()$otherPkgs$CIMseq$Version
last3 <- paste(strsplit(algoV, "\\.")[[1]][2:4], collapse = "")
if(!as.numeric(last3) >= 100) {
  stop("sp.scRNAseq package version too low. Must be >= 0.2.0.0")
}

currPath <- getwd()

#setup spCounts
keep.plates.SI <- c("NJA01202", "NJA01301", "NJA01302", "NJA01501")
keep.plates.colon <- c(
  "NJA01303", "NJA01401", "NJA01205", "NJA01203", "NJA00609", "NJA01503"
)
keep <- c(keep.plates.SI, keep.plates.colon)

s <- str_detect(colnames(MGA.Counts), "^s")
samples <- filter(MGA.Meta, !filtered & unique_key %in% keep)$sample
e <- colnames(MGA.Counts) %in% samples
boolSng <- s & e
boolMul <- !s & e
singlets <- MGA.Counts[, boolSng]
singletERCC <- MGA.CountsERCC[, boolSng]
multiplets <- MGA.Counts[, boolMul]
multipletERCC <- MGA.CountsERCC[, boolMul]
meta <- filter(MGA.Meta, sample %in% colnames(singlets))

#Dimensionality reduction and classification
print(paste0("Starting all cells analysis at ", Sys.time()))
mca <- CreateSeuratObject(
  raw.data = singlets, meta.data = meta, project = "Enge_only"
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
  object = mca, pc.genes = mca@var.genes, pcs.compute = 100, do.print = FALSE
)
mca <- JackStraw(object = mca, num.replicate = 100, display.progress = TRUE, num.pc = 50)
mca <- JackStrawPlot(object = mca, PCs = 1:20)
PCp <- mca@dr$pca@jackstraw@overall.p.values
pcs <- PCp[PCp[, 2] < 0.001, 1]
print(paste0("Using ", max(pcs), " principal components."))
#PCElbowPlot(object = mca, num.pc = 100) + scale_x_continuous(breaks = seq(0, 100, 5))

mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.5,
  n_neighbors = 15, seed.use = 79356
)
mca <- FindClusters(
  object = mca, reduction.type = "pca", dims.use = pcs, resolution = 2,
  save.SNN = TRUE, n.start = 100, nn.eps = 0.5, print.output = FALSE,
  force.recalc = TRUE, random.seed = 93820
)

# DimPlot(
#   object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
#   vector.friendly = FALSE, pt.size = 1
# ) + scale_colour_manual(values = col40())
#
# FeaturePlot(
#   mca,
#   c("Lgr5", "Ptprc", "Chga", "Dclk1", "Slc26a3", "Atoh1", "Lyz1", "Alpi", "Hoxb13", "Mki67"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.1,
#   vector.friendly = FALSE
# )

markers <- FindAllMarkers(
  object = mca, only.pos = TRUE, min.diff.pct = 0.4, logfc.threshold = log(2),
  test.use = "roc"
)
# n <- min(table(markers$cluster))
# toPlot <- markers %>%
#   group_by(cluster) %>%
#   top_n(max(c(n, 100)), myAUC) %>%
#   pull(gene)
#
# DoHeatmap(
#   object = mca, genes.use = toPlot, slim.col.label = TRUE, remove.key = TRUE,
#   group.label.rot = TRUE
# )

print(paste0("Done all cells analysis at ", Sys.time()))

#extract data
if(identical(colnames(singlets), rownames(FetchData(mca, "ident")))) {
  classes <- as.character(FetchData(mca, "ident")[[1]])
} else {
  stop("Naming mismatch")
}
select <- match(markers$gene, rownames(singlets))
names(select) <- markers$cluster
dim.red <- mca@dr$umap@cell.embeddings
colnames(dim.red) <- NULL

#setup CIMseqData objects
cObjSng <- CIMseqSinglets(singlets, singletERCC, dim.red, classes)

####!!!!!!
#Maybe I should use nTopVar for feature selection instead of Seurat for deconvolution
cObjMul <- CIMseqMultiplets(multiplets, multipletERCC, select)

#save
save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))
