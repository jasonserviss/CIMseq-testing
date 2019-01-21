packages <- c("CIMseq", "CIMseq.data", "Seurat", "stringr", "dplyr")
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
#keep.plates.SI <- c(
#  "NJA01201", "NJA01202", "NJA01301", "NJA01302", "NJA01501, NJA01504"
#)
s <- str_detect(colnames(MGA.Counts), "^s")
samples <- filter(MGA.Meta,
  !filtered &
  !str_detect(unique_key, "NJA004") &
  sub_tissue == "small_intestine" &
  subject_strain == "C57BL/6J"
)$sample
e <- colnames(MGA.Counts) %in% samples
boolSng <- s & e
boolMul <- !s & e
singlets <- MGA.Counts[, boolSng]
singletERCC <- MGA.CountsERCC[, boolSng]
multiplets <- MGA.Counts[, boolMul]
multipletERCC <- MGA.CountsERCC[, boolMul]

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
#mca <- JackStraw(object = mca, num.replicate = 100, display.progress = TRUE, num.pc = 50)
#mca <- JackStrawPlot(object = mca, PCs = 1:50)
#PCp <- mca@dr$pca@jackstraw@overall.p.values
#pcs <- PCp[PCp[, 2] < 10^-9, 1]
pcs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20)
print(paste0("Using ", max(pcs), " principal components."))

mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.5,
  n_neighbors = 10, seed.use = 79356
)

mca <- FindClusters(
  object = mca, reduction.type = "pca", dims.use = pcs, resolution = 0.6,
  n.start = 100, n.iter = 1000, nn.eps = 0, k.param = 20, prune.SNN = 1/15,
  algorithm = 1, save.SNN = TRUE, print.output = FALSE, plot.SNN = FALSE,
  force.recalc = TRUE, random.seed = 93820
)

# VlnPlot(object = mca, features.plot = c("nGene", "nUMI"), nCol = 2)
#
DimPlot(
  object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
  vector.friendly = FALSE, pt.size = 1
) + scale_colour_manual(values = col40())
#
#
# FeaturePlot(
#   mca,
#   c("Lgr5", "Ptprc", "Chga", "Dclk1", "Alpi", "Atoh1", "Lyz1", "Mki67"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
#   vector.friendly = FALSE
# )

#find differentially expressed genes
markers <- FindAllMarkers(
  object = mca, only.pos = TRUE, min.diff.pct = 0.2, logfc.threshold = log(2),
  test.use = "roc"
)
# table(markers$cluster)
# DoHeatmap(
#   object = mca, genes.use = unique(markers$gene), slim.col.label = TRUE,
#   remove.key = TRUE, group.label.rot = TRUE
# )

# < 5 markers for class 1. Merge with nearest cluster.

#calculate centroids
# s.data <- RNAseqFunctions::cpm(singlets)
# ident <- as.character(FetchData(mca, "ident")$ident)
# centroids <- sapply(unique(ident), function(i) {
#   matrixStats::rowSums2(s.data[, ident == i])
# })
# #calculate distances between all centroids
# dists <- as.matrix(dist(t(centroids), diag = TRUE, upper = TRUE))
# class1.dists <- dists[, "1"]
# class1.dists <- class1.dists[names(class1.dists) != "1"]
# mergeInto <- names(class1.dists)[class1.dists == min(class1.dists)]
#
# #manually merge clusters 1 and 2. They have few DE genes
# old <- mca@ident
# n <- names(old)
# old <- as.character(old)
# new <- case_when(
#   old %in% c("1", mergeInto) ~ "1",
#   old == "0" ~ "0",
#   TRUE ~ as.character(as.numeric(old) - 1)
# )
# names(new) <- n
# mca@ident <- as.factor(new)

print(paste0("Done all cells analysis at ", Sys.time()))

#extract data
singlets <- singlets[, colnames(singlets) %in% colnames(mca@data)]
singletERCC <- singletERCC[, colnames(singletERCC) %in% colnames(singlets)]
idx <- match(rownames(FetchData(mca, "ident")), colnames(singlets))
classes <- as.character(FetchData(mca, "ident")[[1]])[idx]
names(classes) <- rownames(FetchData(mca, "ident"))[idx]
var.genes <- unique(markers$gene)
select <- which(rownames(singlets) %in% var.genes)
dim.red <- mca@dr$umap@cell.embeddings
colnames(dim.red) <- NULL

#setup CIMseqData objects
cObjSng <- CIMseqSinglets(singlets, singletERCC, dim.red, classes)
cObjMul <- CIMseqMultiplets(multiplets, multipletERCC, select)

#save
if(!"data" %in% list.dirs(currPath, full.names = FALSE)) system('mkdir data')
save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))
