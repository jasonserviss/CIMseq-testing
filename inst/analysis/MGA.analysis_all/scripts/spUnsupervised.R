packages <- c("CIMseq", "CIMseq.data", "seqTools", "tidyverse", "Seurat", "harmony")
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
s <- str_detect(colnames(MGA.Counts), "^s")
keep.plates <- c(
  "NJA01203", "NJA01205","NJA01303", "NJA01401", "NJA01503", "NJA01504",
  "NJA01801", "NJA01803", "NJD00101", "NJD00102", "NJD00103", "NJD00104",
  "NJA01201", "NJA01202", "NJA01301", "NJA01302", "NJA01501", "NJA01901",
  "NJA02001"
)
samples <- filter(MGA.Meta, !filtered & unique_key %in% keep.plates)$sample

e <- colnames(MGA.Counts) %in% samples
boolSng <- s & e
boolMul <- !s & e

iGenes <- intersect(intersect(rownames(RSI.Counts), rownames(MGA.Counts)), rownames(TMD.Counts))

singlets <- cbind(MGA.Counts[iGenes, boolSng], RSI.Counts[iGenes, ], TMD.Counts[iGenes, ])
singletERCC <- cbind(MGA.CountsERCC[, boolSng], RSI.CountsERCC, TMD.CountsERCC)
multiplets <- MGA.Counts[, boolMul]
multipletERCC <- MGA.CountsERCC[, boolMul]

#Dimensionality reduction and classification
print(paste0("Starting all cells analysis at ", Sys.time()))
mca <- CreateSeuratObject(raw.data = singlets)
mca@meta.data$source <- case_when(
  str_detect(rownames(mca@meta.data), "SRR654") | str_detect(rownames(mca@meta.data), "SRR510") ~ "External",
  str_detect(rownames(mca@meta.data), "NJA") | str_detect(rownames(mca@meta.data), "NJD") ~ "Enge",
  TRUE ~ "error"
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
# mca <- JackStraw(object = mca, num.replicate = 100, display.progress = TRUE, num.pc = 50)
# mca <- JackStrawPlot(object = mca, PCs = 1:50)
# PCp <- mca@dr$pca@jackstraw@overall.p.values
# pcs <- PCp[PCp[, 2] < 10^-6, 1]
pcs <- c(
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
  18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
  34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 47
)
print(paste0("Using ", max(pcs), " principal components."))
#PCElbowPlot(object = mca, num.pc = 100) + scale_x_continuous(breaks = seq(0, 100, 5))

mca <- RunHarmony(
  object = mca, group.by.vars = "source", theta = 3, plot_convergence = TRUE, 
  nclust = 50, max.iter.cluster = 50, max.iter.harmony = 50
)
DimPlot(object = mca, reduction.use = "harmony", pt.size = 0.1, group.by = "source", do.return = TRUE)

mca <- RunUMAP(
  object = mca, reduction.use = "harmony", dims.use = pcs, min_dist = 0.5,
  n_neighbors = 20, seed.use = 79356
)

#mca <- RunTSNE(mca, dims.use = pcs, gene.use = mca@var.genes, seed.use = 78239)

mca <- FindClusters(
  object = mca, reduction.type = "harmony", dims.use = pcs, resolution = 2,
  n.start = 100, n.iter = 1000, nn.eps = 0, k.param = 30, prune.SNN = 1/15,
  algorithm = 1, save.SNN = TRUE, print.output = FALSE, plot.SNN = FALSE,
  force.recalc = TRUE, random.seed = 93820
)

DimPlot(
  object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
  vector.friendly = FALSE, pt.size = 1
) + scale_colour_manual(values = col40())

mca@dr$umap@cell.embeddings %>%
  matrix_to_tibble("sample") %>%
  inner_join(rownames_to_column(mca@meta.data, "sample")) %>%
  ggplot() +
  geom_point(aes(UMAP1, UMAP2, colour = source), alpha = 0.5)

FeaturePlot(
  mca,
  c("Lgr5", "Ptprc", "Chga", "Dclk1", "Slc26a3", "Atoh1", "Lyz1", "Alpi", "Hoxb13", "Mki67"),
  reduction.use = "umap", dark.theme = FALSE, pt.size = 0.1,
  vector.friendly = FALSE
)

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


#save
#save(uObjC, file = file.path(currPath, "data/uObjC.rda"))
#save(uObjSi, file = file.path(currPath, "data/uObjSi.rda"))
save(uObj, file = file.path(currPath, "data/uObj.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))
