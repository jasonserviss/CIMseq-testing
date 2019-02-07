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
s <- str_detect(colnames(MGA.Counts), "^s")
si.samples <- filter(
  MGA.Meta,
  !filtered &
  !str_detect(unique_key, "NJA004") &
  sub_tissue == "small_intestine" &
  subject_strain == "C57BL/6J"
)$sample
keep.plates.colon <- c(
  "NJA01203", "NJA01205","NJA01303", "NJA01401", "NJA01503", "NJA01504",
  "NJA01801", "NJA01803"
)
c.samples <- filter(
  MGA.Meta, !filtered & unique_key %in% keep.plates.colon
)$sample

e <- colnames(MGA.Counts) %in% c(si.samples, c.samples)
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
pcs <- 1:26

var <- mca@dr$pca@sdev^2
var.percent <- var / sum(var) * 100
barplot(
  var.percent, xlab = "PC", ylab = "Percent Variance", 
  names.arg = 1:length(var.percent), las = 1, ylim = c(0, max(var.percent) + 2),
  col = "gray"
)

print(paste0("Using ", max(pcs), " principal components."))

mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.5,
  n_neighbors = 10, seed.use = 984324
)

mca <- RunTSNE(mca, dims.use = pcs, seed.use = 2387, perplexity = 30)

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
  c("Lgr5", "Ptprc", "Chga", "Dclk1", "Alpi", "Slc26a3", "Atoh1", "Lyz1", "Mki67", "Hoxb13"),
  reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
  vector.friendly = FALSE
)
FeaturePlot(
  mca,
  c("Lgr5", "Slc26a3", "Alpi", "Mki67", "Hoxb13"),
  reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
  vector.friendly = FALSE
)

matrix_to_tibble(mca@dr$umap@cell.embeddings, "sample") %>%
  inner_join(MGA.Meta, by = "sample") %>%
  ggplot() +
  geom_point(aes(UMAP1, UMAP2, colour = unique_key)) +
  scale_colour_manual(values = c(col40(), "black"))
