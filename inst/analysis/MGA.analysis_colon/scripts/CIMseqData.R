packages <- c("CIMseq", "CIMseq.data", "tidyverse", "Seurat")
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
keep.plates.colon <- c(
  "NJA01203", "NJA01205","NJA01303", "NJA01401", "NJA01503", "NJA01504",
  "NJA01801", "NJA01803", "NJD00101", "NJD00102", "NJD00103", "NJD00104",
  "NJA01901", "NJA02001"
)

s <- str_detect(colnames(MGA.Counts), "^s")
samples <- filter(MGA.Meta, !filtered & unique_key %in% keep.plates.colon)$sample
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
mca@meta.data <- column_to_rownames(merge(mca@meta.data, column_to_rownames(MGA.Meta, "sample"), by = 0), "Row.names")
mca <- NormalizeData(
  object = mca, normalization.method = "LogNormalize", scale.factor = 1e6
)
mca <- FindVariableGenes(
  object = mca, mean.function = ExpMean, dispersion.function = LogVMR,
  do.plot = FALSE, x.low.cutoff = 0.5, x.high.cutoff = Inf
)
load('../TMD.analysis/data/DE.rda')
genes <- de %>% 
  mutate(mean = map_dbl(gene, function(g) mean(singlets[g, ]))) %>%
  filter(mean > 1) %>%
  mutate(pct.diff = abs(pct.1 - pct.2)) %>%
  group_by(cluster) %>% 
  top_n(100, pct.diff) %>% 
  pull(gene) %>% 
  unique() 
#"Tgm3"    "Klk1"    "Slc15a1" "Saa2"    "Atp12a" 
mca@var.genes <- genes

mca <- ScaleData(
  object = mca, genes.use = mca@var.genes, display.progress = FALSE, do.par = TRUE,
  num.cores = 4
)
mca <- RunPCA(
  object = mca, pc.genes = mca@var.genes, pcs.compute = 50, do.print = FALSE,
  seed.use = 983209
)

PCAPlot(object = mca, dim.1 = 1, dim.2 = 2) +
  scale_colour_manual(values = col40())
mca <- JackStraw(object = mca, num.replicate = 100, display.progress = TRUE, num.pc = 50, do.par = TRUE, num.cores = 8)
mca <- JackStrawPlot(object = mca, PCs = 1:50)
PCp <- mca@dr$pca@jackstraw@overall.p.values
pcs <- PCp[PCp[, 2] < 10^-6, 1]
pcs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
        18, 19, 20, 21, 22, 24, 25, 26)
# VizPCA(object = mca, pcs.use = pcs[1:10])
var <- mca@dr$pca@sdev^2
var.percent <- var / sum(var) * 100
barplot(
  var.percent, xlab = "PC", ylab = "Percent Variance",
  names.arg = 1:length(var.percent), las = 1, ylim = c(0, max(var.percent) + 2),
  col = "gray"
)

pcGeneCorr <- function(PCs, gene, cpm) {
  c <- cor(PCs, cpm[gene, ])
  plot(c)
  abline(h = 0, col = "red", lty = 2)
  imax <- which(c[, 1] == max(c[, 1]))
  print(paste0(rownames(c)[imax], " highest correlation at ", c[imax, ]))
}
pcGeneCorr(mca@dr$pca@cell.embeddings, "Hoxb13", singlets)
print(paste0("Using ", max(pcs), " principal components."))

#umap representation at: mca@dr$umap@cell.embeddings
# library(uwot)
# 
# dim.code <- GetDimReduction(object = mca, reduction.type = "pca", slot = "key")
# dim.codes <- paste0(dim.code, pcs[1:10])
# data.use <- GetDimReduction(object = mca, reduction.type = "pca", slot = "cell.embeddings")
# data.use <- data.use[, dim.codes, drop = FALSE]
# #23823
# set.seed(87856)
# u <- umap(
#   data.use, n_neighbors = 15, metric = "euclidean", n_epochs = 500, 
#   spread = 1.65, min_dist = 0.1
# )
# rownames(u) <- rownames(data.use)
# colnames(u) <- c("UMAP1", "UMAP2")
# 
# #add classes and plot
# matrix_to_tibble(u, "sample") %>%
#   full_join(tibble(sample = names(mca@ident), class = as.character(mca@ident))) %>%
#   ggplot() +
#   geom_point(aes(UMAP1, UMAP2, colour = class)) +
#   scale_colour_manual(values = col40())
# 
# genes <- c("Lgr5", "Ptprc", "Chga", "Dclk1", "Slc26a3", "Atoh1", "Mki67", "Hoxb13")
# genes <- c("Lgr5", "Slc26a3", "Mki67", "Hoxb13")
# 
# matrix_to_tibble(t(CIMseq:::.norm.log.counts(singlets)[genes, ]), "sample") %>%
#   full_join(matrix_to_tibble(u, "sample")) %>%
#   gather(gene, value, -sample, -UMAP1, -UMAP2) %>%
#   ggplot() +
#   geom_point(aes(UMAP1, UMAP2, colour = value), size = 1) +
#   facet_wrap(~gene) +
#   scale_colour_viridis_c()
# 
# matrix_to_tibble(u, "sample") %>%
#   inner_join(MGA.Meta) %>%
#   ggplot() +
#   geom_point(aes(V1, V2, colour = unique_key)) +
#   scale_colour_manual(values = c(col40(), "black")) +
#   theme(legend.position = "top")

mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs[1:10], min_dist = 0.5,
  n_neighbors = 15, seed.use = 7937395
)
#mca@dr$umap@cell.embeddings <- u
# 
#mca <- RunTSNE(mca, dims.use = pcs, seed.use = 2387, perplexity = 30)

mca <- FindClusters(
  object = mca, reduction.type = "pca", dims.use = pcs[1:10], resolution = 1,
  n.start = 100, n.iter = 1000, nn.eps = 0, k.param = 20, prune.SNN = 1/15,
  algorithm = 1, save.SNN = TRUE, print.output = FALSE, plot.SNN = FALSE,
  force.recalc = TRUE, random.seed = 93820
)

# heatmap(as.matrix(mca@snn), scale = "none")
#VlnPlot(object = mca, features.plot = c("nGene", "nUMI"), nCol = 2)

DimPlot(
  object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
  vector.friendly = FALSE, pt.size = 1
) + scale_colour_manual(values = col40())

FeaturePlot(
  mca,
  c("Lgr5", "Ptprc", "Chga", "Dclk1", "Slc26a3", "Atoh1", "Mki67", "Hoxb13"),
  reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
  vector.friendly = FALSE
)

FeaturePlot(
  mca,
  c("Lgr5", "Slc26a3", "Mki67", "Hoxb13"),
  reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
  vector.friendly = FALSE
)

matrix_to_tibble(mca@dr$umap@cell.embeddings, "sample") %>%
  inner_join(MGA.Meta, by = "sample") %>%
  ggplot() +
  geom_point(aes(UMAP1, UMAP2, colour = unique_key)) +
  scale_colour_manual(values = c(col40(), "black"))

#find differentially expressed genes
markers <- FindAllMarkers(
  object = mca, only.pos = TRUE, min.diff.pct = 0.2, logfc.threshold = log(2),
  test.use = "roc"
)

DoHeatmap(
  object = mca, genes.use = unique(markers$gene), slim.col.label = TRUE,
  remove.key = TRUE, group.label.rot = TRUE
)

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
save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))
