packages <- c("CIMseq", "CIMseq.data", "tidyverse", "Seurat", "harmony", "future.apply")
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
keep.plates.SI <- c(
  "NJA01201", "NJA01202", "NJA01301", "NJA01302", "NJA01501"
)
s <- str_detect(colnames(MGA.Counts), "^s")
samples <- filter(MGA.Meta, !filtered & unique_key %in% keep.plates.SI)$sample
e <- colnames(MGA.Counts) %in% samples
boolSng <- s & e
boolMul <- !s & e
boolRSI <- colnames(RSI.Counts) %in% filter(RSI.Meta, !filtered)$sample

iGenes <- intersect(intersect(rownames(RSI.Counts), rownames(MGA.Counts)), rownames(TMD.Counts))

singlets <- cbind(MGA.Counts[iGenes, boolSng], RSI.Counts[iGenes, boolRSI])
singletERCC <- cbind(MGA.CountsERCC[, boolSng], RSI.CountsERCC[, boolRSI])
multiplets <- MGA.Counts[iGenes, boolMul]
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
  do.plot = FALSE, x.low.cutoff = 1, y.cutoff = 1
)

mca <- ScaleData(
  object = mca, display.progress = FALSE, do.par = TRUE,
  num.cores = 4, vars.to.regress = "nUMI"
)
norm <- as.matrix(mca@data)
top2000 <- CIMseq::selectTopMax(norm, 2000)
VlnPlot(object = mca, features = c("nGene", "nUMI"), group.by = "source")

pcor <- as.matrix(CIMseq::pearsonsDist(norm[select, ]))
pca.results <- irlba::irlba(A = pcor, nv = 100)
gene.loadings <- pca.results$v
sdev <- pca.results$d/sqrt(max(1, ncol(pcor) - 1))
mca <- RunPCA(object = mca)

mca@dr$pca@cell.embeddings <- pca$x
mca@dr$pca@gene.loadings <- pca$rotation
mca@dr$pca@sdev <- pca$sdev

DimPlot(
  object = mca, reduction.use = "pca", dim.1 = 1, dim.2 = 2, 
  no.legend = FALSE, do.return = TRUE, group.by = "source",
  vector.friendly = FALSE, pt.size = 1
)
#PCElbowPlot(object = mca, num.pc = 100) + scale_x_continuous(breaks = seq(0, 100, 5))
pcs <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11)

######################################################################################
#check for PC 
pc.cor <- cor(mca@dr$pca@cell.embeddings, as.numeric(as.factor(mca@meta.data$source)))[pcs ,]
#PC9 seems to be the culprit
pcs <- pcs[which(pc.cor < 0.25)]
######################################################################################

print(paste0("Using ", max(pcs), " principal components."))

mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.3,
  n_neighbors = 15, seed.use = 9823493
)

#mca <- RunTSNE(mca, dims.use = pcs, gene.use = mca@var.genes, seed.use = 78239)

mca <- FindClusters(
  object = mca, reduction.type = "pca", dims.use = pcs, resolution = 0.65,
  n.start = 100, n.iter = 1000, nn.eps = 0, k.param = 30, prune.SNN = 1/15,
  algorithm = 1, save.SNN = TRUE, print.output = FALSE, plot.SNN = FALSE,
  force.recalc = TRUE, random.seed = 93820
)

DimPlot(
  object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
  vector.friendly = FALSE, pt.size = 1
) + scale_colour_manual(values = col40())

mca@meta.data %>%
  group_by(source, res.0.6) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(source) %>%
  mutate(`%` = n / sum(n) * 100) %>%
  ggplot() +
  geom_bar(aes(res.0.6, `%`, fill = source), stat = "identity", position = position_dodge(width = 1)) +
  facet_wrap(~res.0.6, scales = "free") +
  labs(x = "Class", y = "% of dataset")

mca@meta.data %>%
  count(source, res.0.6) %>%
  ggplot() +
  geom_bar(aes(res.0.6, n, fill = source), stat = "identity", position = position_dodge(width = 1)) +
  facet_wrap(~res.0.6, scales = "free") +
  labs(x = "Class", y = "Count")

mca@dr$umap@cell.embeddings %>%
  matrix_to_tibble("sample") %>%
  mutate(source = case_when(
    str_detect(sample, "SRR654") ~ "Tabula Muris",
    str_detect(sample, "SRR510") ~ "Regev",
    TRUE ~ "Enge"
  )) %>%
  sample_n(nrow(.), FALSE) %>%
  ggplot() +
  geom_point(aes(UMAP1, UMAP2, colour = source), alpha = 0.75)

FeaturePlot(
  mca,
  c("Lgr5", "Ptprc", "Chga", "Dclk1", "Atoh1", "Lyz1", "Alpi", "Mki67"),
  reduction.use = "umap", dark.theme = FALSE, pt.size = 0.1,
  vector.friendly = FALSE
)

FeaturePlot(
  mca,
  c("Lgr5", "Alpi", "Mki67"),
  reduction.use = "umap", dark.theme = FALSE, pt.size = 0.1,
  vector.friendly = FALSE
)
# 
# FeaturePlot(
#   mca,
#   c("Plet1"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.1,
#   vector.friendly = FALSE
# )

test <- FindMarkers(
  mca, ident.1 = 1, ident.2 = NULL, 
  only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = log(1.5),
  test.use = "wilcox"
)

FeaturePlot(
  mca,
  c("Kpna2", "Tuba1c", "Msmo1"),
  reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
  vector.friendly = FALSE
)

markers <- FindAllMarkers(
  object = mca, only.pos = TRUE, min.diff.pct = 0.25, logfc.threshold = log(1.5),
  test.use = "roc"
)

DoHeatmap(
  object = mca, genes.use = unique(markers$gene), slim.col.label = TRUE, remove.key = TRUE,
  group.label.rot = TRUE, cex.row = 1
)

print(paste0("Done all cells analysis at ", Sys.time()))

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
print(paste0("saving data to ", currPath, "."))
save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_CIMseqData.txt"))