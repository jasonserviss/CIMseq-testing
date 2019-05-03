packages <- c("CIMseq", "CIMseq.data", "tidyverse", "Seurat", "harmony", "future.apply")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

devtools::install_github("immunogenomics/harmony", ref = "6f07162", dependencies = FALSE)
library(harmony)

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
boolRSI <- colnames(RSI.Counts) %in% filter(RSI.Meta, !filtered)$sample
boolTMD <- colnames(TMD.Counts) %in% filter(TMD.Meta, !filtered)$sample

iGenes <- intersect(intersect(rownames(RSI.Counts), rownames(MGA.Counts)), rownames(TMD.Counts))


TMsub <- function(data, reference){
  cs <- matrixStats::colSums2(reference)
  #if(!all(cs > matrixStats::colSums2(data)))
  expanded <- apply(data, 2, function(r) rep(rownames(data), r))
  sapply(seq_along(expanded), function(i) {
    out <- rep(0, nrow(data))
    names(out) <- rownames(data)
    n <- sample(cs, 1, FALSE)
    s <- sample(expanded[[i]], n, TRUE)
    c <- table(s)
    out[match(names(c), names(out))] <- c
    out
  })
}

# .check <- function(data, reference) {
#   tibble(
#     type = c(rep("data", ncol(data)), rep("ref", ncol(reference))),
#     cs = matrixStats::colSums2(cbind(data, reference))
#   ) %>%
#     ggplot() +
#     geom_histogram(aes(cs, fill = type))
# }

tmSub <- TMsub(TMD.Counts[iGenes, boolTMD], MGA.Counts[iGenes, boolSng])
colnames(tmSub) <- paste0(colnames(TMD.Counts)[boolTMD], ".sub")
singlets <- cbind(MGA.Counts[iGenes, boolSng], RSI.Counts[iGenes, boolRSI], tmSub[iGenes, ])
singletERCC <- cbind(MGA.CountsERCC[, boolSng], RSI.CountsERCC[, boolRSI], TMD.CountsERCC[, boolTMD])
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
pcs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
         18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
         34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46)

print(paste0("Using ", max(pcs), " principal components."))
#PCElbowPlot(object = mca, num.pc = 100) + scale_x_continuous(breaks = seq(0, 100, 5))

mca <- RunHarmony(
  object = mca, group.by.vars = "source", theta = 3, plot_convergence = TRUE,
  nclust = 50, max.iter.cluster = 50, max.iter.harmony = 50
)
# DimPlot(object = mca, reduction.use = "harmony", pt.size = 0.1, group.by = "source", do.return = TRUE)

mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.5,
  n_neighbors = 20, seed.use = 79354
)

#mca <- RunTSNE(mca, dims.use = pcs, gene.use = mca@var.genes, seed.use = 78239)

mca <- FindClusters(
  object = mca, reduction.type = "pca", dims.use = pcs, resolution = 2.3,
  n.start = 100, n.iter = 1000, nn.eps = 0, k.param = 30, prune.SNN = 1/15,
  algorithm = 1, save.SNN = TRUE, print.output = FALSE, plot.SNN = FALSE,
  force.recalc = TRUE, random.seed = 93820
)

# DimPlot(
#   object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
#   vector.friendly = FALSE, pt.size = 1
# ) + scale_colour_manual(values = col40())
# 
# mca@meta.data %>%
#   group_by(source, res.2.3) %>%
#   summarize(n = n()) %>%
# ungroup() %>%
# group_by(source) %>%
# mutate(`%` = n / sum(n) * 100) %>%
# ggplot() +
# geom_bar(aes(res.2.3, `%`, fill = source), stat = "identity", position = position_dodge(width = 1)) +
# facet_wrap(~res.2.3, scales = "free") +
# labs(x = "class", y = "% of dataset")
# 
# mca@dr$umap@cell.embeddings %>%
#   matrix_to_tibble("sample") %>%
#   mutate(source = case_when(
#     str_detect(sample, "SRR654") ~ "Tabula Muris",
#     str_detect(sample, "SRR510") ~ "Regev",
#     TRUE ~ "Enge"
#   )) %>%
#   sample_n(nrow(.), FALSE) %>%
#   ggplot() +
#   geom_point(aes(UMAP1, UMAP2, colour = source), alpha = 0.75)
# 
# FeaturePlot(
#   mca,
#   c("Lgr5", "Ptprc", "Chga", "Dclk1", "Slc26a3", "Atoh1", "Lyz1", "Alpi", "Hoxb13", "Mki67"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.1,
#   vector.friendly = FALSE
# )
# 
# FeaturePlot(
#   mca,
#   c("Lgr5", "Slc26a3", "Hoxb13", "Mki67"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.1,
#   vector.friendly = FALSE
# )
# 
# FeaturePlot(
#   mca,
#   c("Plet1"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.1,
#   vector.friendly = FALSE
# )

markers <- FindAllMarkers(
  object = mca, only.pos = TRUE, min.diff.pct = 0.4, logfc.threshold = log(2),
  test.use = "roc"
)

# DoHeatmap(
#   object = mca, genes.use = unique(markers$gene), slim.col.label = TRUE, remove.key = TRUE,
#   group.label.rot = TRUE, cex.row = 1
# )

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