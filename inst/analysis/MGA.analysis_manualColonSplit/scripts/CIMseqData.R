packages <- c("CIMseq", "CIMseq.data", "Seurat", "stringr", "dplyr", "readr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

#check package version
algoV <- sessionInfo()$otherPkgs$CIMseq$Version
last3 <- paste(strsplit(algoV, "\\.")[[1]][2:4], collapse = "")
if(!as.numeric(last3) >= 100) {
  stop("sp.scRNAseq package version too low. Must be >= 0.2.0.0")
}

currPath <- getwd()

#setup CIMseqData
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
singlets <- MGA.Counts[, boolSng]
singletERCC <- MGA.CountsERCC[, boolSng]
multiplets <- MGA.Counts[, boolMul]
multipletERCC <- MGA.CountsERCC[, boolMul]

colon.sub <- singlets[, colnames(singlets) %in% filter(MGA.Meta, sub_tissue == "colon")$sample]

#####split prox and dist pre-classification
loc <- apply(colon.sub, 2, function(c) if_else(c['Hoxb13'] > c['Osr2'], "distal", "proximal"))

##########################################################################
#proximal
##########################################################################

col.p <- CreateSeuratObject(raw.data = colon.sub[, loc == "proximal"])
col.p <- NormalizeData(
  object = col.p, normalization.method = "LogNormalize", scale.factor = 1e6
)
x.low.cutoff <- 1
y.cutoff <- 1
col.p <- FindVariableGenes(
  object = col.p, mean.function = ExpMean, dispersion.function = LogVMR,
  do.plot = FALSE, x.low.cutoff = x.low.cutoff, y.cutoff = y.cutoff
)

# col.p@hvg.info %>% 
#   mutate(selected = if_else(gene.mean > x.low.cutoff & gene.dispersion.scaled > y.cutoff, TRUE, FALSE)) %>% 
#   ggplot() + 
#   geom_point(aes(gene.mean, gene.dispersion.scaled, colour = selected)) +
#   geom_hline(yintercept = y.low.cutoff, lty = 2)

col.p <- ScaleData(
  object = col.p, genes.use = col.p@var.genes, display.progress = FALSE, do.par = TRUE,
  num.cores = 4
)
col.p <- RunPCA(
  object = col.p, pc.genes = col.p@var.genes, pcs.compute = 100, do.print = FALSE,
  seed.use = 983209
)

# PCAPlot(object = col.p, dim.1 = 1, dim.2 = 2) +
#   scale_colour_manual(values = col40())
# col.p <- JackStraw(object = col.p, num.replicate = 100, display.progress = TRUE, num.pc = 50)
# col.p <- JackStrawPlot(object = col.p, PCs = 1:50)
# PCp <- col.p@dr$pca@jackstraw@overall.p.values
# pcs <- PCp[PCp[, 2] < 10^-6, 1]
pcs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)

col.p <- RunUMAP(
  object = col.p, reduction.use = "pca", dims.use = pcs, min_dist = 0.5,
  n_neighbors = 10, seed.use = 984335
)

col.p <- FindClusters(
  object = col.p, reduction.type = "pca", dims.use = pcs[1:5], resolution = 0.35,
  n.start = 100, n.iter = 1000, nn.eps = 0, k.param = 20, prune.SNN = 1/15,
  algorithm = 1, save.SNN = TRUE, print.output = FALSE, plot.SNN = FALSE,
  force.recalc = TRUE, random.seed = 93820
)

# DimPlot(
#   object = col.p, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
#   vector.friendly = FALSE, pt.size = 1
# ) + scale_colour_manual(values = col40())
# 
# FeaturePlot(
#   col.p,
#   c("Lgr5", "Ptprc", "Chga", "Dclk1", "Slc26a3", "Atoh1", "Mki67", "Hoxb13"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
#   vector.friendly = FALSE
# )
# 
# FeaturePlot(
#   col.p,
#   c("Lgr5", "Slc26a3", "Mki67", "Hoxb13", "Atoh1"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
#   vector.friendly = FALSE
# )
# FeaturePlot(
#   col.p,
#   c("Lgr5", "Slc26a3", "Mboat1"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
#   vector.friendly = FALSE
# )

col.p.ident <- as.numeric(col.p@ident)
names(col.p.ident) <- names(col.p@ident)

###########################################################################
#distal
##########################################################################

col.d <- CreateSeuratObject(raw.data = colon.sub[, loc == "distal"])
col.d <- NormalizeData(
  object = col.d, normalization.method = "LogNormalize", scale.factor = 1e6
)
x.low.cutoff <- 1
y.cutoff <- 1
col.d <- FindVariableGenes(
  object = col.d, mean.function = ExpMean, dispersion.function = LogVMR,
  do.plot = FALSE, x.low.cutoff = x.low.cutoff, y.cutoff = y.cutoff
)

# col.d@hvg.info %>% 
#   mutate(selected = if_else(gene.mean > x.low.cutoff & gene.dispersion.scaled > y.cutoff, TRUE, FALSE)) %>% 
#   ggplot() + 
#   geom_point(aes(gene.mean, gene.dispersion.scaled, colour = selected)) +
#   geom_hline(yintercept = y.low.cutoff, lty = 2)

col.d <- ScaleData(
  object = col.d, genes.use = col.d@var.genes, display.progress = FALSE, do.par = TRUE,
  num.cores = 4
)
col.d <- RunPCA(
  object = col.d, pc.genes = col.d@var.genes, pcs.compute = 100, do.print = FALSE,
  seed.use = 983209
)

# PCAPlot(object = col.d, dim.1 = 1, dim.2 = 2) +
#   scale_colour_manual(values = col40())
# col.d <- JackStraw(object = col.d, num.replicate = 100, display.progress = TRUE, num.pc = 50)
# col.d <- JackStrawPlot(object = col.d, PCs = 1:50)
# PCp <- col.d@dr$pca@jackstraw@overall.p.values
# pcs <- PCp[PCp[, 2] < 10^-6, 1]
c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

col.d <- RunUMAP(
  object = col.d, reduction.use = "pca", dims.use = pcs, min_dist = 0.5,
  n_neighbors = 10, seed.use = 984335
)

col.d <- FindClusters(
  object = col.d, reduction.type = "pca", dims.use = pcs, resolution = 0.55,
  n.start = 100, n.iter = 1000, nn.eps = 0, k.param = 20, prune.SNN = 1/15,
  algorithm = 1, save.SNN = TRUE, print.output = FALSE, plot.SNN = FALSE,
  force.recalc = TRUE, random.seed = 93820
)

# DimPlot(
#   object = col.d, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
#   vector.friendly = FALSE, pt.size = 1
# ) + scale_colour_manual(values = col40())
# 
# FeaturePlot(
#   col.d,
#   c("Lgr5", "Slc26a3", "Mki67", "Hoxb13", "Atoh1"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
#   vector.friendly = FALSE
# )

col.d.ident <- as.numeric(col.d@ident)
names(col.d.ident) <- names(col.d@ident)

##########################################################################
#everything
##########################################################################

mca <- CreateSeuratObject(raw.data = singlets)
mca <- NormalizeData(
  object = mca, normalization.method = "LogNormalize", scale.factor = 1e6
)
mca <- FindVariableGenes(
  object = mca, mean.function = ExpMean, dispersion.function = LogVMR,
  do.plot = TRUE, x.low.cutoff = 0.5, y.cutoff = 0.5
)
mca <- ScaleData(
  object = mca, genes.use = mca@var.genes, display.progress = FALSE, do.par = TRUE,
  num.cores = 4
)
mca <- RunPCA(
  object = mca, pc.genes = mca@var.genes, pcs.compute = 100, do.print = FALSE,
  seed.use = 983209
)
# mca <- JackStraw(object = mca, num.replicate = 100, display.progress = TRUE, num.pc = 50)
# mca <- JackStrawPlot(object = mca, PCs = 1:50)
# PCp <- mca@dr$pca@jackstraw@overall.p.values
# pcs <- PCp[PCp[, 2] < 10^-6, 1]
pcs <- c(
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
  18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 31, 34
)

mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.45,
  n_neighbors = 10, seed.use = 984335
)

mca <- FindClusters(
  object = mca, reduction.type = "pca", dims.use = pcs, resolution = 1.8,
  n.start = 100, n.iter = 1000, nn.eps = 0, k.param = 20, prune.SNN = 1/15,
  algorithm = 1, save.SNN = TRUE, print.output = FALSE, plot.SNN = FALSE,
  force.recalc = TRUE, random.seed = 93820
)

# DimPlot(
#   object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
#   vector.friendly = FALSE, pt.size = 1
# ) + scale_colour_manual(values = col40())

old <- as.numeric(mca@ident)
names(old) <- names(mca@ident)

###########################################################################
#transfer
##########################################################################

old[match(names(col.p.ident), names(old))] <- col.p.ident + max(old)
old[match(names(col.d.ident), names(old))] <- col.d.ident + max(old)
new.ident <- parse_factor(as.character(old), levels = unique(old))
names(new.ident) <- names(old)
new <- mca
new@ident <- new.ident

# DimPlot(
#   object = new, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
#   vector.friendly = FALSE, pt.size = 1
# ) + scale_colour_manual(values = col40())
#
# singlets[c("Lgr5", "Ptprc", "Chga", "Dclk1", "Slc26a3", "Atoh1", "Mki67", "Hoxb13"), ] %>%
#   matrix_to_tibble("gene") %>%
#   gather(sample, value, -gene) %>%
#   inner_join(rownames_to_column(FetchData(object = new, vars.all = "ident"), "sample")) %>%
#   ggplot() +
#   geom_violin(
#     data = . %>% filter(value != 0),
#     aes(ident, value, fill = ident, colour = ident)
#   ) +
#   geom_jitter(
#     data = . %>% filter(value != 0),
#     aes(ident, value), size = 0.1, height = 0, width = 0.5
#   ) +
#   geom_text(
#     data = . %>% group_by(gene, ident) %>% summarize(p = paste0(round(100 * (length(which(value == 0)) / length(value)), digits = 2), "%")),
#     aes(ident, 0, label = p), size = 2
#   ) +
#   facet_wrap(~gene, scales = "free_y") +
#   scale_fill_manual(values = col40()) +
#   scale_colour_manual(values = col40()) +
#   guides(fill = FALSE, colour = FALSE) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))


# FeaturePlot(
#   new,
#   c("Hoxb13", "Osr2"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.5,
#   vector.friendly = FALSE
# )

#markers
markers <- FindAllMarkers(
  object = new, only.pos = TRUE, min.diff.pct = 0.2, logfc.threshold = log(2),
  test.use = "roc"
)
# markers %>% group_by(gene) %>% filter(n() == 1) %>% ungroup() %>% count(cluster) %>% print(n = nrow(.))
# DoHeatmap(
#   object = new, genes.use = unique(markers$gene), slim.col.label = TRUE,
#   remove.key = TRUE, group.label.rot = TRUE
# )

singlets <- singlets[, colnames(singlets) %in% colnames(new@data)]
singletERCC <- singletERCC[, colnames(singletERCC) %in% colnames(singlets)]
idx <- match(rownames(FetchData(new, "ident")), colnames(singlets))
classes <- as.character(FetchData(new, "ident")[[1]])[idx]
names(classes) <- rownames(FetchData(new, "ident"))[idx]
var.genes <- unique(markers$gene)
select <- which(rownames(singlets) %in% var.genes)
dim.red <- new@dr$umap@cell.embeddings
colnames(dim.red) <- NULL

#setup CIMseqData objects
cObjSng <- CIMseqSinglets(singlets, singletERCC, dim.red, classes)
cObjMul <- CIMseqMultiplets(multiplets, multipletERCC, select)

#save
if(!"data" %in% list.dirs(currPath, full.names = FALSE)) system('mkdir data')
print(paste0("saving data to ", currPath, "."))
save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spUnsupervised.txt"))