packages <- c("CIMseq", "sp.scRNAseqData", "tidyverse", "Seurat")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

##########################!!!!

#The script is incomplete at this time

##########################!!!!

#check package version
algoV <- sessionInfo()$otherPkgs$CIMseq$Version
last3 <- paste(strsplit(algoV, "\\.")[[1]][2:4], collapse = "")
if(!as.numeric(last3) >= 100) {
  stop("sp.scRNAseq package version too low. Must be >= 0.2.0.0")
}

currPath <- getwd()

#only tumor cells
st <- str_detect(colnames(ADM.Counts), "^s")
singlets <- ADM.Counts[, st]
singletsERCC <- ADM.CountsERCC[, st]
multiplets <- ADM.Counts[, !st]
multipletsERCC <- ADM.CountsERCC[, !st]
meta <- filter(ADM.Meta, sample %in% colnames(singlets))

#Dimensionality reduction and classification
print(paste0("Starting only tumor analysis at ", Sys.time()))
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
  object = mca, pc.genes = mca@var.genes, do.print = FALSE, pcs.compute = 100
)
mca <- JackStraw(object = mca, num.replicate = 100, display.progress = TRUE, num.pc = 50)
mca <- JackStrawPlot(object = mca, PCs = 1:50)
PCp <- mca@dr$pca@jackstraw@overall.p.values
pcs <- PCp[PCp[, 2] < 0.001, 1]
#PCElbowPlot(object = mca, num.pc = 100) + scale_x_continuous(breaks = seq(0, 100, 5))

set.seed(7239832)
mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.3,
  n_neighbors = 5
)
mca <- FindClusters(
  object = mca, reduction.type = "pca", dims.use = pcs, resolution = 1,
  save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE,
  force.recalc = TRUE
)

DimPlot(
  object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
  vector.friendly = FALSE, pt.size = 1
) + scale_colour_manual(values = col40())

FeaturePlot(
  mca,
  c(c("Lgr5", "Ptprc", "Atoh1", "Lyz1", "Hoxb13", "Mki67"), c("Tbx21", "Itgam", "Arg2", "Cd9")),
  reduction.use = "umap", dark.theme = FALSE, pt.size = 1,
  vector.friendly = FALSE
)
print(paste0("Done only tumor analysis at ", Sys.time()))

#all cells
sn <- str_detect(colnames(MGA.Counts), "^s")
commonGenes <- intersect(rownames(MGA.Counts), rownames(ADM.Counts))
sng <- cbind(MGA.Counts[commonGenes, sn], ADM.Counts[commonGenes, st])
erccSng <- cbind(MGA.CountsERCC[, sn], ADM.CountsERCC[, st])

meta = tibble(origin = if_else(str_detect(colnames(sng), "^NJA"), "normal", "tumor"))


#save
save(uObj_tumor, uObj_all, file = file.path(currPath, "data/uObjs.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))
