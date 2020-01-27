library(tidyverse)
library(CIMseq.data)

#get names of cells present post-filtering
matrix <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109774&format=file&file=GSE109774%5FColon%2Etar%2Egz"
target <- tempdir()
download.file(url = matrix, destfile = file.path(target, "colon.tar.gz"))
files <- basename(untar(file.path(target, "colon.tar.gz"), list = TRUE))
names <- tools::file_path_sans_ext(files)

#get metadata
link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109774/suppl//GSE109774_list_of_SRR_accessions_and_raw_filenames.txt.gz"
meta <- read_delim(link, col_names = FALSE, delim = "\t") %>%
  mutate(base = str_replace(X3, "(.*)_C._R1.*", "\\1"))
all(names %in% pull(meta, base))
colonID <- "GSM2967048"
srrid <- pull(filter(meta, base %in% names & X1 == colonID), X2)

#plot
colon <- filter(meta, base %in% names & X1 == colonID)
tmp2 <- colon %>% count(base) %>% filter(n == 2) %>% pull(base)
tmp3 <- filter(tmp, base %in% tmp2) %>% pull(X2)
TMD.Meta %>% 
  mutate(sample = str_replace(sample, "..(.*)", "\\1")) %>% 
  inner_join(colon, by = c("sample" = "X2")) %>% 
  mutate(tmp = str_replace(X3, ".*(C.)_R1.fastq.gz", "\\1")) %>% 
  inner_join(tibble(sample = str_replace(colnames(TMD.Counts), "..(.*)", "\\1"), cs = colSums(TMD.Counts))) %>% 
  filter(sample %in% tmp3) %>%
  ggplot() +
  geom_boxplot(aes(tmp, cs)) 

TMD.Meta %>% 
  mutate(sample = str_replace(sample, "..(.*)", "\\1")) %>% 
  inner_join(colon, by = c("sample" = "X2")) %>% 
  mutate(tmp = str_replace(X3, ".*(C.)_R1.fastq.gz", "\\1")) %>% 
  inner_join(tibble(sample = str_replace(colnames(TMD.Counts), "..(.*)", "\\1"), cs = colSums(TMD.Counts))) %>% 
  filter(tmp == "C1") %>%
  mutate(c2 = if_else(sample %in% tmp3, "C2", "no C2")) %>%
  ggplot() + 
  geom_boxplot(aes(c2, cs))

#merge
pairs <- colon %>%
  group_by(base) %>%
  summarize(srr = paste(X2, collapse = "__", sep = "__")) %>%
  mutate(srr = map(srr, ~str_split(.x, "__")[[1]])) %>%
  pull(srr)

projectName <- "Tabula.muris.data_TMD"
countData <- getCountsData(projectName)
countData <- moveGenesToRownames(countData)
colnames(countData) <- str_replace(colnames(countData), "(.*)\\..*", "\\1")
ercc <- detectERCCreads(countData)
countsERCC <- countData[ercc, ]
counts <- countData[!ercc, ]
counts <- countData[!detectNonGenes(counts), ]
counts <- as.matrix(counts)

test <- sapply(pairs, function(p) {
  print(p)
  if(length(p) == 1) {
    if(p %in% colnames(counts)) {
      counts[, p]
    } else {
      rep(0, 24514)
    }
    
  } else if(all(p %in% colnames(counts))) {
    matrixStats::rowSums2(counts[, p])
  } else {
    counts[, p[p %in% colnames(counts)]]
  }
})

test <- test[, colSums(test) > 1e4]
test <- test[, apply(test, 2, function(c) length(which(c > 0)) > 500)]
colnames(test) <- paste0("SRR654", 1:ncol(test))

#dim red
cpm <- CIMseq:::.norm.counts(test)
lcpm <- CIMseq:::.norm.log.counts(test)
features <- CIMseq::selectTopMax(cpm, 2000)
p.cor <- CIMseq::pearsonsDist(lcpm, features)
tsne <- CIMseq::runTsne(p.cor)

genes <- c("Lgr5", "Ptprc", "Chga", "Dclk1", "Slc26a3", "Atoh1", "Alpi", "Lyz1", "Mki67", "Hoxb13")
CIMseq::matrix_to_tibble(tsne, "sample") %>%
  inner_join(matrix_to_tibble(t(lcpm[genes, ]), "sample")) %>%
  gather(gene, lcpm, -(sample:V2)) %>%
  ggplot() +
  geom_point(aes(V1, V2, colour = lcpm), size = 0.5) +
  scale_colour_viridis_c() + 
  facet_wrap(~gene) +
  ggthemes::theme_few()

#merge with others and dim red
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
iGenes <- intersect(intersect(rownames(RSI.Counts), rownames(MGA.Counts)), rownames(test))
singlets <- cbind(MGA.Counts[iGenes, boolSng], RSI.Counts[iGenes, boolRSI], test[iGenes, ])

library(Seurat)
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
pcs <- 1:30
#PCElbowPlot(object = mca, num.pc = 100) + scale_x_continuous(breaks = seq(0, 100, 5))

mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.5,
  n_neighbors = 20, seed.use = 79354
)

#mca <- RunTSNE(mca, dims.use = pcs, gene.use = mca@var.genes, seed.use = 78239)

mca <- FindClusters(
  object = mca, reduction.type = "harmony", dims.use = pcs, resolution = 2.3,
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
# 
# FeaturePlot(
#   mca,
#   c("Lgr5", "Ptprc", "Chga", "Dclk1", "Slc26a3", "Atoh1", "Lyz1", "Alpi", "Hoxb13", "Mki67"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 0.1,
#   vector.friendly = FALSE
# )

