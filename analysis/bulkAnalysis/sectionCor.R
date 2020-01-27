library(sp.scRNAseqData)
library(tidyverse)
si <- MGAB.Counts[, 1:8]
si <- si[matrixStats::rowSums2(si) > 0, ]
si <- RNAseqFunctions::cpm(si)
si <- log2(si + 1)
# DGE <- edgeR::DGEList(si)
# DGE <- edgeR::calcNormFactors(DGE, method = c("TMM"))
# v <- limma::voom(DGE)
# si <- as.matrix(v)

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

lr <- map(1:nrow(si), function(i) lm(si[i, ] ~ seq(1, 8, 1)))
cr <- tibble(
  gene = rownames(si),
  adj.r.squared = map_dbl(lr, function(i) summary(i)$adj.r.squared),
  coef = map_dbl(lr, function(m) coef(m)[2]),
  p = map_dbl(lr, lmp)
)

filter(cr, p < 0.05 & (coef > 0.5 | coef < -0.5)) %>%
  arrange(coef) %>% 
  pull(gene) %>%
  {si[match(., rownames(si)), ]} %>%
  matrix_to_tibble("gene") %>%
  gather(sample, cpm, -gene) %>%
  mutate(sample = parse_factor(sample, levels = paste("SI", 1:8, sep = "."))) %>%
  mutate(gene = parse_factor(gene, levels = unique(gene))) %>%
  ggplot() +
  geom_tile(aes(sample, gene, fill = cpm)) +
  scale_fill_viridis_c() +
  theme(axis.text.y = element_blank())


lr <- map(1:nrow(si), function(i) lm(si[i, ] ~ poly(seq(1, 8, 1), 2)))
cr <- tibble(
  gene = rownames(si),
  adj.r.squared = map_dbl(lr, function(i) summary(i)$adj.r.squared),
  coef = map_dbl(lr, function(m) coef(m)[2]),
  p = map_dbl(lr, lmp)
)

filter(cr, p < 0.05 & (coef > 0.7 | coef < -0.7)) %>%
  arrange(coef) %>% 
  pull(gene) %>%
  {si[match(., rownames(si)), ]} %>%
  matrix_to_tibble("gene") %>%
  gather(sample, cpm, -gene) %>%
  mutate(sample = parse_factor(sample, levels = paste("SI", 1:8, sep = "."))) %>%
  mutate(gene = parse_factor(gene, levels = unique(gene))) %>%
  ggplot() +
  geom_tile(aes(sample, gene, fill = cpm)) +
  scale_fill_viridis_c() +
  theme(axis.text.y = element_blank())

lr <- map(1:nrow(si), function(i) nlm(si[i, ] ~ seq(1, 8, 1), p = 0))
cr <- tibble(
  gene = rownames(si),
  adj.r.squared = map_dbl(lr, function(i) summary(i)$adj.r.squared),
  coef = map_dbl(lr, function(m) coef(m)[2]),
  p = map_dbl(lr, lmp)
)

#Tmprss15, Pdx1, Nkx6-3
find <- function(cpm, index, sds) {
  m <- matrixStats::rowMeans2(cpm)
  s <- matrixStats::rowSds(cpm)
  cut.h <- m + (s * sds)
  cut.l <- m - (s * sds)
  diff <- cpm[, index] - m
  bool <- cpm[, index] > cut.h | cpm[, index] < cut.l
  diff[bool]
}

genes <- map_dfr(tmp, function(i) tibble(gene = names(i), diff = i), .id = "id") %>%
  mutate(col =col40()[as.numeric(id)])

gplots::heatmap.2(
  si[match(genes$gene, rownames(si)), ],
  RowSideColors = genes$col,
  Colv=FALSE, 
  dendrogram="row",
  trace="none",
  scale = "row",
  key = FALSE
)

time <- 1:8
size <- rep(1, nrow(si))
out <- timecourse::mb.long(si, times = 8, reps = size)
