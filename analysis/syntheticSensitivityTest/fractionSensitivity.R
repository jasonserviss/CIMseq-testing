library(CIMseq)
library(tidyverse)
load('../SCM.analysis/CIMseqData.rda')

#generate sets of singlets where there are 2 cell types in each set and
#where in each set the cell types are progressivley more similar to one another.

#subset the two cell types to use as a start point
cpm <- getData(cObjSng, "counts.cpm")
classes <- getData(cObjSng, "classification")
ct1 <- cpm[, classes == 0]
ct2 <- cpm[, classes == 1]

#select a subset of the genes that differ between the two cell types
m.ct1 <- matrixStats::rowMeans2(ct1)
m.ct2 <- matrixStats::rowMeans2(ct2)
fc <- log2((m.ct1 / m.ct2) + 1)
genes <- rownames(cpm)[(fc > 2 | fc < -2) & !is.infinite(fc) & !is.nan(fc)]
select <- which(rownames(cpm) %in% genes)

#synthesize sets
reps <- 30
proportions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
idx <- map(c("0", "1"), function(c) {
  which(classes == c)
})

seed <- 8230
set <- map(proportions, function(p) {
  map_dfc(1:reps, function(i) {
    set.seed(seed + i)
    mat <- matrix(c(cpm[select, sample(idx[[1]], 1)], cpm[select, sample(idx[[2]], 1)]), ncol = 2)
    fractions <- c(p, 1-p)
    adj <- adjustAccordingToFractions(fractions, mat)
    m <- multipletSums(adj)
    r <- matrix(
      rpois(n = nrow(m), lambda = c(m)),
      dimnames = list(rownames(m.ct1), paste("ct2", p, i, sep = "-"))
    )
    as.data.frame(r)
  }) %>%
    as.matrix()
})

#calculate pearsons correlation
cor.ranges <- map(set, function(s) {
  range(map_dbl(1:ncol(s), ~cor(s[, .x], m.ct2[select], method = "p")))
})
c.ranges <- unlist(map(cor.ranges, ~paste(round(.x, digits = 2), collapse = "-")))
names(set) <- c.ranges

#plot correlations
total.set <- cbind(ct2[select, ], reduce(set, cbind))
cplot <- cor(total.set, method = "pearson")
cplot %>%
  matrix_to_tibble("from") %>%
  gather(to, cor, -from) %>%
  filter(str_detect(from, "NJB")) %>%
  filter(!str_detect(to, "NJB")) %>%
  separate(to, into = c("ct", "prob", "rep"), sep = "-", remove = FALSE) %>%
  mutate(prob = case_when(
    prob == proportions[1] ~ c.ranges[1],
    prob == proportions[2] ~ c.ranges[2],
    prob == proportions[3] ~ c.ranges[3],
    prob == proportions[4] ~ c.ranges[4],
    prob == proportions[5] ~ c.ranges[5],
    prob == proportions[6] ~ c.ranges[6],
    prob == proportions[7] ~ c.ranges[7],
    prob == proportions[8] ~ c.ranges[8],
    prob == proportions[9] ~ c.ranges[9],
    prob == proportions[10] ~ c.ranges[10]
  )) %>%
  ggplot() +
  geom_tile(aes(from, to, fill = cor)) +
  scale_fill_viridis_c() +
  facet_wrap(~prob, scales = "free_y") +
  theme(axis.text = element_blank()) +
  labs(x = "Cell type 1", y = "Cell type 2")

#show t-SNE
t <- RNAseqFunctions::runTsne(1-cplot)
rownames(t) <- colnames(cplot)
t %>%
  matrix_to_tibble("sample") %>%
  separate(sample, into = c("name", "prob", "rep"), sep = "-", remove = FALSE, fill = "right") %>%
  mutate(prob = if_else(str_detect(sample, "NJB"), "real", prob)) %>%
  ggplot() +
  geom_point(aes(V1, V2, colour = prob)) +
  scale_colour_manual(values = col40())


#combine into multiplets
reps <- 10
proportions <- c(0.001, 0.005, 0.01, 0.05, 0.1)
seed <- 8230
mul.set <- map(1:length(set), function(u) {
  s <- set[[u]]
  set.name <- names(set)[u]
  map_dfc(proportions, function(p) {
    map_dfc(1:reps, function(i) {
      set.seed(seed + i)
      mat <- matrix(c(s[, sample(1:ncol(s), 1)], ct2[select, sample(1:ncol(ct2), 1)]), ncol = 2)
      fractions <- c(p, 1-p)
      adj <- adjustAccordingToFractions(fractions, mat)
      m <- multipletSums(adj)
      r <- matrix(
        rpois(n = nrow(m), lambda = c(m)),
        dimnames = list(rownames(m.ct1), paste(set.name, p, i, sep = "_"))
      )
      as.data.frame(r)
    })
  }) %>%
    as.matrix()
})

#setup CIMseq classes
sng.objs <- map(set, function(s) {
  CIMseqSinglets(
    cbind(ct2[select, ], s), 
    matrix(),
    matrix(),
    c(rep("A", ncol(ct2)), rep("B", ncol(s)))
  )
})

mul.objs <- map(mul.set, function(s) {
  CIMseqMultiplets(s, matrix(), 1:length(select))
})

save(sng.objs, mul.objs, file = "data/CIMseqData.rda")

#Deconvolute
library(future.apply)
sObjs <- map2(sng.objs, mul.objs, function(s, m) {
  plan(multiprocess)
  CIMseqSwarm(s, m, maxiter = 2, swarmsize = 10, nSyntheticMultiplets = 10)
})

#analyze. start with set 1: pearson = 0.66-0.88
fractions <- map(sObjs, ~getData(.x, "fractions"))
res <- map_dbl(fractions, function(f) {
  d <- f %>%
    matrix_to_tibble("info") %>%
    separate(info, into = c("pearson", "expected.fraction.B", "rep"), sep = "_") %>%
    mutate(expected.fraction.B = as.numeric(expected.fraction.B)) %>%
    mutate(expected.fraction.A = 1 - expected.fraction.B) %>%
    mutate(
      diff.A = expected.fraction.A - A,
      diff.B = expected.fraction.B - B
    )
    mean(c(pull(d, diff.A), pull(d, diff.B)))
})

tibble(
  cor = names(res),
  diff = res
) %>%
  mutate(cor = parse_factor(cor, levels = c(
    "0.46-0.65", "0.48-0.66", "0.49-0.67", "0.51-0.68", "0.54-0.71", 
    "0.58-0.74", "0.62-0.78", "0.67-0.83", "0.69-0.88", "0.66-0.88"
  ))) %>%
  ggplot() +
  geom_point(aes(cor, res)) +
  labs(y = "mean(expected fraction - detected fraction)", x = "Pearsons correlation")

