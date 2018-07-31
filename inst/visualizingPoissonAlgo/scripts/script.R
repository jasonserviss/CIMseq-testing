packages <- c(
  "sp.scRNAseq","sp.scRNAseqData", "sp.scRNAseqTesting", "printr",
  "ggthemes", "tidyverse", "future"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

#############
#FUNCTIONS
############

syntheticMultipletsFromCounts <- function(
  singlets, classes, fractions,
  seed, idx, ...
){
  set.seed(seed)
  exCounts <- singlets[, idx]
  
  #adjust the counts by multiplying with the corresponding values adjustment
  adjusted <- round(t(t(exCounts) * fractions))
  
  #calculate the sum of counts for each gene
  rs <- matrixStats::rowSums2(adjusted)
  matrix(rs)
}

generateSyntheticMultiplets <- function(
  singlets, classes, fractions,
  seed, n, idx = NULL, ...
){
  mat <- purrr::map(1:n, function(y) {
    syntheticMultipletsFromCounts(
      singlets = singlets, classes = classes,
      fractions = fractions, seed = seed + y, idx = idx[[y]]
    )
  }) %>%
    do.call(cbind, .)
  
  colnames(mat) <- paste0("multiplet.", 1:n)
  rownames(mat) <- rownames(singlets)
  mat
}
############

##RUN ALGO
s <- grepl("^s", colnames(countsSorted2))
cObjSng <- spCounts(countsSorted2[,s], countsSortedERCC2[,s])
bool <- colnames(countsSorted2) %in% c("m.NJB00204.A02", "m.NJB00204.B01")
cObjMul <- spCounts(countsSorted2[,!s & bool], countsSortedERCC2[,!s & bool])
uObj <- spUnsupervised(cObjSng)

#rename classes
midx <- match(rownames(getData(uObj, "tsne")), countsSortedMeta2$sample)
classification(uObj) <- countsSortedMeta2$cellTypes[midx]
groupMeans(uObj) <- averageGroupExpression(
  cObjSng, getData(uObj, "classification"), FALSE
)
tsneMeans(uObj) <- tsneGroupMeans(
  getData(uObj, "tsne"), getData(uObj, "classification")
)

plan(sequential)
sObj <- spSwarm(cObjSng, cObjMul, uObj)

save(cObjSng, cObjMul, uObj, sObj, file = "data/algoOut.rda")

##GENERATE SYNTHETIC MULTIPLETS
#generate synthetic multiplets
g <- c("ANXA3", "CD74", "ACTG2")
class <- getData(uObj, "classification")

#A375-HCT116 multiplet: m.NJB00204.A02
m <- "m.NJB00204.A02"
mul <- getData(cObjMul, "counts.cpm")
oneMultiplet <- ceiling(mul[, colnames(mul) == m])
n.multiplets <- 1000
idx <- map(1:n.multiplets, function(n) {
  purrr::map_int(unique(class), ~sample(which(class == .x), 1))
})

#generate synthetic multiplets
#use correct fractions
adjust.right <- as.numeric(getData(sObj, "spSwarm")[m, unique(class)]) #A375-HCT116
synthetic.mul.right <- generateSyntheticMultiplets(
  getData(cObjSng, "counts.cpm"), class, fractions = adjust.right, 
  seed = 8923, n = n.multiplets, idx = idx
) %>%
  matrix_to_tibble("gene") %>%
  add_column(real.multiplet = oneMultiplet) %>%
  filter(gene %in% g) %>%
  gather(sample, value, -gene, -real.multiplet) %>%
  mutate(type = paste0(
    "Correct solution: A375 (", round(adjust.right[3], digits = 2), 
    "), HCT116 (", round(adjust.right[2], digits = 2), 
    "), HOS (", round(adjust.right[1], digits = 2), ")"
  ))

#use incorrect fractions
adjust.wrong <- rep(1/3, 3)
synthetic.mul.wrong <- generateSyntheticMultiplets(
  getData(cObjSng, "counts.cpm"), class, fractions = adjust.wrong, 
  seed = 8923, n = n.multiplets, idx = idx
) %>%
  matrix_to_tibble("gene") %>%
  add_column(real.multiplet = oneMultiplet) %>%
  filter(gene %in% g) %>%
  gather(sample, value, -gene, -real.multiplet) %>%
  mutate(type = paste0(
    "Incorrect solution: A375 (", round(adjust.wrong[3], digits = 2), 
    "), HCT116 (", round(adjust.wrong[2], digits = 2), 
    "), HOS (", round(adjust.wrong[1], digits = 2), ")"
  ))

#bind results into one structure
synthetic.mul <- bind_rows(synthetic.mul.right, synthetic.mul.wrong) %>%
  mutate(gene = case_when(
    gene == "ANXA3" ~ "HCT116: ANXA3",
    gene == "CD74" ~ "A375: CD74",
    gene == "ACTG2" ~ "HOS: ACTG2",
    TRUE ~ "error"
  ))

save(synthetic.mul, file = "data/syntheticMultiplets.rda")
writeLines(capture.output(sessionInfo()), "logs/sessionInfo.txt")
