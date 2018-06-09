#PACKAGES
packages <- c(
  "sp.scRNAseq",
  "sp.scRNAseqData",
  "sp.scRNAseqTesting",
  "printr",
  "ggthemes",
  "tidyverse"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

s <- str_detect(colnames(countsMgfp), "^s")
commonGenes <- intersect(rownames(countsMgfp), rownames(countsRegev))

sng <- cbind(countsMgfp[commonGenes, s], countsRegev[commonGenes, ])
mul <- countsMgfp[commonGenes, !s]

erccSng <- cbind(
  countsMgfpERCC[, s], 
  matrix(NA, nrow = nrow(countsMgfpERCC), ncol = ncol(countsRegev))
)
erccMul <- cbind(countsMgfpERCC[, !s])

boolMulC <- colnames(mul) %in% filter(countsMgfpMeta, tissue == "colon")$sample
boolSngC <- colnames(sng) %in% filter(countsMgfpMeta, tissue == "colon")$sample

#setup spCounts
cObjSng <- spCounts(sng, erccSng)
cObjMul <- spCounts(mul, erccMul)

uObj <- spUnsupervised(cObjSng, max_iter = 7000, initial_dims = ncol(sng), seed = 87689)
save(uObj, file = "./uObjMouse.rda")

sObj <- spSwarm(
  cObjSng, cObjMul, uObj, maxiter = 100, swarmsize = 500, 
  nSyntheticMultiplets = 500, cores = 16
)
save(sObj, file = "./newAlgoSpSwarmMouse.rda")
