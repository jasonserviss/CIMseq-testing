
packages <- c("CIMseq", "sp.scRNAseqData", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

load('~/Github/sp.scRNAseqTesting/analysis/SCM.analysis/data/CIMseqData.rda')

#extract sample names for HCT116 singlets
s <- filter(
  SCM.Meta, 
  cellNumber == "Singlet" & cellTypes == "HCT116" & !filtered
)$sample

#extract sample names for HOS singlets
s2 <- filter(
    SCM.Meta,
    cellNumber == "Singlet" & cellTypes == "HOS" & !filtered
)$sample

#setup 3 random classes for the HCT116 singlets 
#there are 27 HCT116 singlets in the dataset
class <- rep(LETTERS[1:3], each = length(s) / 3)
class <- c(class, rep("HOS", length = length(s2)))

#make CIMseqSinglets object
cObjSng2 <- CIMseqSinglets(
  getData(cObjSng, "counts")[, c(s, s2)],
  getData(cObjSng, "counts.ercc")[, c(s, s2)],
  getData(cObjSng, "dim.red")[c(s, s2), ],
  class
)

#only include one gene with all 0 in HCT116 and some != 0 in HOS
idx <- apply(getData(cObjSng, "counts")[, c(s, s2)], 1, function(r) {
  zero <- sum(r[class %in% LETTERS[1:3]])
  notzero <- sum(r[class == "HOS"])
  if_else(zero == 0 & notzero != 0, TRUE, FALSE)
})
select <- as.integer(which(idx)[1:2] - 1) #0-based indices

#setup CIMseqMultiplets object using HCT116-HOS doublets
s <- filter(
  SCM.Meta,
  cellNumber == "Multiplet" & cellTypes == "HCT116-HOS" & !filtered
)$sample

cObjMul2 <- CIMseqMultiplets(
  getData(cObjMul, "counts")[, s],
  matrix(), 
  select
)

#deconvolute
plan(multiprocess)
sObj2 <- CIMseqSwarm(
  cObjSng2, cObjMul2,  maxiter = 100, swarmsize = 500,
  nSyntheticMultiplets = 400
)

#save
save(cObjSng2, cObjMul2, sObj2, file = './output4.rda')
