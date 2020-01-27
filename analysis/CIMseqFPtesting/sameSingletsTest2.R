
packages <- c("CIMseq", "sp.scRNAseqData", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

load('~/Github/sp.scRNAseqTesting/analysis/SCM.analysis/data/CIMseqData.rda')

#extract sample names for HCT116 singlets
s <- filter(
  SCM.Meta, 
  cellNumber == "Singlet" & cellTypes == "HCT116" & !filtered
)$sample

#setup 3 random classes for the HCT116 singlets 
#there are 27 HCT116 singlets in the dataset
class <- rep(LETTERS[1:3], each = length(s) / 3)

#make CIMseqSinglets object
cObjSng2 <- CIMseqSinglets(
  getData(cObjSng, "counts")[, s],
  getData(cObjSng, "counts.ercc")[, s],
  getData(cObjSng, "dim.red")[s, ],
  class
)

#setup CIMseqMultiplets object using HCT116-HOS doublets
s <- filter(
  SCM.Meta,
  cellNumber == "Multiplet" & cellTypes == "HCT116-HOS" & !filtered
)$sample

cObjMul2 <- CIMseqMultiplets(
  getData(cObjMul, "counts")[, s],
  matrix(), 
  getData(cObjMul, "features")
)

#deconvolute
plan(multiprocess)
sObj2 <- CIMseqSwarm(
  cObjSng2, cObjMul2,  maxiter = 100, swarmsize = 500,
  nSyntheticMultiplets = 400
)

#save
save(cObjSng2, cObjMul2, sObj2, file = './output2.rda')
