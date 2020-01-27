
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

#gradually change the gene expression in the features of 2 of the HCT116 singlets classes
#to determine when the algo can differentiate between them
features <- getData(cObjMul, "features")
change <- getData(cObjSng, "counts")[, c(s, s2)]
changed <- expand.grid(
  nGenesToChange = c(1, 10, 100, 1000),
  degreeToChange = c(1.1, 1.25, 1.5, 2)
) %>%
  as_tibble() %>%
  mutate(counts = map2(nGenesToChange, degreeToChange, function(g, d) {
    idx <- features[sample(1:length(features), g)]
    change[idx, s] <- change[idx, s] * d
    change
  }))

#make CIMseqSinglets object
cObjSngs <- map(changed, function(c) {
  CIMseqSinglets(
    c,
    getData(cObjSng, "counts.ercc")[, c(s, s2)],
    getData(cObjSng, "dim.red")[c(s, s2), ],
    class
  )
})

#setup CIMseqMultiplets object using HCT116-HOS doublets
s <- filter(
  SCM.Meta,
  cellNumber == "Multiplet" & cellTypes == "HCT116-HOS" & !filtered
)$sample[1:8]

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
save(cObjSng2, cObjMul2, sObj2, file = './output3.rda')
