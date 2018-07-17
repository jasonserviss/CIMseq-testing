
packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse", "future", "future.apply")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

##spCounts
s <- str_detect(colnames(countsMgfp), "^s")
commonGenes <- intersect(rownames(countsMgfp), rownames(countsRegev))

sng <- cbind(countsMgfp[commonGenes, s], countsRegev[commonGenes, ])
mul <- countsMgfp[commonGenes, !s]

erccSng <- cbind(
  countsMgfpERCC[, s], 
  matrix(NA, nrow = nrow(countsMgfpERCC), ncol = ncol(countsRegev))
)
erccMul <- cbind(countsMgfpERCC[, !s])

#setup spCounts
cObjSng <- spCounts(sng, erccSng)
#cObjMul <- spCounts(mul, erccMul)

testSamples <- c(
  "m.NJA00107.G09", "m.NJA00107.D12", "m.NJA00107.A02",
  "m.NJA00107.A10", "m.NJA00107.C08"
)
cObjMul <- spCounts(mul[, testSamples], erccMul[, testSamples])

load('../testingPoissonMouse/data/uObj.rda')
load('../testingPoissonMouse/data/sObj.rda')

nPerms <- 10
perms <- map(1:nPerms, function(x) {
  plan(multiprocess)
  spSwarm(
    cObjSng, cObjMul, uObj, maxiter = 10, swarmsize = 150,
    nSyntheticMultiplets = 400, permute = TRUE
  )
})

save(perms, file = file.path(currPath, "permutations.rda"))

