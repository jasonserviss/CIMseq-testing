
packages <- c("sp.scRNAseq", "sp.scRNAseqData", "purrr", "future", "future.apply")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

##spCounts
s <- grepl("^s", colnames(countsSorted2))
cObjSng <- spCounts(countsSorted2[, s], countsSortedERCC2[, s])
#cObjMul <- spCounts(countsSorted2[, !s], countsSortedERCC2[, !s])
cObjMul <- spCounts(countsSorted2[, c("m.NJB00204.C04", "m.NJB00204.E04")], countsSortedERCC2[, c("m.NJB00204.C04", "m.NJB00204.E04")])

load('../testingPoissonSorted/data/uObj.rda')
load('../testingPoissonSorted/data/sObj.rda')

nPerms <- 10
perms <- map(1:nPerms, function(x) {
  plan(multiprocess)
  spSwarm(
    cObjSng, cObjMul, uObj, maxiter = 10, swarmsize = 150,
    nSyntheticMultiplets = 400, permute = TRUE
  )
})



