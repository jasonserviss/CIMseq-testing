
#PACKAGES
packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

##spCounts
s <- grepl("^s", colnames(countsSorted2))
cObjSng <- spCounts(countsSorted2[,s], countsSortedERCC2[,s])
cObjMul <- spCounts(countsSorted2[,!s], countsSortedERCC2[,!s])

##spUnsupervised
uObj <- spUnsupervised(cObjSng)

#rename classes
midx <- match(rownames(getData(uObj, "tsne")), countsSortedMeta2$sample)
classification(uObj) <- countsSortedMeta2$cellTypes[midx]
groupMeans(uObj) <- averageGroupExpression(
  getData(cObjSng, "counts.cpm"), getData(uObj, "classification"), FALSE
)
tsneMeans(uObj) <- tsneGroupMeans(
  getData(uObj, "tsne"), getData(uObj, "classification")
)

save(uObj, file = "./data/uObj.rda")

##spSwarm
future.batchtools::plan(
  batchtools_slurm,
  template = "/crex/proj/snic2018-8-151/private/batchtools.slurm.tmpl",
  resources = list(
    account = "snic2018-8-151", partition = "core", ntasks = 1L, 
    time = "00:00:10", jobname = "testingPoissonSorted",
    modules = "R_packages/3.4.3", R = "R/3.4.3"
  )
)

sObj <- spSwarm(cObjSng, cObjMul, uObj, maxiter = 100, swarmsize = 500, nSyntheticMultiplets = 1000)
save(sObj, file = "./data/sObj.rda")

writeLines(capture.output(sessionInfo()), "./data/sessionInfo.txt")

