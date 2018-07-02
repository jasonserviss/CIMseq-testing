
#PACKAGES
packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse", "future", "future.apply")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

##spCounts
s <- grepl("^s", colnames(countsSorted2))
cObjSng <- spCounts(countsSorted2[, s], countsSortedERCC2[, s])
cObjMul <- spCounts(countsSorted2[, !s], countsSortedERCC2[, !s])

print("spCounts done")

##spUnsupervised
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

save(uObj, file = file.path(currPath, "data/uObj.rda"))
print("spUnsupervised done")

##spSwarm

tests <- seq(50, 1000, 100)
results <- map(tests, function(n) {
  future::plan(multiprocess)
  print(paste0("Starting deconvolution at ", Sys.time()))
  sObj <- spSwarm(cObjSng, cObjMul, uObj, maxiter = 10, swarmsize = 150, nSyntheticMultiplets = n)
  print(paste0("Finished deconvolution at ", Sys.time()))
  return(sObj)
})

save(results, file = file.path(currPath, "data/results.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))
print("finished")
