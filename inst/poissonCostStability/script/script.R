#PACKAGES
packages <- c(
  "sp.scRNAseq", "sp.scRNAseqData", "printr", "ggthemes", "tidyverse", "future"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

##DATA
s <- grepl("^s", colnames(countsSorted2))
cObjSng <- spCounts(countsSorted2[,s], countsSortedERCC2[,s])
cObjMul <- spCounts(countsSorted2[, !s][, 1:2], countsSortedERCC2[, !s][, 1:2])

#spUnsupervised
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

#spSwarm
plan(multiprocess)
sObj <- spSwarm(
  cObjSng, cObjMul, uObj, maxiter = 10, 
  swarmsize = 150, nSyntheticMultiplets = 200
)
save(sObj, file = file.path(currPath, "data/sObj.rda"))

##TEST
singlets <- getData(cObjSng, "counts.cpm")[getData(uObj, "selectInd"), ]
oneMultiplet <- getData(cObjMul, "counts.cpm")[getData(uObj, "selectInd"), 1]
fractions <- unlist(getData(sObj, "spSwarm")[1, ])
classes <- getData(uObj, "classification")

tests <- seq(50, 1000, 50)
reps <- 50
data <- map_dfc(tests, function(x) {
  tibble(
    Cost = map_dbl(1:reps, function(n) {
      singletSubset <- sp.scRNAseq:::.subsetSinglets(classes, singlets, x)
      sp.scRNAseq::calculateCost(ceiling(oneMultiplet), singletSubset, fractions, x)
    })
  )
}) 

save(data, file = file.path(currPath, "data/costs.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))
