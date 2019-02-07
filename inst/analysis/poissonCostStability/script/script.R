#PACKAGES
packages <- c("CIMseq", "CIMseq.data", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

##DATA
if(file.exists('../SCM.analysis/data/CIMseqData.rda')) {
  load('../SCM.analysis/data/CIMseqData.rda')
}
if(file.exists('../SCM.analysis/data/sObj.rda')) {
  load('../SCM.analysis/data/sObj.rda')
}

##TEST
features <- getData(cObjMul, "features")
singlets <- getData(cObjSng, "counts.cpm")[features, ]
oneMultiplet <- getData(cObjMul, "counts.cpm")[features, 1]
fractions <- unlist(getData(sObj, "fractions")[1, ])
classes <- getData(cObjSng, "classification")

nSyntheticMultiplets <- seq(50, 1000, 50)
reps <- 50
data <- map_dfc(nSyntheticMultiplets, function(x) {
  tibble(
    Cost = map_dbl(1:reps, function(n) {
      singletIdx <- purrr::map(1:x, ~sampleSinglets(classes))
      singletSubset <- appropriateSinglets(cObjSng, singletIdx, features)
      CIMseq::calculateCost(ceiling(oneMultiplet), singletSubset, fractions, x)
    })
  )
}) 

save(data, file = file.path(currPath, "data/costs.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))