packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse", "e1071", "caret", "matrixStats", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

##FUNCTIONS
nTopDeltaCV <- function(counts, n) {
  valid <- matrixStats::rowSums2(counts) > 0
  mu <- matrixStats::rowMeans2(counts)
  sd <- matrixStats::rowSds(counts)
  ok <- mu > 0 & sd > 0
  cv <- sd[ok] / mu[ok]
  
  log2_m <- log2(mu[ok])
  log2_cv <- log2(cv)
  
  svr_gamma <- 1000 / length(mu[ok])
  modelsvm <- svm(log2_cv ~ log2_m, gamma = svr_gamma)
  score <- log2_cv - predict(modelsvm, log2_m)
  score <- score * valid[ok]
  names(score) <- rownames(counts)[ok]
  s <- sort(score, decreasing = TRUE)[1:n]
  which(rownames(counts) %in% names(s))
}

recursiveFE <- function(cpm, classes, nCV = 10, s) {
  normalized <- scale(cpm, center = TRUE, scale = TRUE)
  control <- rfeControl(functions = rfFuncs, method = "cv", number = nCV)
  results <- rfe(
    as.data.frame(t(normalized)),
    as.factor(classes),
    sizes = s,
    rfeControl = control
  )
  p <- predictors(results)
  which(rownames(cpm) %in% p)
}

findHighlyCorrelated <- function(cpm, cut) {
  sd <- matrixStats::rowSds(cpm)
  c <- cpm[sd != 0, ]
  correlationMatrix <- cor(t(c))
  hc <- findCorrelation(correlationMatrix, cutoff = cut)
  n <- 1:nrow(cpm)
  n[!n %in% hc]
}

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

#setup tests
selectList <- list(
  var = spTopVar(cObjSng, 2000),
  max = spTopMax(cObjSng, 2000),
  nTopDeltaCV = nTopDeltaCV(getData(cObjSng, "counts.cpm"), 2000),
  recursiveFE = recursiveFE(
    cpm = getData(cObjSng, "counts.cpm"), 
    classes = getData(uObj, "classification"), 
    nCV = 10, 
    s = seq(100, nrow(getData(cObjSng, "counts")), 1000)
  ),
  nonCorr = findHighlyCorrelated(getData(cObjSng, "counts.cpm"), 0.5)
)

#look at overlap
cmb <- combn(names(selectList), 2)
int <- apply(cmb, 2, function(x) {
  length(intersect(selectList[[x[1]]], selectList[[x[2]]]))
})
names(int) <- apply(cmb, 2, paste, collapse = "-")

##spSwarm
print(paste0("Starting deconvolution at ", Sys.time()))
sObjs <- purrr::map(selectList, function(x) {
  future::plan(multiprocess)
  spSwarm(
    cObjSng, cObjMul, uObj, maxiter = 10, swarmsize = 150, 
    nSyntheticMultiplets = 400, selectInd = x
  )
})
print(paste0("Finished deconvolution at ", Sys.time()))

save(sObjs, file = file.path(currPath, "data/sObjs.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))
print("finished")

