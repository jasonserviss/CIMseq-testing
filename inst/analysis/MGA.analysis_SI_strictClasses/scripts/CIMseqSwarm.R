#ARGS
args <- commandArgs(TRUE)

#PACKAGES
packages <- c("CIMseq", "CIMseq.data", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#DATA
if(file.exists(file.path(currPath, '../MGA.analysis_SI/data/CIMseqData.rda'))) {
  load(file.path(currPath, '../MGA.analysis_SI/data/CIMseqData.rda'))
}

#FUNCTIONS
averageGroupExpression <- function(singlets) {
  classes <- getData(singlets, "classification")
  u.classes <- unique(classes)
  exp <- getData(singlets, "counts.cpm")
  
  sapply(u.classes, function(x) {
    matrixStats::rowMeans2(exp[, classes == x])
  })
}

medianGroupExpression <- function(singlets) {
  classes <- getData(singlets, "classification")
  u.classes <- unique(classes)
  exp <- getData(singlets, "counts.cpm")
  
  sapply(u.classes, function(x) {
    matrixStats::rowMedians(exp[, classes == x])
  })
}

distToCenter <- function(singlets, method) {
  if(method == "mean") {
    center <- averageGroupExpression(singlets)
  } else if(method == "median") {
    center <- medianGroupExpression(singlets)
  } else {
    stop("Enter a valid method.")
  }
  
  cpm <- getData(singlets, "counts.cpm")
  classes <- getData(singlets, "classification")
  
  dists <- sapply(1:ncol(cpm), function(i) {
    d <- as.matrix(dist(t(cbind(cpm[, i], center))))[, 1]
    unname(d[classes[i]])
  })
  tibble(
    sample = colnames(cpm),
    class = classes,
    distToCenter = dists
  )
}

pickBestClassRepresentation <- function(
  singlets, method, topPercent, minPerGroup = NULL, safe = TRUE
){
  filtered <- distToCenter(singlets, method) %>%
    group_by(class) %>%
    arrange(class, distToCenter) %>%
    filter(distToCenter < quantile(distToCenter, topPercent)) 
  
  if(safe & !is.null(minPerGroup)) {
    n <- filtered %>% summarize(n = n()) %>% pull(n)
    if(any(n < minPerGroup)) stop("Too few in group. Increase top percent.")
  }
  pull(filtered, sample)
}  

#SCRIPT
keep <- pickBestClassRepresentation(cObjSng, "median", 0.65, 5, TRUE)
keep.idx <- colnames(getData(cObjSng, "counts")) %in% keep
cObjSng.new <- CIMseqSinglets(
  counts = getData(cObjSng, "counts")[, keep.idx],
  counts.ercc = getData(cObjSng, "counts.ercc")[, keep.idx],
  dim.red = getData(cObjSng, "dim.red")[keep.idx, ],
  classification = getData(cObjSng, "classification")[keep.idx]
)

plan(sequential)
i <- as.integer(args[1])
print(paste0("Running multiplet ", i))
counts <- getData(cObjMul, "counts")
counts.ercc <- getData(cObjMul, "counts.ercc")
cObjMul <- CIMseqMultiplets(
  matrix(
    counts[, i], ncol = 1, 
    dimnames = list(rownames(counts), colnames(counts)[i])
  ),
  matrix(
    counts.ercc[, i], ncol = 1, 
    dimnames = list(rownames(counts.ercc), colnames(counts.ercc)[i])
  ),
  getData(cObjMul, "features")
)

#run deconvolution
print(paste0("Starting deconvolution at ", Sys.time()))
sObj <- CIMseqSwarm(
  cObjSng.new, cObjMul, maxiter = 100, swarmsize = 500, 
  nSyntheticMultiplets = 400
)
print(paste0("Finished deconvolution at ", Sys.time()))
save(sObj, file = file.path(currPath, paste0("tmp/sObj_", i, "_uppmax.rda")))

writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spSwarm.SI.txt"))
print("finished")
