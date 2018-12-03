#PACKAGES
packages <- c("CIMseq", "sp.scRNAseqData", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
if(file.exists(file.path(currPath, '../MGA.analysis_SI/data/CIMseqData.rda'))) {
  load(file.path(currPath, '../MGA.analysis_SI/data/CIMseqData.rda'))
}

cObjMul <- CIMseqMultiplets(
 getData(cObjMul, "counts")[, 1:2],
 getData(cObjMul, "counts.ercc")[, 1:2],
 getData(cObjMul, "features")
)

calculateCentroids <- function(mat, select = NULL) {
  classes <- getData(cObjSng, "classification")
  if(is.null(select)) select <- 1:nrow(mat)
  mat <- mat[select, ]
  map_dfc(unique(classes), function(c) {
    tibble(mean = matrixStats::rowMeans2(mat[, classes == c])) %>%
      setNames(c)
  }) %>%
    add_column(gene = rownames(mat), .before = 1)
}

calculateCentroidKNN <- function(mat, select, k) {
  classes <- getData(cObjSng, "classification")
  centroid <- calculateCentroids(mat, select) %>%
    select(-gene) %>%
    as.matrix()
  if(is.null(select)) select <- 1:nrow(mat)
  if(length(k) == 1) k <- rep(k, ncol(centroid))
  mat <- mat[select, ]
  map(1:ncol(centroid), function(i) {
    c <- colnames(centroid)[i]
    knn <- nn2(t(mat[, classes == c]), matrix(centroid[, i], ncol = nrow(centroid)), k = k[i])[[1]]
    colnames(cpm)[classes == c][c(knn)]
  })
}

classes <- getData(cObjSng, "classification")
select <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 20)
cpm <- getData(cObjSng, "counts.cpm")
geneIdx <- getData(cObjMul, "features")
pca <- gmodels::fast.prcomp(t(cpm[geneIdx, ]), scale. = TRUE)$x
t <- as.numeric(table(getData(cObjSng, "classification"))[unique(classes)])
k <- if_else(t > 30, 30, t)
nn <- calculateCentroidKNN(t(pca), select, k)

cObjSng.hq <- cObjSng
keep <- unlist(nn)
counts <- getData(cObjSng.hq, "counts")
counts.ercc <- getData(cObjSng.hq, "counts.ercc")
dim.red <- getData(cObjSng.hq, "dim.red")
class <- getData(cObjSng.hq, "classification")
cObjSng.hq@counts <- counts[, colnames(counts) %in% keep]
cObjSng.hq@counts.ercc <- counts.ercc[, colnames(counts.ercc) %in% keep]
cObjSng.hq@dim.red <- dim.red[rownames(dim.red) %in% keep, ]
cObjSng.hq@classification <- class[colnames(counts) %in% keep]

#future::plan(multiprocess)
options(future.wait.interval = 10000.0)
options(future.wait.timeout = 1e9)
future::plan(
  future.batchtools::batchtools_slurm,
  template = "/crex/proj/snic2018-8-151/private/batchtools.slurm.tmpl",
  resources = list(
    account = "snic2018-8-151", partition = "core", ntasks = 1L,
    time = "24:00:00", jobname = "testingPoissonSorted",
    modules = "R_packages/3.5.0", R = "R/3.5.0", log.file = file.path(currPath, "logs/slurm.txt")
  ),
  workers = 100
)

#run deconvolution
print(paste0("Starting deconvolution at ", Sys.time()))
sObj <- CIMseqSwarm(
  cObjSng.hq, cObjMul, maxiter = 100, swarmsize = 500, nSyntheticMultiplets = 400
)
print(paste0("Finished deconvolution at ", Sys.time()))

save(cObjSng.hq, sObj, file = file.path(currPath, "data/output.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spSwarm.SI.txt"))
print("finished")
