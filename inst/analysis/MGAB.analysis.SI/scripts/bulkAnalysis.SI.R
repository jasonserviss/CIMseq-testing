#PACKAGES
packages <- c("CIMseq", "CIMseq.data", "tidyverse", "future.apply", "pso")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#FUNCTIONS
#Calculates counts per million
cpm <- CIMseq:::.norm.counts

#Calculates the mean expression of each gene in each class
averageGroupExpression <- function(exp, classes) {
  c <- unique(classes)
  means <- purrr::map_dfc(c, function(x) {
    data.frame(rowMeans(exp[, classes == x]))
  })
  means <- as.matrix(means)
  colnames(means) <- c
  return(means)
}

#calculates the sum for each gene after each cell type has been multiplied by
#fractions
.makeSyntheticSlice <- function(blueprint, fractions) {
  return(colSums(t(blueprint) * fractions))
}

#cost function
distToSlice <- function(fractions, blueprint, oneMultiplet) {
  if(sum(fractions) == 0) {
    return(999999999)
  }
  normFractions <- fractions / sum(fractions)
  a <- .makeSyntheticSlice(blueprint, normFractions)
  sum(abs(a - oneMultiplet))
}

distToSliceNorm <- function(fractions, blueprint, oneMultiplet) {
  if(sum(fractions) == 0) {
    return(999999999)
  }
  normFractions <- fractions / sum(fractions)
  blueprint <- blueprint / mean(blueprint)
  a = .makeSyntheticSlice(blueprint, normFractions)
  a <- a/mean(a)
  sum(abs((oneMultiplet - a) / (a+1)))
}

euclid <- function(fractions, blueprint, oneMultiplet) {
  if(sum(fractions) == 0) {
    return(999999999)
  }
  normFractions <- fractions / sum(fractions)
  a = .makeSyntheticSlice(blueprint, normFractions)
  dist(rbind(a, oneMultiplet))
}

pdist <- function(fractions, blueprint, oneMultiplet) {
  if(sum(fractions) == 0) {
    return(999999999)
  }
  normFractions <- fractions / sum(fractions)
  a = .makeSyntheticSlice(blueprint, normFractions)
  1 - cor(a, oneMultiplet, method = "pearson")
}

#optimization function for 1 individual sample to be deconvoluted
.optim.fn <- function(
  i, fractions, blueprint, multiplets, control
){
  psoptim(
    par = fractions, fn = distToSlice, blueprint = blueprint,
    oneMultiplet = multiplets[, i], lower = 0, upper = 1, control = control
  )
}

#wrapper for deconvoluting all samples
deconvoluteAll <- function(
  to, fractions, multiplets, blueprint, control
){
  future_lapply(
    X = 1:to, FUN = function(i) {
      .optim.fn(
        i, fractions = fractions, multiplets = multiplets,
        blueprint = blueprint, control = control
      )
    })
}

#SCRIPT

#setup fractions, multiplets, blueprint, and control variables
#load singlet data to be used as a blueprint
load('../MGA.analysis_SI/data/CIMseqData.rda')

features <- getData(cObjMul, "features")
shared.features <- intersect(
  rownames(MGAB.Counts), 
  rownames(getData(cObjSng, "counts"))[features]
)

blueprint.cpm <- getData(cObjSng, "counts.cpm")
blueprint <- averageGroupExpression(
  blueprint.cpm[which(rownames(blueprint.cpm) %in% shared.features), ], 
  getData(cObjSng, "classification")
)
blueprint <- blueprint[, order(colnames(blueprint))]

#setup fractions variable
nuc <- length(unique(getData(cObjSng, "classification")))
fractions <- rep(1.0 / nuc, nuc)

#setup multiplets
MGAB.Cpm <- cpm(MGAB.Counts[, filter(MGAB.Meta, sub_tissue == "small_intestine")$sample])

features.bulk <- which(rownames(MGAB.Cpm) %in% shared.features)
multiplets <- matrix(
  MGAB.Cpm[features.bulk, ],
  ncol = ncol(MGAB.Cpm),
  dimnames = list(rownames(MGAB.Cpm)[features.bulk], colnames(MGAB.Cpm))
)

#setup control
to <- dim(multiplets)[2]
control <- list(maxit = 100, s = 500)

#Run deconvolution
plan(multiprocess)
opt.out <- deconvoluteAll(to, fractions, multiplets, blueprint, control)

#process data
cn <- colnames(blueprint)
rn <- colnames(multiplets)

fractions <- CIMseq:::.processSwarm(opt.out, cn, rn, TRUE)
costs <- purrr::map_dbl(opt.out, 2)
names(costs) <- rownames(fractions)

#save
save(fractions, costs, file = file.path(currPath, "data/output.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))


