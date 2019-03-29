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
load('../MGA.analysis_colonProximal/data/CIMseqData.rda')
cObjSng.p <- cObjSng
cObjMul.p <- cObjMul
rm(cObjSng, cObjMul)
load('../MGA.analysis_colonDistal/data/CIMseqData.rda')
cObjSng.d <- cObjSng
cObjMul.d <- cObjMul
rm(cObjSng, cObjMul)

#rename classes
renameClasses.d <- function(class) {
  n <- case_when(
    class == "0" ~ "Stem",
    class == "1" ~ "TA",
    class == "2" ~ "Goblet.Junb",
    class == "3" ~ "Goblet",
    class == "4" ~ "Goblet.Plet1",
    class == "5" ~ "Lgr5+.Mki67",
    class == "6" ~ "Colonocytes",
    class == "7" ~ "Goblet.Mki67",
    class == "8" ~ "Tufft",
    class == "9" ~ "Enteroendocrine",
    TRUE ~ "error"
  )
  paste0(n, ".distal")
}
getData(cObjSng.d, "classification") <- renameClasses.d(getData(cObjSng.d, "classification"))

renameClasses.p <- function(class) {
  n <- case_when(
    class == "0" ~ "Stem",
    class == "1" ~ "TA",
    class == "2" ~ "Colonocytes",
    class == "3" ~ "Goblet",
    class == "4" ~ "Lgr5+.Mki67",
    class == "5" ~ "Enteroendocrine",
    class == "6" ~ "Tufft",
    TRUE ~ "error"
  )
  paste0(n, ".proximal")
}
getData(cObjSng.p, "classification") <- renameClasses.p(getData(cObjSng.p, "classification"))

#get features
stopifnot(identical(rownames(getData(cObjSng.p, "counts")), rownames(getData(cObjSng.d, "counts"))))
features <- unique(c(getData(cObjMul.p, "features"), getData(cObjMul.d, "features")))
shared.features <- intersect(
  rownames(MGAB.Counts), 
  rownames(getData(cObjSng.p, "counts"))[features]
)

blueprint.cpm <- cbind(getData(cObjSng.p, "counts.cpm"), getData(cObjSng.d, "counts.cpm"))
blueprint <- averageGroupExpression(
  blueprint.cpm[which(rownames(blueprint.cpm) %in% shared.features), ], 
  c(getData(cObjSng.p, "classification"), getData(cObjSng.d, "classification"))
)
blueprint <- blueprint[, order(colnames(blueprint))]

#setup fractions variable
nuc <- length(unique(c(getData(cObjSng.p, "classification"), getData(cObjSng.d, "classification"))))
fractions <- rep(1.0 / nuc, nuc)
stopifnot(identical(length(fractions), ncol(blueprint)))

#setup multiplets
MGAB.Cpm <- cpm(MGAB.Counts[, filter(MGAB.Meta, sub_tissue == "colon")$sample])

features.bulk <- which(rownames(MGAB.Cpm) %in% shared.features)
multiplets <- matrix(
  MGAB.Cpm[features.bulk, ],
  ncol = ncol(MGAB.Cpm),
  dimnames = list(rownames(MGAB.Cpm)[features.bulk], colnames(MGAB.Cpm))
)

#setup control
to <- dim(multiplets)[2]
control <- list(maxit = 100, s = 1000)

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
