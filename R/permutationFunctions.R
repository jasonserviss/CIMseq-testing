
#' permuteGroupMeans
#'
#' Permutes the matrix in the groupMeans slot of the spUnsupervised object that
#' is input. To be called repeatedly to generate multiple permutations and used
#' for evaluating the performance of the spSwarm method. Permutation is per row
#' and, thus, the values in each row are randomly shuffled.
#'
#' @name permuteGroupMeans
#' @rdname permuteGroupMeans
#' @aliases permuteGroupMeans
#' @param spUnsupervised An spUnsupervised object.
#' @param ... additional arguments to pass on
#' @return Outputs one matrix of permuted group means with dims equal to the
#' origional group means matrix.
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export

permuteGroupMeans <- function(
  groupMeans
){
  t(apply(groupMeans, 1, sample))
}

#' permuteSwarm
#'
#' Write something
#'
#' Description
#'
#' @name permuteSwarm
#' @rdname permuteSwarm
#' @aliases permuteSwarm
#' @param spCountsMul An spCounts object containing multiplets.
#' @param spUnsupervised An spUnsupervised object.
#' @param distFun The distance function used to calculate the cost. Either the
#'    name of a custom function in the local environment or one of the included
#'    functions, i.e. \code{distToSlice, distToSliceNorm, distToSliceTop,
#'    distToSliceEuclid, distToSlicePearson, bic}.
#' @param maxiter pySwarm argument indicating maximum optimization iterations.
#' @param swarmsize pySwarm argument indicating the number of swarm particals.
#' @param cores The number of cores to be used while running spSwarm.
#' @param norm Logical indicating if the sum of fractions should equal 1.
#' @param iter The number of permutations to perform.
#' @param cellNumbers Tibble; Output from \code{estimateCells} function.
#' @param e Numeric; The epsilon value for the .complexityPenilty unit.
#' @param ... additional arguments to pass on
#' @return ?.
#' @author Jason T. Serviss
#' @keywords permuteSwarm
#' @examples
#'
#' #use demo data
#'
#'
#' #run function
#'
#'
NULL

#' @rdname permuteSwarm
#' @export

setGeneric("permuteSwarm", function(
    spCountsMul,
    ...
){
    standardGeneric("permuteSwarm")
})

#' @rdname permuteSwarm
#' @export
#' @import sp.scRNAseq
#' @importFrom dplyr pull
#' @importFrom parallel mclapply

setMethod("permuteSwarm", "spCounts", function(
  spCountsMul,
  spUnsupervised,
  distFun = sp.scRNAseq:::dtsnCellNum,
  maxiter = 10,
  swarmsize = 50,
  cores = 1,
  norm = TRUE,
  iter,
  cellNumbers,
  e = 0.0025,
  ...
){
  distFun <- match.fun(distFun)
  geneIdx <- getData(spUnsupervised, "selectInd")
  groupMeans <- getData(spUnsupervised, "groupMeans")[geneIdx, ]
  multiplets <- getData(spCountsMul, "counts.cpm")[geneIdx, ]
  
  fractions <- rep(1.0 / ncol(groupMeans), ncol(groupMeans))
  to <- if(ncol(multiplets) == 1) {to <- 1} else {to <- ncol(multiplets)}
  control <- list(maxit = maxiter, s = swarmsize)
  matchIdx <- match(colnames(multiplets), pull(cellNumbers, sampleName))
  cellNumbers <- pull(cellNumbers, cellNumberMedian)[matchIdx]
  
  permData <- .runPermutations(
    iter,
    multiplets,
    groupMeans,
    distFun = distFun,
    fractions,
    to,
    control,
    cores = cores,
    cellNumbers,
    e,
    ...
  )
  
  return(permData)
})

.runPermutations <- function(
  iter,
  multiplets,
  groupMeans,
  distFun = distFun,
  fractions,
  to,
  control,
  cores,
  cellNumbers,
  e,
  ...
){
  permData <- list()
  
  for(y in 1:iter) {
    set.seed(y)
    permData[[y]] <- mclapply(1:to, function(i, ...) {
      sp.scRNAseq:::.optim.fn(
        i,
        fractions,
        distFun,
        permuteGroupMeans(groupMeans),
        control,
        multiplets,
        cellNumbers,
        e,
        ...
      )
    }, mc.cores = cores)
  }
  return(permData)
}

#sObj should be the spSwarm object with the real data
#permData should be the output of the permuteSwarm function

#processPermData <- function(sObj, permData) {
  #setup real data
#  real <- tibble(
#    sample = rownames(getData(sObj, "spSwarm"))
#    cost <- getData(sObj, "cost")
#
#  )
  #get cost and fractions (normalized) and ? from permData
  #cost (called "value" in permData object)
#  costs <- map(permData, function(x) {
#    map_dbl(x, function(y) {
#      y[['value']]
#    })
#  }) %>%
#  map(., as_tibble) %>%
#  bind_rows() #output needs to be named
  
#  fracs <- map(permData, function(x) {
#    map_dfr(x, function(y) {
#      tibble(list(y[['par']]))
#    })
#  }) %>%
#  bind_rows() #output needs to be named
#}
