#' permuteMeans
#'
#' Sets up permutations of the classes for evaluating the performance
#' of various distance functions used with the spSwarm method. Subsequently
#' used with the \code{permutationMeans} function.
#'
#' @name permuteMeans
#' @rdname permuteMeans
#' @aliases permuteMeans
#' @param spCountsMul A spCounts object containing multiplets.
#' @param spUnsupervised An spUnsupervised object.
#' @param nPerms Number of total permutations. If NULL set to 100 per multiplet.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq
#' @importFrom purrr map

permuteMeans <- function(
  spCountsMul,
  spUnsupervised,
  nPerms = NULL
){
  means <- getData(spUnsupervised, "groupMeans")
  idx <- getData(spUnsupervised, "selectInd")
  l <- ncol(means)
  if(is.null(nPerms)) {nPerms <- 100 * ncol(getData(spCountsMul, "counts"))}
  
  pIdx <- map(1:(nPerms * length(idx)), function(x) {
    sample(seq(1, l, 1), l, replace = FALSE)
  }) %>%
  unlist()
  
  sIdx <- matrix(
    c(rep(sort(rep(1:length(idx), l)), nPerms), pIdx),
    ncol = 2
  )
  
  matrix(
    means[sIdx],
    ncol = l,
    dimnames = list(sort(rep(1:nPerms, length(idx))), colnames(means)),
    byrow = TRUE
  )
}

#' permutationCosts
#'
#' Calculates the costs for each permutation using the provided function and
#' returns these together with the real costs. Note that the real costs in the
#' spSwarm object should be calculated with the SAME distance function as
#' provided in the \code{distFun} argument.
#'
#' @name permutationCosts
#' @rdname permutationCosts
#' @aliases permutationCosts
#' @param permutationMeans Output from permutationMeans function.
#' @param spCountsMul An spCounts object containing multiplets.
#' @param spUnsupervised An spUnsupervised object.
#' @param nPermsPerMul Numeric; Number of permutations per multiplet. If null
#'  this is set to 100.
#' @param Function; A distance function that returns the cost.
#' @param ... additional arguments to pass on.
#' @return matrix
#' @author Jason T. Serviss
#'
#'
#'
NULL

#' @export
#' @import sp.scRNAseq
#' @importFrom purrr map_dfc
#' @importFrom dplyr mutate
#' @importFrom purrr map2_int pmap_int
#' @importFrom tibble add_column
#' @importFrom tidyr gather spread
#' @importFrom magrittr "%>%"

permutationCosts <- function(
  permutationMeans,
  spCountsMul,
  spUnsupervised,
  spSwarm,
  nPermsPerMul = NULL,
  distFun = sp.scRNAseq:::distToSlice,
  ...
){
  idx <- getData(spUnsupervised, "selectInd")
  multiplets <- getData(spCountsMul, "counts.cpm")[idx, ]
  if(is.null(nPermsPerMul)) {nPermsPerMul <- 100}
  
  l <- getData(spUnsupervised, "classification") %>%
  unique() %>%
  length()
  
  real <- tibble(
    key = rownames(getData(spSwarm, "spSwarm")),
    cost = getData(spSwarm, "costs")
  )
  
  map_dfc(
    1:ncol(multiplets),
    ~.cost(., multiplets[, .], nPermsPerMul, permutationMeans, l, distFun)
  ) %>%
  setNames(colnames(multiplets)) %>%
  add_column(tmp = 1:nPermsPerMul) %>%
  gather(key, value, -tmp) %>%
  spread(tmp, value) %>%
  left_join(real, by = "key") %>%
  gather(iteration, cost, -key) %>%
  mutate(Type = if_else(iteration == "cost", "real", "permuted")) %>%
  select(-iteration)
}

.cost <- function(i, multiplet, nperms, means, l, distFun) {
  pIdx <- ((i * nperms) - (nperms - 1)):(i * nperms)
  
  map_dbl(1:nperms, function(x) {
    distFun(rep(1/l, l), means[rownames(means) == pIdx[x], ], multiplet)
  }) %>%
  as_tibble()
}

#' plotDistPermutations
#'
#' Plots the results from the \code{permutationCosts} function.
#'
#' @name plotDistPermutations
#' @rdname plotDistPermutations
#' @aliases plotDistPermutations
#' @param permCosts An spCounts object containing multiplets.
#' @param facet An spUnsupervised object.
#' @param ... additional arguments to pass on.
#' @return matrix
#' @author Jason T. Serviss
#'
#'
#'
NULL

#' @export
#' @import sp.scRNAseq
#' @import ggplot2
#' @importFrom ggthemes theme_few
#' @importFrom magrittr "%>%"

plotDistPermutations <- function(
  permCosts,
  facet = FALSE,
  ...
){
  p <- permCosts %>%
  ggplot() +
  geom_histogram(aes(x = cost, fill = Type), binwidth = 10000) +
  theme_few() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "top"
  )
  
  if(facet) {
    p <- p + facet_wrap(~key)
    p
    return(p)
  } else {
    p
    return(p)
  }
}
