
#' Generate a summary when testing different distance functions.
#'
#' Generates a full summary from a list of spSwarm objects that compares results
#' to the expected results and generates TRP, TNR, and accuracy metrics.
#'
#' @name generateTestSummary
#' @rdname generateTestSummary
#' @aliases generateTestSummary
#' @param spSwarms Output from testDistanceFunctions function.
#' @param metadata A tibble including sample metadata. Provided to
#'  \code{setupPlate} function.
#' @param ... additional arguments to pass on
#' @return mean accuracy
#' @author Jason T. Serviss
#'
#'
#'
NULL

#' @export
#' @importFrom dplyr "%>%" mutate case_when select
#' @importFrom purrr map_dfr
#' @importFrom stringr str_extract

generateTestSummary <- function(spSwarms, metadata, ...) {
  known <- setupPlate(metadata)
  spSwarms %>%
    map_dfr(
      ., ~checkResults(.x, known, edge.cutoff = 0),
      .id = "costFunction"
    ) %>%
    mutate(cellsInWell = case_when(
      str_extract(multiplet, "..$") %in% c("01", "02", "03", "04") ~ 2L,
      str_extract(multiplet, "..$") %in% c("05", "06", "07", "08") ~ 3L,
      str_extract(multiplet, "..$") %in% c("09", "10", "11", "12") ~ 4L
    )) %>%
    select(costFunction, multiplet, cellsInWell, data.detected:ACC)
}

#' Tests multiple distance functions for comparison.
#'
#' A wrapper for testing improvments in distance/cost functions.
#'
#' @name testDistanceFunctions
#' @rdname testDistanceFunctions
#' @aliases testDistanceFunctions
#' @param spCountsMul An spCounts object containing multiplets.
#' @param spUnsupervised An spUnsupervised object.
#' @param distFuns A named list of distance functions to test.
#' @param params A list with the same number of elements as distFuns containing
#'  additional parameters that should be passed to each distFun.
#' @param ... additional arguments to pass on
#' @return The summary.
#' @author Jason T. Serviss
NULL

#' @export
#' @import CIMseq
#' @importFrom purrr map map2 invoke_map
#' @importFrom dplyr "%>%"

testDistanceFunctions <- function(
  spCountsMul,
  spUnsupervised,
  distFuns,
  params,
  ...
){
  params <- map(params, function(x) {
    c(list(spCounts = spCountsMul, spUnsupervised = spUnsupervised), x)
  })
  
  map2(distFuns, params, ~c(distFun =.x, .y)) %>%
  invoke_map(spSwarm, .) %>%
  #map(~do.call("spSwarm", args = .x)) %>%
    setNames(names(distFuns))
}
