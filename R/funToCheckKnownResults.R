
#' checkResults
#'
#' Untested for self connections.
#'
#' @name checkResults
#' @rdname checkResults
#' @aliases checkResults
#' @param sObj A spSwarm object.
#' @param known A spSwarm object of the known result. The results should be
#' represented in the spSwarm slot according to the edge.cutoff provided.
#' @param edge.cutoff Standard definition.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq

checkResults <- function(
  sObj,
  known,
  edge.cutoff,
  ...
){
  
  detected <- getEdgesForMultiplet(
    sObj,
    edge.cutoff,
    rownames(getData(sObj, "spSwarm"))
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  expected <- getEdgesForMultiplet(
    known,
    edge.cutoff,
    rownames(getData(sObj, "spSwarm"))
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  full_join(detected, expected, by = "multiplet") %>%
  mutate(
    tp = .tp(.),
    fp = .fp(.),
    fn = .fn(.),
    tn = .tn(.),
    TPR = .TPR(tp, fn),
    TNR = .TNR(tn, fp),
    ACC = .ACC(tp, tn, fp, fn),
  )
}

#calculates the true positives
.tp <- function(data) {
  data %>%
  {map2_int(.$data.x, .$data.y, function(detected, expected) {
    sum(unlist(expected) %in% unlist(detected))
  })}
}

#calculates the false positives
.fp <- function(data) {
  data %>%
  {map2_int(.$data.x, .$data.y, function(detected, expected) {
    sum(!unlist(detected) %in% unlist(expected))
  })}
}

#calculates the false negatives
.fn <- function(data) {
  data %>%
  {map2_int(.$data.x, .$data.y, function(detected, expected) {
    sum(!unlist(expected) %in% unlist(detected))
  })}
}

#calculates the true negatives
#in order to know which comdinations are possible and which are missing from
# both expected and detected, the .possibleCombs function is used.
.tn <- function(data) {
  possibleCombs <- .possibleCombs(data)
  
  data %>%
  add_column(possibleCombs = list(possibleCombs)) %>%
  {pmap_int(
    list(.$data.x, .$data.y, .$possibleCombs),
    function(detected, expected, possibleCombs) {
      bool1 <- !unlist(possibleCombs) %in% unlist(expected)
      bool2 <- !unlist(possibleCombs) %in% unlist(detected)
      sum(bool1 & bool2)
    }
  )}
}

#Calculates all possible connections with all cell types
#Helper for .tn
.possibleCombs <- function(data) {
  data %>%
  select(data.x, data.y) %>%
  unlist() %>%
  str_split("-") %>%
  unlist() %>%
  unique() %>%
  combn(., 2) %>%
  t() %>%
  as_tibble() %>%
  unite(connections, V1, V2, sep = "-")
}

#Calculates the true positive rate (sensitivity)
.TPR <- function(tp, fn) {
  tp / (tp + fn)
}

#Calculates the true negative rate (specificity)
.TNR <- function(tn, fp) {
  tn / (tn + fp)
}

#Calculates the accuracy
.ACC <- function(tp, tn, fp, fn) {
  (tp + tn) / (tp + fn + fp + tn)
}


