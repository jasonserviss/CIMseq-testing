
#' checkResults
#'
#' Untested for self connections.
#'
#' @name checkResults
#' @rdname checkResults
#' @aliases checkResults
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param known A CIMseqSwarm object of the known result. The results should be
#' represented in the fractions slot according to the edge.cutoff provided.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param edge.cutoff Standard definition.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export
#' @import CIMseq
#' @importFrom tidyr unite nest
#' @importFrom dplyr full_join mutate
#' @importFrom purrr map2_int pmap_int
#' @importFrom tibble as_tibble
#' @importFrom magrittr "%>%"

checkResults <- function(
  swarm, known, singlets, edge.cutoff, ...
){

  detected <- CIMseq::getEdgesForMultiplet(
    swarm, singlets, edge.cutoff, rownames(getData(swarm, "fractions"))
  ) %>%
    unite(connections, from, to, sep = "-") %>%
    distinct() %>%
    nest(-multiplet)

  expected <- CIMseq::getEdgesForMultiplet(
    known, singlets, edge.cutoff, rownames(getData(swarm, "fractions"))
  ) %>%
    unite(connections, from, to, sep = "-") %>%
    distinct() %>%
    nest(-multiplet)

  full_join(detected, expected, by = "multiplet") %>%
    mutate(
      tp = .tp(.),
      fp = .fp(.),
      fn = .fn(.),
      tn = .tn(., sObj, known),
      TPR = .TPR(tp, fn),
      TNR = .TNR(tn, fp),
      ACC = .ACC(tp, tn, fp, fn)
    ) %>%
    rename(data.detected = data.x, data.expected = data.y)
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
.tn <- function(data, sObj, known) {
  possibleCombs <- .possibleCombs(sObj, known)

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
.possibleCombs <- function(sObj, known) {
  ctKnown <- colnames(getData(known, "fractions"))
  ctDetected <- colnames(getData(sObj, "fractions"))

  c(ctKnown, ctDetected) %>%
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

#' printResults
#'
#' Takes the output from checkResults and prints it in a more verbose form.
#' (note: rewrite as general to take an arg for nested columns to expand)
#'
#' @name printResults
#' @rdname printResults
#' @aliases printResults
#' @param data Output from checkResults function.
#' @param addFractions Logical indicating if fractions should be added to the
#'  output.
#' @param spSwarm An spSwarm object. Required if \code{addFractions} is TRUE.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq
#' @importFrom dplyr pull
#' @importFrom stringr str_split
#' @importFrom purrr map_dfr
#' @importFrom tibble add_column column_to_rownames
#' @importFrom magrittr "%>%"

printResults <- function(
  data,
  addFractions = TRUE,
  spSwarm = NULL,
  ...
){
  #expand connections
  data.detected <- select(data, multiplet, data.detected) %>%
  .expandNested()

  data.expected <- select(data, multiplet, data.expected) %>%
  .expandNested()

  expanded <- select(data, -data.detected, -data.expected) %>%
  full_join(data.detected, by = "multiplet") %>%
  full_join(data.expected, by = "multiplet") %>%
  rename(data.detected = connections.x, data.expected = connections.y) %>%
  select(costFunction:cellsInWell, data.expected, data.detected, tp:ACC) %>%
  arrange(costFunction, multiplet)

  #add fractions
  if(addFractions) {
    return(.addFractions(expanded, spSwarm))
  } else {
    return(expanded)
  }
}

.expandNested <- function(contracted) {
  contracted %>%
  unnest() %>%
  bind_rows() %>%
  group_by(multiplet) %>% #change in generalized solution
  summarize(connections = paste(connections, collapse = ", "))
}

.addFractions <- function(expanded, spSwarm) {
  fracs <- getData(spSwarm, "spSwarm")
  names <- paste(colnames(fracs), collapse = ", ")
  formated <- paste("frac (", names, ")", sep = "")

  fracs %>%
  round(digits = 2) %>%
  rownames_to_column("multiplet") %>%
  as_tibble() %>%
  full_join(expanded, by = "multiplet") %>%
  unite(fractions, 2:(ncol(fracs) + 1), sep = ", ") %>%
  setNames(c(colnames(.)[1], formated, colnames(.)[3:ncol(.)])) %>%
  select(multiplet, cellsInWell, 2, data.expected:ACC) %>%
  arrange(multiplet)
}


#' resultsInPlate
#'
#' Takes the output from checkResults and setupPlate and displays results in
#' plate format.
#'
#' @name resultsInPlate
#' @rdname resultsInPlate
#' @aliases resultsInPlate
#' @param results Output from checkResults function.
#' @param plate Output from setup plate or a tibble including the columns \code{
#'  row, column, and multipletName} indicating the position of each multiplet in
#'  the plate.
#' @param var The variable to visualize. May be "tp, fp, tn, fn, TPR, TNR, ACC".
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss

NULL

#' @export
#' @import sp.scRNAseq
#' @importFrom dplyr select full_join
#' @import ggplot2
#' @importFrom ggthemes theme_few
#' @importFrom magrittr "%>%"

resultsInPlate <- function(
  results,
  plate,
  var,
  ...
){

  results %>%
  select(multiplet, var) %>%
  full_join(plate, by = c("multiplet" = "multipletName")) %>%
  ggplot(
    data = .,
    aes(x = column, y = ordered(row, levels = rev(LETTERS[1:8])))
  ) +
  geom_tile(aes_string(fill = var)) +
  scale_x_discrete(position = "top") +
  theme_few() +
  theme(
    axis.title = element_blank(),
    legend.position = "top"
  ) +
  guides(fill = guide_colourbar(barwidth = 13))
}

#' setupPlate
#'
#' Sets up the CIMseqSwarm object for the \code{known} argument to the
#' \code{checkResults} function.
#'
#' @name setupPlate
#' @rdname setupPlate
#' @aliases setupPlate
#' @param plateData tibble; Includes the columns \code{row, column,
#'  multipletName, multipletComposition, connections} where \code{row} and
#'  \code{column} includes a letter or numerical value, respectivley, indicating
#'  the index in the plate of each multiplet, \code{multipletName} corresponds
#'  to the name of the multiplet and matches the names in the subsequent
#'  detected results, \code{multipletComposition} is the cell types included in
#'  the multiplet seperated by a "-", and \code{connections} is a list for each
#'  multiplet containing a matrix with all possible combinations of 2 of the
#'  cell types included in the multiplet.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss
#'
#'
#'
NULL

#' @export
#' @import CIMseq
#' @importFrom dplyr pull
#' @importFrom stringr str_split
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble add_column column_to_rownames
#' @importFrom magrittr "%>%"

setupPlate <- function(
  plateData,
  ...
){

  #make spSwarm slot for spSwarm object
  spSwarm <- plateData %>%
    filter(cellNumber == "Multiplet") %>%
    mutate(connections = str_split(cellTypes, "-")) %>%
    mutate(connections = {map(.$connections, combn, 2)}) %>%
    {map_dfr(.$connections, function(x) {
      cellTypes <- .getCellTypes(plateData)
      vec <- vector(mode = "numeric", length = length(cellTypes))
      names(vec) <- cellTypes
      vec[names(vec) %in% unique(as.character(x))] <- 1 / length(unique(as.character(x)))
      as.data.frame(t(as.data.frame(vec)))
    })} %>%
    add_column(rowname = filter(plateData, cellNumber == "Multiplet")$sample) %>%
    column_to_rownames()

  #create spSwarm object
  new("CIMseqSwarm",
    fractions = as.matrix(spSwarm),
    costs = vector(mode = "numeric"),
    convergence = vector(mode = "character"),
    stats = tibble(),
    singletIdx = list(),
    arguments = tibble()
  )
}

.getCellTypes <- function(plateData) {
  plateData %>%
    pull(cellTypes) %>%
    str_split("-") %>%
    unlist() %>%
    unique() %>%
    sort()
}

#' viewAsPlate
#'
#' Takes the output from setupPlate and reformats it to a 96-well plate format.
#'
#' @name viewAsPlate
#' @rdname viewAsPlate
#' @aliases viewAsPlate
#' @param plate Output from setup plate.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss

NULL

#' @export
#' @importFrom dplyr select
#' @importFrom magrittr "%>%"
#' @importFrom tidyr spread

viewAsPlate <- function(plate) {
  plate %>%
  select(row, column, multipletComposition) %>%
  spread(column, multipletComposition)
}
