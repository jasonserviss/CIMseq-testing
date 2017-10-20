#' countsSortedTest
#'
#' Subtitle
#'
#'
#' @name countsSortedTest
#' @rdname countsSortedTest
#' @aliases countsSortedTest
#' @param data .
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @examples
#'
#' #use demo data
#' \dontrun{countsSortedTest(data = "countsSorted1.rda")}
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq

countsSortedTest <- function(
    data = countsSorted1,
    ...
){
    
    #make spCounts objects
    data(data) #fix this
    s <- grepl("s.", colnames(countsSorted1))
    cObjSng <- spCounts(countsSorted1[, s], countsSortedERCC1[, s])
    cObjMul <- spCounts(countsSorted1[, !s], countsSortedERCC1[, !s])
    
    #make spUnsupervised object
    select <- spTopMax(cObjSng, 2000)
    my.dist <- pearsonsDist(cObjSng, select)
    my.tsne <- runTsne(my.dist, max_iter = 10000, perplexity = 30)
    my.mclust <- runMclust(my.tsne)
    
}