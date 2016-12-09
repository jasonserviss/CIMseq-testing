
#' experimentalDataTest
#'
#' Subtitle
#'
#'
#' @name experimentalDataTest
#' @rdname experimentalDataTest
#' @aliases experimentalDataTest
#' @param outPath Path to output results to.
#' @param cores Number of cores to be used for spSwarm.
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords experimentalDataTest
#' @examples
#'
#' #use demo data
#' expDataTest()
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq

experimentalDataTest <- function(outPath = 'data', cores=5, ...) {
    #load test data
    data(expData)
    
    #make spCounts object
    cObj <- expData
    
    #run unsupervised learning
    experimentalDataUnsupervised <- spUnsupervised(
        cObj,
        max=2000,
        max_iter=2*10^4
    )
    
    #run pySwarm
    experimentalDataSwarm <- spSwarm(
        experimentalDataUnsupervised,
        limit="none",
        maxiter=10,
        swarmsize=250,
        cores=cores
    )
    
    #save
    save(
        experimentalDataUnsupervised,
        file=paste(outPath, "experimentalDataUnsupervised.rda", sep="/"),
        compress="bzip2"
    )
    save(
        experimentalDataSwarm,
        file=paste(outPath, "experimentalDataSwarm.rda", sep="/"),
        compress="bzip2"
    )
}