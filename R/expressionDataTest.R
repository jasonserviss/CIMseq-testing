
#' expressionDataTest
#'
#' Subtitle
#'
#'
#' @name expressionDataTest
#' @rdname expressionDataTest
#' @aliases expressionDataTest
#' @param outPath Path to output results to.
#' @param cores Number of cores to be used for spSwarm.
#' @param n Number of multuplets.
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords expressionDataTest
#' @examples
#'
#' #use demo data
#' expressionDataTest()
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq

expressionDataTest <- function(outPath = 'data', cores=5, n=15, ...) {
    
    #make test data
    tmp <- .expressionTestData(n)
    expressionTestData <- tmp[[1]]
    testExpTestData <- tmp[[2]]
    expressionTestUnsupervised <- tmp[[3]]
    expressionTestTable <- tmp[[4]]
    
    #make spCounts object
    expressionTestDataCounts <- spCounts(testExpTestData, matrix(), "m.")
    
    #run pySwarm
    expressionTestSwarm <- spSwarm(expressionTestUnsupervised, limit="none", maxiter=10, swarmsize=250, cores=cores)
    
    #save
    save(
        expressionTestData,
        file=paste(outPath, "expressionTestData.rda", sep="/"),
        compress="bzip2"
    )
    save(
        expressionTestDataCounts,
        file=paste(outPath, "expressionTestDataCounts.rda", sep="/"),
        compress="bzip2"
    )
    save(
        expressionTestUnsupervised,
        file=paste(outPath, "expressionTestUnsupervised.rda", sep="/"),
        compress="bzip2"
    )
    save(
        expressionTestSwarm,
        file=paste(outPath, "expressionTestSwarm.rda", sep="/"),
        compress="bzip2"
    )
    save(
        expressionTestTable,
        file=paste(outPath, "expressionTestTable.rda", sep="/"),
        compress="bzip2"
    )
}

.expressionTestData <- function(n) {
    
    #load data
    counts <- getData(expData, "counts")
    sampleType <- getData(expData, "sampleType")
    sng <- counts[ ,sampleType == "Singlet"]
    
    unsupervised <- spUnsupervised(expData, max=2500, max_iter=10^5)
    
    classification <- getData(unsupervised, "classification")
    means <- getData(unsupervised, "groupMeans")
    
    dataset <- data.frame(row.names=1:nrow(sng))
    names <- c()
    
    #same cell type
    tmp <- .assemble(means, classification, dataset, names, 2, type="homo")
    dataset <- tmp[[1]]
    names <- tmp[[2]]

    #doublets
    tmp <- .assemble(means, classification, dataset, names, 2, type="hetero")
    dataset <- tmp[[1]]
    names <- tmp[[2]]
    
    #triplets
    tmp <- .assemble(means, classification, dataset, names, 3, type="hetero")
    dataset <- tmp[[1]]
    names <- tmp[[2]]
    
    #quadruplets
    tmp <- .assemble(means, classification, dataset, names, 4, type="hetero")
    dataset <- tmp[[1]]
    names <- tmp[[2]]
    
    #adjust colnames
    colnames(dataset) <- paste("m.", names, sep="")
    colnames(sng) <- paste("s.", classification, sep="")
    
    #select multuplets for current test
    testMultuplets <- dataset[ ,sample(1:ncol(dataset), n, replace=FALSE)]
    table <- calculateConnections(testMultuplets, type="multuplets")
    
    expTestData <- as.matrix(cbind(sng, dataset))
    testExpTestData <- as.matrix(cbind(sng, testMultuplets))

    #add new data to uObj
    counts(unsupervised) <- testExpTestData
    counts.log(unsupervised) <- sp.scRNAseq:::.norm.log.counts(testExpTestData)
    sampleType(unsupervised) <- c(rep("Singlet", ncol(sng)), rep("Multuplet", n))
    
    return(list(expTestData, testExpTestData, unsupervised, table))
    
}

.assemble <- function(means, classification, dataset, names, x, type) {
    
    if(type=="hetero") {
        comb <- combn(unique(classification), x)
    } else {
        comb <- t(matrix(rep(unique(classification), 2), ncol=2))
    }
    
    for( i in 1:ncol(comb)) {
        currentMult <- rowMeans(means[,comb[,i]])
        name <- paste(comb[,i], sep="", collapse="")
        
        dataset <- cbind(dataset, currentMult)
        names <- c(names, name)
        
    }
    return(list(dataset, names))
}