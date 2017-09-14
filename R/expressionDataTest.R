
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
#' @param ngenes Specifies the \code{n} argument to the \code{spTopMax}
#'    function.
#' @param distFun The distance function to use during swarm optimization.
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

expressionDataTest <- function(
    outPath = './data',
    cores = 6,
    n = 30,
    ngenes = 2000,
    distFun = bic,
    ...
){
    
    #make test data
    tmp <- .expressionTestData(n, ngenes)
    expressionTestData <- tmp[[1]]
    testExpTestData <- tmp[[2]]
    expressionTestUnsupervised <- tmp[[3]]
    expressionTestTable <- tmp[[4]]
    
    #make spCounts object
    single <- grepl("^s", colnames(testExpTestData))
    sng <- testExpTestData[ ,single]
    mul <- testExpTestData[ ,!single]
    cObjSng <- spCounts(sng, matrix(NA, ncol=ncol(sng)))
    cObjMul <- spCounts(mul, matrix(NA, ncol=ncol(mul)))
    
    #run pySwarm
    expressionTestSwarm <- spSwarm(
        cObjMul,
        expressionTestUnsupervised,
        maxiter = 10^4,
        swarmsize = 250,
        cores = cores,
        distFun = distFun
    )
    
    #save
    save(
        expressionTestData,
        file = paste(outPath, "expressionTestData.rda", sep = "/"),
        compress = "bzip2"
    )
    save(
        testExpTestData,
        expressionTestUnsupervised,
        expressionTestTable,
        expressionTestSwarm,
        file = paste(outPath, "expressionTest.rda", sep = "/"),
        compress = "bzip2"
    )
}

.expressionTestData <- function(n, ngenes) {
    
    #load data
    counts <- expCounts
    single <- grepl("^s", colnames(counts))
    sng <- counts[ ,single]
    mul <- counts[ ,!single]
    cObjSng <- spCounts(sng, matrix(NA, ncol = ncol(sng)))
    cObjMul <- spCounts(mul, matrix(NA, ncol = ncol(mul)))
    
    uObj <- spUnsupervised(
        cObjSng,
        theta = 0,
        k = 2,
        max_iter = 10^4,
        perplexity = 10,
        initial_dims = ncol(sng),
        Gmax = 50,
        seed = 15,
        type = "max",
        max = ngenes
    )
    
    classification <- getData(uObj, "classification")
    means <- getData(uObj, "groupMeans")
    
    dataset <- data.frame(row.names = 1:nrow(sng))
    names <- c()
    
    #same cell type
    tmp <- .assemble(means, classification, dataset, names, 2, type = "homo")
    dataset <- tmp[[1]]
    names <- tmp[[2]]

    #doublets
    tmp <- .assemble(means, classification, dataset, names, 2, type = "hetero")
    dataset <- tmp[[1]]
    names <- tmp[[2]]
    
    #triplets
    tmp <- .assemble(means, classification, dataset, names, 3, type = "hetero")
    dataset <- tmp[[1]]
    names <- tmp[[2]]
    
    #quadruplets
    tmp <- .assemble(means, classification, dataset, names, 4, type = "hetero")
    dataset <- tmp[[1]]
    names <- tmp[[2]]
    
    #adjust colnames
    colnames(dataset) <- paste("m.", names, sep = "")
    colnames(sng) <- paste("s.", classification, sep = "")
    
    #select multuplets for current test
    testMultuplets <- dataset[ ,sample(1:ncol(dataset), n, replace = FALSE)]
    table <- calculateConnections(testMultuplets, type = "multuplets")
    
    expTestData <- as.matrix(cbind(sng, dataset))
    testExpTestData <- as.matrix(cbind(sng, testMultuplets))

    #add new data to uObj
    #counts(unsupervised) <- testExpTestData
    #counts.log(unsupervised) <- sp.scRNAseq:::.norm.log.counts(testExpTestData)
    #sampleType(unsupervised) <- c(rep("Singlet", ncol(sng)), rep("Multuplet", n))
    
    return(list(expTestData, testExpTestData, uObj, table))
    
}

#.assemble <- function(means, classification, dataset, names, x, type) {
#
#    if(type == "hetero") {
#        comb <- combn(unique(classification), x)
#    } else {
#        comb <- t(matrix(rep(unique(classification), 2), ncol = 2))
#    }
#
#    for( i in 1:ncol(comb)) {
#        currentMult <- rowMeans(means[ ,comb[ ,i]])
#        name <- paste(comb[ ,i], sep = "", collapse = "")
#
#        dataset <- cbind(dataset, currentMult)
#        names <- c(names, name)
#
#    }
#    return(list(dataset, names))
#}

################################################################################
#                                                                              #
# Big dataset                                                                  #
#                                                                              #
################################################################################

#' expressionBigTest
#'
#' Subtitle
#'
#'
#' @name expressionBigTest
#' @rdname expressionBigTest
#' @aliases expressionBigTest
#' @param outPath Path to output results to.
#' @param cores Number of cores to be used for spSwarm.
#' @param n Number of multuplets.
#' @param ngenes Specifies the \code{n} argument to the \code{spTopMax}
#'    function.
#' @param distFun The distance function to use during swarm optimization.
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords expressionBigTest
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

expressionBigTest <- function(
    outPath = './data',
    cores = 6,
    n = 30,
    ngenes = 2000,
    distFun = bic,
    ...
){
    
    #make test data
    tmp <- .expressionTestData(n, ngenes)
    expressionTestData <- tmp[[1]]
    testExpTestData <- tmp[[2]]
    expressionTestUnsupervised <- tmp[[3]]
    expressionTestTable <- tmp[[4]]
    
    #make spCounts object
    single <- grepl("^s", colnames(testExpTestData))
    sng <- testExpTestData[ ,single]
    mul <- testExpTestData[ ,!single]
    cObjSng <- spCounts(sng, matrix(NA, ncol=ncol(sng)))
    cObjMul <- spCounts(mul, matrix(NA, ncol=ncol(mul)))
    
    #run pySwarm
    expressionTestSwarm <- spSwarm(
    cObjMul,
    expressionTestUnsupervised,
    maxiter = 10^4,
    swarmsize = 250,
    cores = cores,
    distFun = distFun
    )
    
    #save
    save(
    expressionTestData,
    file = paste(outPath, "expressionTestData.rda", sep = "/"),
    compress = "bzip2"
    )
    save(
    testExpTestData,
    expressionTestUnsupervised,
    expressionTestTable,
    expressionTestSwarm,
    file = paste(outPath, "expressionTest.rda", sep = "/"),
    compress = "bzip2"
    )
}

.expressionTestData <- function(n, ngenes) {
    
    #load data
    counts <- bigPancreasDataset
    classification <- pull(bigPancreasClasses, class)
    
    #calculate group means
    c <- unique(classification)
    means <- lapply(c, function(x) {
        rowMeans(counts[, classification == x])
    })
    means <- as.matrix(as.data.frame(means))
    colnames(means) <- c
    
    dataset <- data.frame(row.names = 1:nrow(counts))
    names <- c()
    
    #same cell type
    tmp <- .assemble(counts, classification, dataset, names, 2, type = "homo")
    dataset <- tmp[[1]]
    names <- tmp[[2]]
    
    #doublets
    tmp <- .assemble(counts, classification, dataset, names, 2, type = "hetero")
    dataset <- tmp[[1]]
    names <- tmp[[2]]
    
    #triplets
    tmp <- .assemble(counts, classification, dataset, names, 3, type = "hetero")
    dataset <- tmp[[1]]
    names <- tmp[[2]]
    
    #quadruplets
    tmp <- .assemble(counts, classification, dataset, names, 4, type = "hetero")
    dataset <- tmp[[1]]
    names <- tmp[[2]]
    
    #adjust colnames
    colnames(dataset) <- paste("m.", names, sep = "")
    colnames(sng) <- paste("s.", classification, sep = "")
    
    #select multuplets for current test
    testMultuplets <- dataset[ ,sample(1:ncol(dataset), n, replace = FALSE)]
    table <- calculateConnections(testMultuplets, type = "multuplets")
    
    expTestData <- as.matrix(cbind(sng, dataset))
    testExpTestData <- as.matrix(cbind(sng, testMultuplets))
    
    #add new data to uObj
    #counts(unsupervised) <- testExpTestData
    #counts.log(unsupervised) <- sp.scRNAseq:::.norm.log.counts(testExpTestData)
    #sampleType(unsupervised) <- c(rep("Singlet", ncol(sng)), rep("Multuplet", n))
    
    return(list(expTestData, testExpTestData, uObj, table))
    
}

.assemble <- function(
    counts,
    classification,
    dataset,
    names,
    x,
    type
){
    
    if(type == "hetero") {
        comb <- combn(unique(classification), x)
    } else {
        comb <- t(matrix(rep(unique(classification), 2), ncol = 2))
    }
    
    #find indexes
    idxes <- lapply(1:ncol(comb), function(x) {
        lapply(1:length(comb[, x]), function(y) {
            which(classification %in% comb[y, x])
        })
    })
    
    #sample indexes
    rand <- lapply(idxes, function(x) {
        sapply(x, function(y) {
            y[sample(1:length(y), 1)]
        })
    })
    
    #calculate means
    mul <- sapply(rand, function(x) {
        rowMeans(counts[, x])
    })
    
    #track names
    names <- sapply(1:ncol(comb), function(x)
        paste(comb[, x], collapse = "")
    )
    
    #bind and return
    dataset <- cbind(dataset, mul)
    return(list(dataset, names))
}





