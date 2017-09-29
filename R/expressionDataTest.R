
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

.assemble <- function(means, classification, dataset, names, x, type) {

    if(type == "hetero") {
        comb <- combn(unique(classification), x)
    } else {
        comb <- t(matrix(rep(unique(classification), 2), ncol = 2))
    }

    for( i in 1:ncol(comb)) {
        currentMult <- rowMeans(means[ ,comb[ ,i]])
        name <- paste(comb[ ,i], sep = "", collapse = "")

        dataset <- cbind(dataset, currentMult)
        names <- c(names, name)

    }
    return(list(dataset, names))
}

################################################################################
#                                                                              #
# Big dataset                                                                  #
#                                                                              #
################################################################################

#' anyDataTest
#'
#' Subtitle
#'
#'
#' @name anyDataTest
#' @rdname anyDataTest
#' @aliases anyDataTest
#' @param counts matrix. Counts data.
#' @param mulFrac numeric. The factions of the submitted data that should be
#'   reserved for synthesizing multiplets.
#' @param outPath Path to output results to.
#' @param cores Number of cores to be used for spSwarm.
#' @param nMultiplets Number of multuplets.
#' @param nGenes Specifies the \code{n} argument to the \code{spTopMax}
#'    function.
#' @param distFun The distance function to use during swarm optimization.
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords anyDataTest
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

################################################################################
#In order to moderate the fraction of contribution from each cell type in a
# multiplet, require a uObj to be submitted as well as a data.frame indicating
# the fraction of contribution. In this way the user is aware of the
# classifications before running the algo.
################################################################################

anyDataTest <- function(
    counts,
    uObj = NULL,
    mulFrac = 0.8,
    outPath = './data',
    cores = 6,
    nMultiplets = 30,
    nGenes = 2000,
    distFun = "bic",
    target = 20,
    frac = NULL
    ...
){
    
    #make test data
    tmp <- .expressionTestData(
        nMultiplets,
        nGenes,
        counts,
        mulFrac,
        target,
        uObj,
        frac
    )
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

.expressionTestData <- function(
    nMultiplets,
    nGenes,
    counts,
    mulFrac,
    target,
    uObj,
    frac
    ...
){
    
    #################!!!!!!!!!!!!!!!!
    #note that if a uObj is provided the names of the multiplets should correspond to the classification names
    tmp <- .expressionTestMultiplets(
        nGenes,
        counts,
        mulFrac,
        target,
        uObj,
        frac
    )
    singlets <- tmp[[1]]
    multiplets <- tmp[[2]]
    uObj <- tmp[[3]]
    
    #select multiplets for current test
    #(limit test multiplets to half of the total multiplet data... good? bad?)
    if(nMultiplets < (ncol(multiplets) * 0.5)) {
        nM <- nMultiplets
    } else {
        nM <- ncol(multiplets) * 0.5
    }
    idx <- sample(1:ncol(multiplets), nM, replace = FALSE)
    testMultiplets <- multiplets[, idx]
    
    #name, calculate connections, bind, and return
    colnames(multiplets) <- gsub(
        "^([A-Z0-9]*)\\..*",
        "\\1",
        colnames(multiplets)
    )
    colnames(testMultiplets) <- gsub(
        "^([A-Z0-9]*)\\..*",
        "\\1",
        colnames(testMultiplets)
    )
    
    names <- c(
        paste("s", colnames(singlets), 1:ncol(singlets), sep = "."),
        paste("m", colnames(multiplets), 1:ncol(multiplets), sep = ".")
    )
    
    testNames <- c(
    paste("s", colnames(singlets), 1:ncol(singlets), sep = "."),
        paste("m", colnames(testMultuplets), 1:ncol(testMultuplets), sep = ".")
    )
    
    testExpTestData <- cbind(singlets, testMultuplets)
    expTestData <- cbind(singlets, multiplets)
    
    colnames(testExpTestData) <- testNames
    colnames(expTestData) <- names
    
    testExpTestData <- as.matrix(testExpTestData)
    expTestData <- as.matrix(expTestData)
    
    s <- grepl("^m", colnames(testExpTestData))
    table <- quantifyConnections(colnames(testExpTestData)[s])
    
    return(list(expTestData, testExpTestData, uObj, table))
}

.expressionTestMultiplets <- function(
    nGenes,
    counts,
    mulFrac,
    target,
    uObj,
    frac
){
    #subset singlets
    subSngIdx <- sample(
        1:ncol(counts),
        ncol(counts) * (1 - mulFrac),
        replace = FALSE
    )
    subSng <- counts[, subSngIdx]
    
    #make cObjs
    cObjSngTotal <- spCounts(counts, matrix(1, ncol = ncol(counts)))
    cObjSngSub <- spCounts(subSng, matrix(1, ncol = ncol(subSng)))
    
    #run spUnsupervised with total dataset and subset of singlets and check
    # that all cell types are discovered with subset
    if(is.null(uObj)) {
        uObjTotal <- spUnsupervised(
            cObjSngTotal,
            initial_dims = ncol(counts),
            max = nGenes
        )
    } else {
        uObjTotal <- uObj
    }
    
    
    uObjSub <- spUnsupervised(
        cObjSngSub,
        initial_dims = ncol(counts),
        max = nGenes
    )
    
    classTotal <- getData(uObjTotal, "classification")
    classSub <- getData(uObjSub, "classification")
    
    if(!all(unique(classTotal) %in% unique(classSub))) {
        m <- "All classes were not detected with the current mulFrac. Lower it."
        stop(m)
    }
    
    #subset cells for multiplets
    subMul <- counts[, -subSngIdx]
    
    #your going to need to rename the dimensions according to the standards used
    colnames(subSng) <- paste(classSub)
    colnames(subMul) <- paste(classTotal[-subSngIdx])
    
    #make base multiplets
    cellTypes <- unique(classTotal[-subSngIdx])
    tmp <- data.frame(row.names = 1:nGenes)
    doublets <- makeMultuplet(2, cellTypes, tmp, subMul, frac)
    triplets <- makeMultuplet(3, cellTypes, doublets, subMul, frac)
    quadruplets <- makeMultuplet(4, cellTypes, triplets, subMul, frac)

    #adjust multiplets
    multiplets <- adjustMultuplets(
        subMul,
        quadruplets,
        nGenes,
        length(cellTypes),
        target,
        frac
    )
    
    #adjust self
    multiplets <- adjustSelf(multiplets, subMul)
    
    return(list(subSng, multiplets, uObjSub))
}
