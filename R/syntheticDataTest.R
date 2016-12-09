
#' syntheticDataTest
#'
#' Subtitle
#'
#'
#' @name syntheticDataTest
#' @rdname syntheticDataTest
#' @aliases syntheticDataTest
#' @param outPath Path to output results to.
#' @param cores Number of cores to be used for spSwarm.
#' @param n Number of multuplets.
#' @param ngenes Number of genes in generated synthetic data.
#' @param ncells Number of cells per cell type in generated synthetic data.
#' @param cellTypes Number of cell types in generated synthetic data.
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords syntheticDataTest
#' @examples
#'
#' #use demo data
#' syntheticDataTest()
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq

syntheticDataTest <- function(outPath = 'data', cores = 5, n=15, ngenes=2000, ncells=100, cellTypes=10, ...) {
    
    #create synthetic data
    tmp <- .syntheticTestData(n, ngenes, ncells, cellTypes)
    syntheticData <- tmp[[1]]
    syntheticDataTest <- tmp[[2]]
    syntheticDataUnsupervised <- tmp[[3]]
    syntheticDataTable <- tmp[[4]]
    
    #make spCounts object
    syntheticDataCounts <- spCounts(syntheticDataTest, matrix(), sampleType = "m.")
    
    #run pySwarm
    syntheticDataSwarm <- spSwarm(
        syntheticDataUnsupervised,
        limit="none",
        maxiter=10,
        swarmsize=250,
        cores=cores
    )
    
    #save
    save(
        syntheticData,
        file=paste(outPath, "syntheticData.rda", sep="/"),
        compress="bzip2"
    )
    save(
        syntheticDataCounts,
        file=paste(outPath, "syntheticDataCounts.rda", sep="/"),
        compress="bzip2"
    )
    save(
        syntheticDataUnsupervised,
        file=paste(outPath, "syntheticDataUnsupervised.rda", sep="/"),
        compress="bzip2"
    )
    save(
        syntheticDataSwarm,
        file=paste(outPath, "syntheticDataSwarm.rda", sep="/"),
        compress="bzip2"
    )
    save(
        syntheticDataTable,
        file=paste(outPath, "syntheticDataTable.rda", sep="/"),
        compress="bzip2"
    )
}

.syntheticTestData <- function(n, ngenes, ncells, cellTypes, perplexity=10) {
    tmp <- .syntheticMultuplets(ngenes, ncells, cellTypes, perplexity)
    singlets <- tmp[[1]]
    multuplets <- tmp[[2]]
    uObj <- tmp[[3]]
    
    #select multuplets for current test
    testMultuplets <- multuplets[ ,sample(1:ncol(multuplets), n, replace=FALSE)]
    table <- calculateConnections(testMultuplets, type="multuplets")
    
    #name, bind, and return
    names <- c(
        paste("s", colnames(singlets), sep="."),
        paste("m", colnames(multuplets), sep=".")
    )
    
    testNames <- c(
        paste("s", colnames(singlets), sep="."),
        paste("m", colnames(testMultuplets), sep=".")
    )
    
    testSyntheticData <- cbind(singlets, testMultuplets)
    syntheticData <- cbind(singlets, multuplets)
    
    colnames(testSyntheticData) <- testNames
    colnames(syntheticData) <- names
    
    testSyntheticData <- as.matrix(testSyntheticData)
    syntheticData <- as.matrix(syntheticData)
    
    #add new counts and sampleType to uObj
    counts(uObj) <- testSyntheticData
    counts.log(uObj) <- sp.scRNAseq:::.norm.log.counts(testSyntheticData)
    sampleType(uObj) <- c(getData(uObj, "sampleType"), rep("Multuplet", n))
    
    return(list(syntheticData, testSyntheticData, uObj, table))
}

.syntheticSinglets <- function(ngenes, ncells, cellTypes) {

    for( i in 1:cellTypes) {
        set.seed(i)
        meanExprs <- 2^runif(ngenes, 0, 5)
        counts <- matrix(rnbinom(ngenes*ncells, mu=meanExprs, size=i), nrow=ngenes)
        if( i == 1 ) {
            singlets <- counts
        } else {
            singlets <- cbind(singlets, counts)
        }
    }
    colnames(singlets) <- paste(sort(rep(letters, ncells))[1:(cellTypes*ncells)], 1:ncells, sep="")
    singlets <- as.data.frame(singlets)
    
    return(singlets)
}

.syntheticMultuplets <- function(ngenes, ncells, cellTypes, perplexity) {
    singlets <- .syntheticSinglets(ngenes, ncells, cellTypes)
    cObj <- spCounts(as.matrix(singlets), counts.ercc=matrix(), sampleType="[A-Z]")
    uObj <- spUnsupervised(cObj, max=ngenes, max_iter = 1000, perplexity=perplexity)
    colnames(singlets) <- getData(uObj, "classification")
    
    mean <- getData(uObj, "groupMeans")
    
    #same cell type
    combos <- t(matrix(rep(unique(colnames(singlets)), 2), ncol=2))
    
    for(y in 1:ncol(combos)) {
        current <- combos[ ,y]
        new <- data.frame(
            rowMeans(
                data.frame(
                    mean[ , colnames(mean) %in% current[1]],
                    mean[ , colnames(mean) %in% current[2]]
                )
            )
        )
        
        if( y == 1 ) {
            multuplets <- new
            names <- paste(current, collapse="")
        } else {
            multuplets <- cbind(multuplets, new)
            names <- c(names, paste(current, collapse=""))
        }
    }
    
    #doublets
    combos <- combn(unique(colnames(singlets)), 2)
    
    for(y in 1:ncol(combos)) {
        current <- combos[ ,y]
        new <- data.frame(rowMeans(mean[ , colnames(mean) %in% current]))
        multuplets <- cbind(multuplets, new)
        names <- c(names, paste(current, collapse=""))
    }
    
    #triplets
    combos <- combn(unique(colnames(singlets)), 3)
    
    for(u in 1:ncol(combos)) {
        current <- combos[ ,u]
        new <- data.frame(rowMeans(mean[ , colnames(mean) %in% current]))
        multuplets <- cbind(multuplets, new)
        names <- c(names,  paste(current, collapse=""))
    }
    
    #quadruplets
    combos <- combn(unique(colnames(singlets)), 4)
    
    for(u in 1:ncol(combos)) {
        current <- combos[ ,u]
        new <- data.frame(rowMeans(mean[ , colnames(mean) %in% current]))
        multuplets <- cbind(multuplets, new)
        names <- c(names,  paste(current, collapse=""))
    }
    
    colnames(multuplets) <- names
    
    return(list(singlets, multuplets, uObj))
}





