
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
#' @param distFun The distance function to use for swarm optimization.
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

syntheticDataTest <- function(
    outPath = './data',
    cores = 6,
    n = 20,
    ngenes = 2000,
    ncells = 100,
    cellTypes = 10,
    distFun = bic,
    ...
){
    
    #create synthetic data
    tmp <- .syntheticTestData(n, ngenes, ncells, cellTypes)
    syntheticData <- tmp[[1]]
    syntheticDataTest <- tmp[[2]]
    syntheticDataUnsupervised <- tmp[[3]]
    syntheticDataTable <- tmp[[4]]
    
    #make spCounts object
    single <- grepl("^s", colnames(syntheticDataTest))
    syntheticDataCountsSng <- spCounts(
        syntheticDataTest[,single],
        matrix(NA, ncol = length(single[single == TRUE]))
    )
    syntheticDataCountsMul <- spCounts(
        syntheticDataTest[,!single],
        matrix(NA, ncol = length(single[single == FALSE]))
    )

    #run pySwarm
    syntheticDataSwarm <- spSwarm(
        syntheticDataCountsMul,
        syntheticDataUnsupervised,
        maxiter = 1000,
        swarmsize = 500,
        cores = cores,
        distFun = bic,
        report = TRUE,
        reportRate = 10
    )
    
    #save
    save(
        syntheticData,
        file = paste(outPath, "syntheticData.rda", sep = "/"),
        compress = "bzip2"
    )
    save(
        syntheticDataTest,
        syntheticDataUnsupervised,
        syntheticDataSwarm,
        syntheticDataTable,
        file = paste(outPath, "syntheticDataTest.rda", sep = "/"),
        compress = "bzip2"
    )
}

.syntheticTestData <- function(
    n,
    ngenes,
    ncells,
    cellTypes,
    perplexity = 10
){
    tmp <- .syntheticMultuplets(ngenes, ncells, cellTypes, perplexity)
    singlets <- tmp[[1]]
    multuplets <- tmp[[2]]
    uObj <- tmp[[3]]
    
    #select multuplets for current test
    testMultuplets <- multuplets[ ,sample(1:ncol(multuplets), n, replace = FALSE)]
    table <- calculateConnections(testMultuplets, type = "multuplets")
    
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
    #counts(uObj) <- testSyntheticData
    #counts.log(uObj) <- sp.scRNAseq:::.norm.log.counts(testSyntheticData)
    #sampleType(uObj) <- c(getData(uObj, "sampleType"), rep("Multuplet", n))
    
    return(list(syntheticData, testSyntheticData, uObj, table))
}

.syntheticSinglets <- function(
    ngenes,
    ncells,
    cellTypes
){

    for(i in 1:cellTypes) {
        set.seed(i)
        meanExprs <- 2^runif(ngenes, 0, 5)
        counts <- matrix(
            rnbinom(ngenes * ncells, mu = meanExprs, size = i),
            nrow=ngenes
        )
        if( i == 1 ) {
            singlets <- counts
        } else {
            singlets <- cbind(singlets, counts)
        }
    }
    colnames(singlets) <- paste(
        sort(rep(letters, ncells))[1:(cellTypes * ncells)],
        1:ncells,
        sep = ""
    )
    singlets <- as.data.frame(singlets)
    
    return(singlets)
}

.syntheticMultuplets <- function(
    ngenes,
    ncells,
    cellTypes,
    perplexity
){
    singlets <- .syntheticSinglets(ngenes, ncells, cellTypes)
    cObjSng <- spCounts(
        as.matrix(singlets),
        counts.ercc = matrix(
            rep(NA, ncol(singlets)),
            ncol = ncol(singlets)
        )
    )
    uObj <- spUnsupervised(
        cObjSng,
        theta = 0,
        k = 2,
        max_iter = 2000,
        perplexity = 10,
        initial_dims = ncol(singlets),
        Gmax = 50,
        seed = 11,
        type = "max"
    )
    
    colnames(singlets) <- getData(uObj, "classification")
    #mean <- getData(uObj, "groupMeans")
    
    #doublets
    cellNames <- unique(colnames(singlets))
    tmp <- .makeMultuplet(2, cellNames, data.frame(), data.frame(), singlets)
    multuplets <- tmp[[1]]
    names <- tmp[[2]]
    
    #triplets
    tmp <- .makeMultuplet(3, cellNames, multuplets, names, singlets)
    multuplets <- tmp[[1]]
    names <- tmp[[2]]
    
    #quadruplets
    tmp <- .makeMultuplet(4, cellNames, multuplets, names, singlets)
    multuplets <- tmp[[1]]
    names <- tmp[[2]]
    
    colnames(multuplets) <- names
    
    return(list(singlets, multuplets, uObj))
}

.makeMultuplet <- function(
    n,
    cellNames,
    multuplets,
    names,
    singlets
){
    switch(n - 1,
        {combos <- expand.grid(cellNames, cellNames)},
        {combos <- expand.grid(cellNames, cellNames, cellNames)},
        {combos <- expand.grid(cellNames, cellNames, cellNames, cellNames)}
    )
    
    dat.sort <- t(apply(combos, 1, sort))
    combos <- t(combos[!duplicated(dat.sort),])
    
    for(u in 1:ncol(combos)) {
        current <- combos[ ,u]
        
        set.seed(2918834)
        idxs <- c()
        for(i in 1:length(current)) {
            idxs[i] <- sample(
                which(colnames(singlets) == current[1]),
                size = 1,
                replace = FALSE
            )
        }
        
        new <- data.frame(rowMeans(singlets[, idxs]))
        
        if(ncol(multuplets) == 0) {
            multuplets <- new
            names <- paste(current, collapse = "")
        } else {
            multuplets <- cbind(multuplets, new)
            names <- c(names,  paste(current, collapse = ""))
        }
    }
    return(list(multuplets, names))
}






.adjustFreq <- function(ncells, cellNames, multuplets) {
    intMat <- .interactionMatrix(ncells, cellNames, 23443)
    
    x <- as.list(gsub("(.{2})", "\\1 ", colnames(multuplets)))
    l <- lapply(x, function(y) strsplit(y, " "))
    c <- lapply(l, function(y) sum(c("A1", "B1") %in% y[[1]])==2)
    n <- lapply(l, function(y) length(grep("A1", y[[1]])))
    
}

#try to write a fucntion to calculate the frequency of interactions from the colnames of a dataframe with multuplets
.countsInteracts <- function(mul, cellNames) {
    names <- colnames(mul)
    x <- as.list(gsub("(.{2})", "\\1 ", colnames(multuplets)))
    l <- lapply(x, function(y) strsplit(y, " "))
    c <- lapply(l, function(j) combn(j[[1]], 2))
    p <- lapply(c, function(m) sapply(1:ncol(m), function(i) paste(m[,i], collapse="-")))
    t <- table(unlist(p))
    
    m <- matrix(NA, ncol=length(cellNames), nrow=length(cellNames), dimnames=list(cellNames, cellNames))
    
    for(i in 1:length(t)) {
        n1 <- strsplit(names(t)[i], "-")[[1]][1]
        n2 <- strsplit(names(t)[i], "-")[[1]][2]
        m[n1, n2] <- t[i]

    }
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
}

.interactionMatrix <- function(ncells, cellNames, seed) {
    set.seed(seed)
    m <- matrix(NA, ncol=ncells, nrow=ncells, dimnames=(list(cellNames, cellNames)))
    diag(m) <- sample(seq(45, 55, 1), size=length(diag(m)))

    .f <- function(x) {
        s <- sum(x, na.rm=TRUE)
        int <- sample(c(1,2), 1)
        x[sample(which(is.na(x) == TRUE), size=int)] <- sample(seq(30, (100-s), 1), size=int)
        x <- x/sum(x, na.rm=TRUE)*90
        x[is.na(x) == TRUE] <- sample(seq(0, 10, 0.1), size=length(x[is.na(x) == TRUE]))
        return(round(x/sum(x)*100))
    }
    
    apply(m, 1, .f)
}
