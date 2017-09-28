
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
#' @param nMultiplets Number of multiplets.
#' @param nGenes Number of genes in generated synthetic data.
#' @param nCells Number of cells per cell type in generated synthetic data.
#' @param nCellTypes Number of cell types in generated synthetic data. For the
#'    frequency adjustment to work this should be an even number.
#' @param distFun The distance function to use for swarm optimization.
#' @param target The target percentage of connections for interacting cells.
#' @param singletExpansion This argument is multiplied by the nCells argument to
#'    give the total number of singlets synthesized per cell type. nCells
#'    singlets are then reserved for the singlets used in testing whereas the
#'    extra singlets created are used to adjust the connection frequencies via
#'    the .adjustMultuplets and .adjustSelf units. Argument must be > 1.
#' @param perplexity The perplexity provided to Rtsne during the spUnsupervised
#'    stage.
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @examples
#'
#' #use demo data
#' \dontrun{syntheticDataTest()}
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq

################################################################################
#                                                                              #
# Design a method to track the singlets to synthesize multiplets and a plot    #
#    allowing tracking to make sure there are enough singlets in the base pool #
#    and that they are not getting re-used many times.                         #
#                                                                              #
################################################################################

syntheticDataTest <- function(
    outPath = './data',
    cores = 6,
    nMultiplets = 500,
    nGenes = 2000,
    nCells = 100,
    nCellTypes = 10,
    distFun = "bic",
    target = 20,
    singletExpansion = 20,
    perplexity = 10,
    ...
){
    
    #create synthetic data
    tmp <- syntheticTestData(
        nMultiplets,
        nGenes,
        nCells,
        nCellTypes,
        target,
        singletExpansion,
        perplexity
    )
    syntheticData <- tmp[[1]]
    syntheticDataTest <- tmp[[2]]
    syntheticDataUnsupervised <- tmp[[3]]
    syntheticDataTable <- tmp[[4]]
    
    
    #make spCounts object
    single <- grepl("^s", colnames(syntheticDataTest))
    syntheticDataCountsSng <- spCounts(
        syntheticDataTest[, single],
        matrix(NA, ncol = length(single[single == TRUE]))
    )
    
    #fix colnames (repeat colnames on multiplets won't work with the spSwarm
    #method)
    mul <- syntheticDataTest[, !single]
    tmpNames <- paste(colnames(mul), 1:nMultiplets, sep = "_")
    colnames(mul) <- tmpNames
    syntheticDataCountsMul <- spCounts(
        mul,
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
        report = FALSE,
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

#' syntheticTestData
#'
#' This unit and the corresponding units synthesize the synthetic test data
#' that is subsequently saved and, as well, is used for testing.
#'
#'
#' @name syntheticTestData
#' @aliases syntheticTestData
#' @param nMultiplets Number of multiplets.
#' @param nGenes Number of genes in generated synthetic data.
#' @param nCells Number of cells per cell type in generated synthetic data.
#' @param nCellTypes Number of cell types in generated synthetic data. For the
#'    frequency adjustment to work this should be an even number.
#' @param perplexity The perplexity parameter supplied to spUnsupervised.
#' @param target The target percentage of connections for interacting cells.
#' @param singletExpansion This argument is multiplied by the nCells argument to
#'    give the total number of singlets synthesized per cell type. nCells
#'    singlets are then reserved for the singlets used in testing whereas the
#'    extra singlets created are used to adjust the connection frequencies via
#'    the .adjustMultuplets and .adjustSelf units. Argument must be > 1.
#' @param ... additional arguments to pass on
#' @return A list containing variables syntheticData (all data generated),
#'    testSyntheticData (data to be used in the testing), uObj (spUnsupervised
#'    output using the generated singlets), table (a table of the connections
#'    synthesized).
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq

syntheticTestData <- function(
    nMultiplets,
    nGenes,
    nCells,
    nCellTypes,
    target,
    singletExpansion,
    perplexity,
    ...
){
    tmp <- syntheticMultuplets(
        nGenes,
        nCells,
        nCellTypes,
        perplexity,
        target,
        singletExpansion
    )
    singlets <- tmp[[1]]
    multiplets <- tmp[[2]]
    uObj <- tmp[[3]]
    
    #select multiplets for current test
    idx <- sample(1:ncol(multiplets), nMultiplets, replace = FALSE)
    testMultuplets <- multiplets[, idx]
    
    #name, calculate connections, bind, and return
    colnames(multiplets) <- gsub(
        "^([A-Z0-9]*)\\..*",
        "\\1",
        colnames(multiplets)
    )
    colnames(testMultuplets) <- gsub(
        "^([A-Z0-9]*)\\..*",
        "\\1",
        colnames(testMultuplets)
    )
    
    table <- calculateConnections(testMultuplets, type = "multiplets")
    
    names <- c(
        paste("s", colnames(singlets), sep = "."),
        paste("m", colnames(multiplets), sep = ".")
    )
    
    testNames <- c(
        paste("s", colnames(singlets), sep = "."),
        paste("m", colnames(testMultuplets), sep = ".")
    )
    
    testSyntheticData <- cbind(singlets, testMultuplets)
    syntheticData <- cbind(singlets, multiplets)
    
    colnames(testSyntheticData) <- testNames
    colnames(syntheticData) <- names
    
    testSyntheticData <- as.matrix(testSyntheticData)
    syntheticData <- as.matrix(syntheticData)
    
    return(list(syntheticData, testSyntheticData, uObj, table))
}

#' syntheticSinglets
#'
#' This unit uses the negative binomial distribution to synthesize the
#' singlets used in the synthetic data test.
#'
#'
#' @name syntheticSinglets
#' @aliases syntheticSinglets
#' @param nGenes Number of genes in generated synthetic data.
#' @param nCells Number of cells per cell type in generated synthetic data.
#' @param nCellTypes Number of cell types in generated synthetic data. For the
#'    frequency adjustment to work this should be an even number.
#' @param singletExpansion This argument is multiplied by the nCells argument to
#'    give the total number of singlets synthesized per cell type. nCells
#'    singlets are then reserved for the singlets used in testing whereas the
#'    extra singlets created are used to adjust the connection frequencies via
#'    the .adjustMultuplets and .adjustSelf units. Argument must be > 1.
#' @param ... additional arguments to pass on
#' @return A list containing the singlets generated for connection adjustments
#'    and the singlets for use with spUnsupervised in the testing procedure.
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq

syntheticSinglets <- function(
    nGenes,
    nCells,
    nCellTypes,
    singletExpansion,
    ...
){
    nCellsIncrease <- nCells * singletExpansion
    synth <- sapply(1:nCellTypes, function(x) {
        set.seed(x)
        rnbinom(nGenes * nCellsIncrease, mu = 2^runif(nGenes, 0, 5), size = x)
    })
    
    singlets <- matrix(as.numeric(synth), nrow = nGenes)
    
    colnames(singlets) <- paste(
        sort(rep(LETTERS, nCellsIncrease))[1:(nCellTypes * nCellsIncrease)],
        1:nCellsIncrease,
        sep = ""
    )
    
    idx <- sapply(1:nCellTypes, function(t) {
        typeIdx <- which(grepl(LETTERS[t], colnames(singlets)))
        sample(typeIdx, size = nCells, replace = FALSE)
    })
    singletsProper <- as.data.frame(singlets[, as.numeric(idx)])
    singletsFull <- as.data.frame(singlets[, -as.numeric(idx)])
    return(list(singletsFull, singletsProper))
}

#' syntheticMultuplets
#'
#' This unit synthesizes "all possible" multiplets combinations from the
#' singlets that were generated. Currently, doublets, triplets, and quadruplets
#' are synthesized. The unit also calls the adjustMultuplets and adjustSelf
#' units to adjust the connections in a manner that 1) self connections are the
#' most prevelant connection type (currently 2/3 connections) and 2) that each
#' cell type has a high percentage of connections to one other cell type.
#'
#' @name syntheticMultuplets
#' @aliases syntheticMultuplets
#' @param nGenes numeric. Number of genes in generated synthetic data.
#' @param nCells numeric. Number of cells per cell type in generated synthetic
#'    data.
#' @param nCellTypes numeric. Number of cell types in generated synthetic data.
#'    For the frequency adjustment to work this should be an even number.
#' @param perplexity numeric. The perplexity parameter supplied to
#'    spUnsupervised.
#' @param target numeric. The target percentage of connections for interacting
#'    cells.
#' @param singletExpansion numeric. This argument is multiplied by the nCells
#'    argument to give the total number of singlets synthesized per cell type.
#'    nCells singlets are then reserved for the singlets used in testing whereas
#'    the extra singlets created are used to adjust the connection frequencies
#'    via the .adjustMultuplets and .adjustSelf units. Argument must be > 1.
#' @param ... additional arguments to pass on
#' @return A list containing singlets, multiplets, and the spUnsupervised
#'    object.
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq


syntheticMultuplets <- function(
    nGenes,
    nCells,
    nCellTypes,
    perplexity,
    target,
    singletExpansion,
    ...
){
    sng <- syntheticSinglets(nGenes, nCells, nCellTypes, singletExpansion)
    singletsFull <- sng[[1]]
    singlets <- sng[[2]]
    
    cObjSngFull <- spCounts(
        as.matrix(singletsFull),
        counts.ercc = matrix(
            rep(NA, ncol(singletsFull)),
            ncol = ncol(singletsFull)
        )
    )
    
    cObjSng <- spCounts(
        as.matrix(singlets),
        counts.ercc = matrix(
            rep(NA, ncol(singlets)),
            ncol = ncol(singlets)
        )
    )
    
    uObjFull <- spUnsupervised(
        cObjSngFull,
        theta = 0,
        k = 2,
        max_iter = 500,
        perplexity = 10,
        initial_dims = ncol(singlets),
        Gmax = nCellTypes,
        seed = 11,
        type = "max"
    )
    
    uObj <- spUnsupervised(
        cObjSng,
        theta = 0,
        k = 2,
        max_iter = 500,
        perplexity = 10,
        initial_dims = ncol(singlets),
        Gmax = nCellTypes,
        seed = 11,
        type = "max"
    )
    
    colnames(singletsFull) <- getData(uObjFull, "classification")
    colnames(singlets) <- getData(uObj, "classification")

    cellTypes <- unique(colnames(singletsFull))
    tmp <- data.frame(row.names = 1:nGenes)
    
    #doublets
    doublets <- makeMultuplet(2, cellTypes, tmp, singletsFull)
    
    #triplets
    triplets <- makeMultuplet(3, cellTypes, doublets, singletsFull)
    
    #quadruplets
    quadruplets <- makeMultuplet(4, cellTypes, triplets, singletsFull)
    
    #adjust perfered connections
    multiplets <- adjustMultuplets(
        singletsFull,
        quadruplets,
        nGenes,
        nCellTypes,
        target
    )
    
    #adjust self connections
    multiplets <- adjustSelf(multiplets, singletsFull)
    
    return(list(singlets, multiplets, uObj))
}

