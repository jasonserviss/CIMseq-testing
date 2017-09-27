
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
    multuplets <- tmp[[2]]
    uObj <- tmp[[3]]
    
    #select multuplets for current test
    idx <- sample(1:ncol(multuplets), nMultiplets, replace = FALSE)
    testMultuplets <- multuplets[, idx]
    
    #name, calculate connections, bind, and return
    colnames(multuplets) <- gsub(
        "^([A-Z0-9]*)\\..*",
        "\\1",
        colnames(multuplets)
    )
    colnames(testMultuplets) <- gsub(
        "^([A-Z0-9]*)\\..*",
        "\\1",
        colnames(testMultuplets)
    )
    
    table <- calculateConnections(testMultuplets, type = "multuplets")
    
    names <- c(
        paste("s", colnames(singlets), sep = "."),
        paste("m", colnames(multuplets), sep = ".")
    )
    
    testNames <- c(
        paste("s", colnames(singlets), sep = "."),
        paste("m", colnames(testMultuplets), sep = ".")
    )
    
    testSyntheticData <- cbind(singlets, testMultuplets)
    syntheticData <- cbind(singlets, multuplets)
    
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
#' @return A list containing singlets, multuplets, and the spUnsupervised
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
    multuplets <- adjustMultuplets(
        singletsFull,
        quadruplets,
        nGenes,
        nCellTypes,
        target
    )
    
    #adjust self connections
    multuplets <- adjustSelf(multuplets, singletsFull)
    
    return(list(singlets, multuplets, uObj))
}

#' makeMultuplet
#'
#' This unit synthesizes "all possible" multiplets combinations from the
#' singlets that were generated. Currently, doublets, triplets, and quadruplets
#' are synthesized. The unit also calls the adjustMultuplets and adjustSelf
#' units to adjust the connections in a manner that 1) self connections are the
#' most prevelant connection type (currently 2/3 connections) and 2) that each
#' cell type has a high percentage of connections to one other cell type.
#'
#' @name makeMultuplet
#' @aliases makeMultuplet
#' @param nCellsInMultiplet numeric. Defines the number of cells to be included
#'    in the multiplet. Currently can only be 2, 3, or 4.
#' @param cellTypes character. Provides the names of the cell types to include
#'    in the multiplet(s). These must be present in the singlets matrix and the
#'    names specified correctly in the colnames.
#' @param multuplets data.frame. The multiplets data.frame or a empty data.frame
#'    with rows equal to the number of genes in the desired multiplets variable.
#' @param singlets data.frame The singlets variable output from the
#'    syntheticSinglets function.
#' @param repetitions numeric. The number of times the connections should be
#'    repeated or NULL if they should only be repeated once.
#' @param self logical. Indicates if self connections should be included.
#' @param ... additional arguments to pass on
#' @return The multuplets variable.
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq


########################
#NOTE!!!
#add an option to weight the contribution of individual cell types in the multiplet!!
#######################

makeMultuplet <- function(
    nCellsInMultiplet,
    cellTypes,
    multuplets,
    singlets,
    repetitions = NULL,
    self = TRUE,
    ...
){
    switch(nCellsInMultiplet - 1,
        {combos <- t(expand.grid(cellTypes, cellTypes))},
        {combos <- t(expand.grid(cellTypes, cellTypes, cellTypes))},
        {combos <- t(expand.grid(cellTypes, cellTypes, cellTypes, cellTypes))}
    )
    
    dat.sort <- apply(combos, 2, sort)
    combos <- matrix(dat.sort[, !duplicated(t(dat.sort))], nrow = nrow(combos))
    
    if(!self) {
        bools <- !sapply(1:ncol(combos), function(c)
            length(unique(combos[, c])) == 1
        )
        combos <- matrix(unlist(combos[, bools]), nrow = nrow(combos))
    }
    
    if(!is.null(repetitions)) {
        combos <- matrix(
            combos[, rep(seq_len(ncol(combos)), repetitions)],
            nrow = nrow(combos)
        )
    }
    
    add <- sapply(1:ncol(combos), function(x) {
        idxs <- sapply(combos[, x], function(y) which(colnames(singlets) == y))
        set.seed(x)
        if(class(idxs) == "integer") {
            pick <- sapply(1:length(idxs), function(x)
                sample(idxs[x], size = 1
            ))
        } else {
            pick <- sapply(1:ncol(idxs), function(x)
                sample(idxs[, x], size = 1
            ))
        }
        
        rowMeans(singlets[, pick])
    })
    
    colnames(add) <- sapply(1:ncol(combos), function(i) {
        paste(combos[, i], collapse = "")
    })
    
    multuplets <- cbind(multuplets, add)
    return(multuplets)
}

#' adjustMultuplets
#'
#' This unit adjusts the multuplets so that each cell type has one prevelant
#' connection. Which connection a cell type prefers is decided randomly. The
#' strength of the preference is determined by the target argument. This works
#' well as long as the number of cell types in the dataset is even (and not
#' odd).
#'
#' @name adjustMultuplets
#' @aliases adjustMultuplets
#' @param singletsFull data.frame The singlets variable output from the
#'    syntheticSinglets function.
#' @param multuplets data.frame. The multiplets data.frame or a empty data.frame
#'    with rows equal to the number of genes in the desired multiplets variable.
#' @param nGenes numeric. Number of genes in generated synthetic data.
#' @param nCellTypes numeric. Number of cell types in generated synthetic data.
#'    For the frequency adjustment to work this should be an even number.
#' @param perplexity numeric. The perplexity parameter supplied to
#'    spUnsupervised.
#' @param target numeric. The target percentage of connections for interacting
#'    cells.
#' @param ... additional arguments to pass on
#' @return The multuplets variable.
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq

adjustMultuplets <- function(
    singletsFull,
    multuplets,
    nGenes,
    nCellTypes,
    target,
    ...
){
    
    #decide which connection preferences should exist
    cellTypes <- unique(colnames(singletsFull))
    targetConnections <- decideConnections(cellTypes, target)
    
    #quantify current connections
    current <- quantifyConnections(colnames(multuplets))
    
    #adjust connection frequency
    
    for(i in 1:nrow(targetConnections)) {
        
        #calculate current percent
        tmp <- calculateFreqAndSum(i, targetConnections, current)
        currentTypes <- tmp[[1]]
        currentFreq <- tmp[[2]]
        currentSum <- tmp[[3]]
        
        #calculate number to add
        add <- calculateNumToAdd(i, currentFreq, currentSum, targetConnections)
        
        #synthesize and add connections to achieve target percentages
        multuplets <- synthesizeAndAdd(
            add,
            currentTypes,
            multuplets,
            singletsFull
        )
        
        #check frequency
        targetConnections <- checkFrequency(
            i,
            multuplets,
            targetConnections,
            currentTypes
        )
    }
    
    #check that targets were met
    sumD <- sum(targetConnections$target - targetConnections$current)
    okError <- sum(nrow(targetConnections) * 5)
    sumThigh <- okError
    sumTlow <- okError * -1
    
    if(!sumD > sumThigh & sumD < sumTlow) {
        stop("target connection percentage was not met. check code")
    }
    
    return(multuplets)
}

#' decideConnections
#'
#' This randomly assigns connections between all cell types and returns a
#' data.frame to allow tracking of the target and current connection percentage
#' per connection.
#'
#' @name decideConnections
#' @aliases decideConnections
#' @param cellTypes character. Provides the names of the cell types to include
#'    in the multiplet(s). These must be present in the singlets matrix and the
#'    names specified correctly in the colnames.
#' @param target numeric. The target percentage of connections for interacting
#'    cells.
#' @param ... additional arguments to pass on.
#' @return data.frame
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq

decideConnections <- function(
    cellTypes,
    target,
    ...
){
    
    done <- c()
    count <- 0
    combos <- c()
    for(i in 1:(length(cellTypes) / 2)) {
        set.seed(988797 + count)
        comb <- sort(sample(cellTypes, size = 2, replace = FALSE))
        while(any(comb %in% done)) {
            set.seed(988797 + count)
            comb <- sort(sample(cellTypes, size = 2, replace = FALSE))
            count <- count + 1
        }
        done <- c(done, comb)
        combos[i] <- paste(comb[1], comb[2], sep = "-")
    }
    
    if(any(!cellTypes %in% unlist(strsplit(combos, "-")))) {
        idx <- which(!cellTypes %in% unlist(strsplit(combos, "-")))
        missing <- cellTypes[idx]
        combos <- c(combos, paste(missing, sample(cellTypes, 1), sep = "-"))
    }
    
    return(data.frame(
        conn = combos,
        target = target,
        current = 0,
        stringsAsFactors = FALSE
    ))
}

#' quantifyConnections
#'
#' Uses the multuplet naming convention to quantify the frequency of each
#' connection.
#'
#' @name quantifyConnections
#' @aliases quantifyConnections
#' @param multupletNames character. Names of the multiplets specified using the
#'    typical LETTER NUMBER convention.
#' @param ... additional arguments to pass on.
#' @return data.frame
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq
#' @export

quantifyConnections <- function(
    multipletNames,
    ...
){
    adjNames <- .adjustNames(multipletNames)
    combs <- .findCombos(adjNames)
    .buildOutput(combs)
}

.adjustNames <- function(
    multupletNames,
    ...
){
    suffix <- gsub("^[m-s]\\.(.*)", "\\1", multupletNames)
    prefix <- gsub("^(.*)_[0-9]*", "\\1", suffix)
    split <- trimws(gsub("(.{2})", "\\1 ", prefix))
    suffixRemove <- gsub("^([A-Z0-9]* [A-Z0-9]*) \\..*$", "\\1", split)

}

.findCombos <- function(
    adjNames,
    ...
) {
    ss <- strsplit(adjNames, " ")
    l <- lapply(ss, function(x) {
        if(length(unique(x)) > 1) {
            combn(unique(x), 2)
        }
    })
    self <- lapply(ss, function(x) which(duplicated(x)))
    add <- lapply(1:length(ss), function(o) unique(ss[[o]][self[[o]]]))
    complete <- lapply(1:length(add), function(x) {
        if(length(add[[x]]) == 0) {
            l[[x]]
        } else if(length(add[[x]]) == 1) {
            cbind(l[[x]], matrix(add[[x]], nrow = 2))
        } else {
            cbind(
                l[[x]],
                matrix(rep(add[[x]], length(add[[x]])), nrow = 2, byrow = TRUE)
            )
        }
    })
    srt <- lapply(complete, function(x) apply(x, 2, sort))
    lapply(srt, function(x) paste(x[1, ], x[2, ], sep = "-"))
}

.buildOutput <- function(
    combs,
    ...
) {
    d <- as.data.frame(table(unlist(combs)), stringsAsFactors = FALSE)
    d$type1 <- gsub("^(..)...$", "\\1", d$Var1)
    d$type2 <- gsub("^...(..)$", "\\1", d$Var1)
    d
}

#' calculateFreqAndSum
#'
#' Calculates statistics concerning the current connection. These are used
#' downstream to calculate the number of connections that need to be added to
#' meet the target percentage.
#'
#' @name calculateFreqAndSum
#' @aliases calculateFreqAndSum
#' @param i integer. Current iteration.
#' @param targetConnections data.frame. Output from the decideConnections unit.
#' @param current data.frame. Output from the quantifyConnections unit.
#' @param ... additional arguments to pass on.
#' @return A list containing the current cell types, the current frequency of
#'    of the prefered connection, and the sum of the current number of
#'    connections between the cell types.
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq

calculateFreqAndSum <- function(
    i,
    targetConnections,
    current
){
    type1 <- gsub("^(..)...$", "\\1", targetConnections[i, "conn"])
    type2 <- gsub("^...(..)$", "\\1", targetConnections[i, "conn"])
    types <- c(type1, type2)
    
    bool <- current$Var1 == targetConnections[i, "conn"]
    currentFreq <- current[bool, "Freq"]
    bool <- current$type1 %in% types | current$type2 %in% types
    currentSum <- sum(current[bool, "Freq"])
    return(list(types, currentFreq, currentSum))
}

#' calculateNumToAdd
#'
#' Calculates the number of connections that need to be added to meet the target
#' percentage for each connection.
#'
#' @name calculateNumToAdd
#' @aliases calculateNumToAdd
#' @param i integer. Current iteration.
#' @param currentFreq numeric. The current frequency of of the prefered
#'    connection.
#' @param currentSum numeric.
#' @param targetConnections data.frame. Output from the decideConnections unit.
#' @param ... additional arguments to pass on.
#' @return Integer specifying the number of connections to add in order to meet
#'    the target percentage.
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq

calculateNumToAdd <- function(
    i,
    currentFreq,
    currentSum,
    targetConnections,
    ...
){
    target <- targetConnections[i, "target"]
    numerator <- (100 * currentFreq) - (currentSum * target)
    denominator <- targetConnections[i, "target"] - 100
    add <- ceiling(numerator / denominator)
    mess <- paste(
        "Adjusting connections percentages by adding ", add, " connections.",
        sep = ""
    )
    print(mess)
    return(add)
}

#' synthesizeAndAdd
#'
#' Synthesizes the number of connections necessary to meet the target
#' percentages according to the add argument. Doublets, triplets, and
#' quadruplets are synthesized.
#'
#' @name synthesizeAndAdd
#' @aliases synthesizeAndAdd
#' @param add integer. The number of connections to add.
#' @param currentTypes data.frame. Output from the quantifyConnections unit.
#' @param multuplets data.frame. The multiplets data.frame or a empty data.frame
#'    with rows equal to the number of genes in the desired multiplets variable.
#' @param singletsFull data.frame The singlets variable output from the
#'    syntheticSinglets function.
#' @param ... additional arguments to pass on.
#' @return The multiplets variable with the additional multiplets added by the
#'    unit.
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq

synthesizeAndAdd <- function(
    add,
    currentTypes,
    multuplets,
    singletsFull,
    ...
){
    tmp <- data.frame(row.names = 1:nrow(multuplets))
    doublets <- makeMultuplet(
        nCellsInMultiplet = 2,
        cellTypes = currentTypes,
        multuplets = tmp,
        singlets = singletsFull,
        repetitions = add,
        self = FALSE
    )
    
    triplets <- makeMultuplet(
        nCellsInMultiplet = 3,
        cellTypes = currentTypes,
        multuplets = tmp,
        singlets = singletsFull,
        repetitions = add,
        self = FALSE
    )
    
    quadruplets <- makeMultuplet(
        nCellsInMultiplet = 4,
        cellTypes = currentTypes,
        multuplets = tmp,
        singlets = singletsFull,
        repetitions = add,
        self = FALSE
    )
    
    #add random quadruplets
    set.seed(2938230)
    quadrupletIdx <- sample(
        1:ncol(quadruplets),
        size = ceiling(add / 6),
        replace = FALSE
    )
    mat <- matrix(unlist(quadruplets[, quadrupletIdx]), nrow = nrow(multuplets))
    names <- colnames(quadruplets)[quadrupletIdx]
    
    #quantify quadruplet connections
    qConn <- quantifyConnections(names)
    bools <- qConn$Var1 == paste(currentTypes, collapse = "-")
    added <- qConn[bools, "Freq"]
    
    #fill remaining connections using doublets and triplets
    if(!added >= add) {
        a <- ceiling((add - added) / 2)
        
        set.seed(938434)
        doubletIdx <- sample(1:ncol(doublets), size = a, replace = FALSE)
        set.seed(872390)
        tripletIdx <- sample(1:ncol(triplets), size = a, replace = FALSE)
        
        mat <- cbind(mat, doublets[, doubletIdx])
        mat <- cbind(mat, triplets[, tripletIdx])
        
        names <- c(names, colnames(doublets)[doubletIdx])
        names <- c(names, colnames(triplets)[tripletIdx])
        colnames(mat) <- names
        
    } else {
        stop("quadruplets added too many connections in .synthesizeAndAdd")
    }
    
    return(cbind(multuplets, mat))
}

#' checkFrequency
#'
#' Performs a final check to see if the target frequencies have been met. The
#' current percentages are added to the targetConnections data.frame which is
#' subsequently checked outside of the unit.
#'
#' @name checkFrequency
#' @aliases checkFrequency
#' @param i integer. Current iteration.
#' @param currentTypes data.frame. Output from the quantifyConnections unit.
#' @param multuplets data.frame. The multiplets data.frame.
#' @param targetConnections data.frame. Output from the decideConnections unit.
#' @param currentTypes The current cell types being examined.
#' @param ... additional arguments to pass on.
#' @return The targetConnections data.frame.
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq

checkFrequency <- function(
    i,
    multuplets,
    targetConnections,
    currentTypes,
    ...
){
    after <- quantifyConnections(colnames(multuplets))
    bool <- after$Var1 == targetConnections[i, "conn"]
    currentFreq <- after[bool, "Freq"]
    bool <- after$type1 %in% currentTypes | after$type2 %in% currentTypes
    currentSum <- sum(after[bool, "Freq"])
    currentPercent <- (currentFreq / currentSum) * 100
    targetConnections[i, "current"] <- round(currentPercent)
    return(targetConnections)
}

#' adjustSelf
#'
#' Adjust the number of self connections. This currently dilutes the current
#' connections into a third of the total connections by adding twice as many
#' self connections.
#'
#' @name adjustSelf
#' @aliases adjustSelf
#' @param multuplets data.frame. The multiplets data.frame.
#' @param singletsFull data.frame The singlets variable output from the
#'    syntheticSinglets function.
#' @param ... additional arguments to pass on.
#' @return The targetConnections data.frame.
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq

adjustSelf <- function(
    multuplets,
    singletsFull
){
    
    #count number of current connections
    conn <- quantifyConnections(colnames(multuplets))
    selfConn <- subset(conn, type1 == type2)
    
    #calculate total connections per cell type
    cellTypes <- unique(c(conn$type1, conn$type2))
    selfConn$totals <- sapply(1:length(cellTypes), function(o) {
        sum(conn[conn$type1 %in% cellTypes[o] | conn$type2 %in% cellTypes[o], "Freq"])
    })
    
    #calculate target
    selfConn$add <- (selfConn$totals * 2) - selfConn$Freq
    
    #calculate connections to add per cell type
    #currently dilutes the perfered connections by 3
    mess <- paste(
        "Adding ",
        sum(selfConn$add),
        " self connections.",
        sep = ""
    )
    print(mess)
    
    self <- data.frame(row.names = 1:nrow(multuplets))
    for(i in 1:nrow(selfConn)) {
        
        #calculate the number of each multuplet type to add
        #doublets = 1 connection
        #triplets = 3 connections
        #quadruplets = 6 connections
        #NOTE THAT THE dToAdd, tToAdd, qToAdd rarely give the exact expected amount of added connections. Fix?
        connPerCellComb <- selfConn[i, "add"] / 3
        dToAdd <- ceiling(connPerCellComb)
        tToAdd <- ceiling(connPerCellComb / 3)
        qToAdd <- ceiling(connPerCellComb / 6)
        
        #generate multuplets
        #doublets
        
        tmp <- data.frame(row.names = 1:nrow(multuplets))
        doublets <- makeMultuplet(
            nCellsInMultiplet = 2,
            cellTypes = c(selfConn[i, "type1"], selfConn[i, "type2"]),
            multuplets = tmp,
            singlets = singletsFull,
            repetitions = dToAdd,
            self = TRUE
        )
        
        triplets <- makeMultuplet(
            nCellsInMultiplet = 3,
            cellTypes = c(selfConn[i, "type1"], selfConn[i, "type2"]),
            multuplets = tmp,
            singlets = singletsFull,
            repetitions = tToAdd,
            self = TRUE
        )
        
        quadruplets <- makeMultuplet(
            nCellsInMultiplet = 4,
            cellTypes = c(selfConn[i, "type1"], selfConn[i, "type2"]),
            multuplets = tmp,
            singlets = singletsFull,
            repetitions = qToAdd,
            self = TRUE
        )
        
        #add
        self <- cbind(self, doublets, triplets, quadruplets, row.names = NULL)
        
        selfConn[i, "current"] <- sum(c(dToAdd, tToAdd * 3, qToAdd * 6))
    }
    
    return(cbind(multuplets, self))
}


################################################################################
#                                                                              #
# Plot distribution of connections per cell type.                              #
#                                                                              #
################################################################################

#' plotConnectionDist
#'
#' Plots the results from the testing of the synthetic dataset.
#'
#' @name plotConnectionDist
#' @aliases plotConnectionDist
#' @param multipletNames character. A character vector of multiplet adhering to
#'    the "m.A1B1_1" convention.
#' @param ... additional arguments to pass on.
#' @return data.frame
#' @author Jason T. Serviss
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename rowwise mutate
#' @importFrom tidyr gather
#' @importFrom ggthemes scale_fill_ptol theme_few
#' @export

plotConnectionDist <- function(
    multupletNames
){
    freq <- quantifyConnections(multupletNames)
    cellTypes <- unique(c(freq$type1, freq$type2))

    percents <- sapply(cellTypes, function(x) {
        curr <- subset(freq, type1 == x | type2 == x)
        sum <- sum(curr$Freq)

        sapply(cellTypes, function(y) {
            if(y == x) {
                100 * (subset(curr, type1 == y & type2 == y)$Freq / sum)
            } else {
                c <- subset(curr, type1 == y | type2 == y)$Freq
                if(length(c) == 0) {
                    0
                } else {
                    100 * (c / sum)
                }
            }
        })
    })
    
    .conn <- function(x, y) {
        data.frame(
            paste(sort(c(x, y)), collapse = "-"),
            stringsAsFactors = FALSE
        )
    }
    
    p <- percents %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            rename(from = rowname) %>%
            gather(to, percent, -from) %>%
            rowwise() %>%
            mutate(connection = paste(sort(c(from, to)), collapse = "-")) %>%
            ggplot(aes(connection, percent)) +
                geom_bar(
                    aes(fill = to),
                    stat = "identity",
                    position = position_dodge(width = 1)
                ) +
                facet_grid(from~to, scales = "free", space = "free") +
                ggthemes::scale_fill_ptol() +
                ggthemes::theme_few() +
                theme(
                    axis.text.x = element_text(angle = 90)
                ) +
                guides(fill = FALSE) +
                labs(
                    x = "Connection",
                    y = "Percent"
                )
    p
    return(p)
}

################################################################################
#                                                                              #
# Plot distribution of doublets, triplets, quadruplets in testing data.        #
#                                                                              #
################################################################################

#' plotMultipletsCellNrs
#'
#'
#' @name plotMultipletsCellNrs
#' @aliases plotMultipletsCellNrs
#' @param multipletNames character. A character vector of multiplet adhering to
#'    the "m.A1B1_1" convention.
#' @param ... additional arguments to pass on.
#' @return data.frame
#' @author Jason T. Serviss
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_fill_ptol
#' @importFrom tibble as_tibble
#' @export

plotMultipletsCellNrs <- function(multipletNames) {
    p <- gsub("^[m-s]\\.(.*)", "\\1", multipletNames) %>%
            gsub("^(.*)_[0-9]*", "\\1", .) %>%
            nchar(.) %>%
            `/`(2) %>%
            table() %>%
            as_tibble() %>%
            setNames(c("Cell number", "n")) %>%
            ggplot(., aes(`Cell number`, n)) +
                geom_bar(stat = "identity", aes(fill = `Cell number`)) +
                theme_few() +
                scale_fill_ptol()
    p
    return(p)
}


################################################################################
#                                                                              #
# Plot expected vs observed connections post-swarm.                            #
#                                                                              #
################################################################################

#' plotExpVSObs
#'
#'
#' @name plotExpVSObs
#' @aliases plotExpVSObs
#' @param spSwarm spSwarm. An spSwarm object.
#' @param edge.cutoff numeric. The minimum fraction to consider (?).
#' @param summary logical. Indicates if edges should be summarized or not.
#' @param ... additional arguments to pass on.
#' @return data.frame
#' @author Jason T. Serviss
#' @examples
#'
#' cat("No example")
#'
NULL
#' @import sp.scRNAseq
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_fill_ptol
#' @importFrom dplyr mutate select rename left_join group_by summarise ungroup
#' @importFrom tidyr gather
#' @importFrom stringr str_extract
#' @export

plotExpVSObs <- function(
    spSwarm,
    edge.cutoff,
    summary = FALSE
){
    p <- spSwarmPoisson(
        spSwarm,
        edge.cutoff = edge.cutoff,
        min.pval = 1,
        min.num.edges = 0
    ) %>%
        mutate(Var1 = paste(from, to, sep = "-")) %>%
        select(Var1, weight) %>%
        rename(Detected = weight) %>%
        left_join(
            quantifyConnections(rownames(getData(spSwarm, "spSwarm"))),
            by = "Var1"
        ) %>%
        select(Var1, Detected, Freq) %>%
        rename(
            Connection = Var1,
            Expected = Freq
        ) %>%
        mutate(Expected = ifelse(is.na(Expected) & Detected == 0, 0, Expected)) %>%
        gather(Phase, n, -Connection) %>%
        {
            
            if(summary) {
                .plotExpVSObs_summary(.)
            } else {
                .plotExpVSObs_noSummary(.)
            }
        }
    
    p
    return(p)
}

.plotExpVSObs_summary <- function(
    data
){
    data %>%
        mutate(`Self connection` = ifelse(
            str_extract(.$Connection, "^..") == str_extract(.$Connection, "..$"),
            TRUE,
            FALSE
        )) %>%
        group_by(., Phase, `Self connection`) %>%
        summarise(Sum = sum(n)) %>%
        ungroup() %>%
        mutate(
            Phase = readr::parse_factor(Phase, c("Expected", "Detected"))
        ) %>%
        ggplot(aes(`Self connection`, Sum)) +
            geom_bar(
                aes(fill = Phase),
                stat = "identity",
                position = position_dodge(width = 1)
            ) +
            geom_label(
                aes(label = Sum, group = Phase),
                vjust = -0.2,
                position = position_dodge(width = 1)
            ) +
            theme_few() +
            scale_fill_ptol() +
            theme(
                axis.text.x = element_text(angle = 90),
                legend.title = element_blank()
            )
}

.plotExpVSObs_noSummary <- function(
    data
){
    data %>%
        ggplot(aes(Connection, n)) +
            geom_bar(
                aes(fill = Phase),
                stat = "identity",
                position = position_dodge(width = 1)
            ) +
            theme_few() +
            scale_fill_ptol() +
            theme(
                axis.text.x = element_text(angle = 90),
                legend.title = element_blank()
            )
}

################################################################################
#                                                                              #
# Get multiplets with undetected connections.                                  #
#                                                                              #
################################################################################

getUndetectedMultiplets <- function(
    spSwarm,
    edge.cutoff
){
    #detected
    getEdgesForMultiplet(
        spSwarm,
        edge.cutoff,
        rownames(getData(spSwarm, "spSwarm"))
    ) %>%
        rowwise() %>%
        mutate(connection = paste(sort(c(from, to)), collapse = "-")) %>%
        ungroup() %>%
        group_by(multiplet, connection) %>%
        summarise(detected = n()) %>%
        ungroup() %>%
        #expected
        full_join(
            .edgesPerMultiplet(rownames(getData(spSwarm, "spSwarm"))) %>%
                group_by(names, variables) %>%
                summarise(expected = n()) %>%
                rename(multiplet = names, connection = variables)
        ) %>%
        #filter(detected != expected) %>%
        mutate(
            from = str_extract(connection, "^.."),
            to = str_extract(connection, "..$"),
            self = ifelse(from == to, TRUE, FALSE)
        ) %>%
        mutate_if(is.numeric, funs(ifelse(is.na(.), 0, .)))

}

#Uses the multiplets names to report the connections expected per multiplet
#in a long form tibble
.edgesPerMultiplet <- function(
    multipletNames
){
    adjNames <- .adjustNames(multipletNames)
    combs <- .findCombos(adjNames)
    names(combs) <- multipletNames
    namedListToTibble(combs)
}

