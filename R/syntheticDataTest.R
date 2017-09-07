
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
#' @param cellTypes Number of cell types in generated synthetic data. For the
#'    frequency adjustment to work this should be an even number.
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
    n = 500,
    ngenes = 2000,
    ncells = 100,
    cellTypes = 10,
    distFun = bic,
    target = 20,
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
    tmp <- .syntheticMultuplets(ngenes, ncells, cellTypes, perplexity, target)
    singlets <- tmp[[1]]
    multuplets <- tmp[[2]]
    uObj <- tmp[[3]]
    
    #select multuplets for current test
    idx <- sample(1:ncol(multuplets), n, replace = FALSE)
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
    synth <- sapply(1:cellTypes, function(x) {
        set.seed(x)
        rnbinom(ngenes * ncells, mu = 2^runif(ngenes, 0, 5), size = i)
    })
    
    singlets <- matrix(as.numeric(synth), nrow = ngenes)
    
    colnames(singlets) <- paste(
        sort(rep(LETTERS, ncells))[1:(cellTypes * ncells)],
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
    perplexity,
    target
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
    cellNames <- unique(colnames(singlets))
    tmp <- data.frame(row.names = 1:ngenes)
    
    #doublets
    tmp <- .makeMultuplet(2, cellNames, tmp, singlets)
    
    #triplets
    tmp <- .makeMultuplet(3, cellNames, tmp, singlets)
    
    #quadruplets
    tmp <- .makeMultuplet(4, cellNames, tmp, singlets)
    
    #adjust perfered connections
    multuplets <- .adjustMultuplets(
        singlets,
        tmp3[[1]],
        ngenes,
        cellTypes,
        target
    )
    
    #adjust self connections
    multuplets <- .adjustSelf(multuplets, singlets)
    
    return(list(singlets, multuplets, uObj))
}


.makeMultuplet <- function(
    n,
    cellNames,
    multuplets,
    singlets
){
    switch(n - 1,
        {combos <- expand.grid(cellNames, cellNames)},
        {combos <- expand.grid(cellNames, cellNames, cellNames)},
        {combos <- expand.grid(cellNames, cellNames, cellNames, cellNames)}
    )
    
    dat.sort <- t(apply(combos, 1, sort))
    combos <- t(combos[!duplicated(dat.sort),])
    
    add <- sapply(1:ncol(combos), function(x) {
        idxs <- sapply(combos[, x], function(y) which(colnames(singlets) == y))
        pick <- sapply(1:ncol(idxs), function(x) sample(idxs[, x], size = 1))
        rowMeans(singlets[, pick])
    })
    
    colnames(add) <- sapply(1:ncol(combos), function(i) {
        paste(combos[, i], collapse = "")
    })
    
    multuplets <- cbind(multuplets, add)
    return(multuplets)
}

.adjustMultuplets <- function(
    singlets,
    multuplets,
    ngenes,
    cellTypes,
    target
){
    
    #decide which connection preferences should exist
    targetConnections <- decideConnections(unique(colnames(singlets)), target)
    
    #quantify current connections
    current <- .quantifyConnections(colnames(multuplets))
    
    #adjust connection frequency
    
    for(i in 1:nrow(targetConnections)) {
        
        #calculate current percent
        type1 <- gsub("^(..)...$", "\\1", targetConnections[i, "conn"])
        type2 <- gsub("^...(..)$", "\\1", targetConnections[i, "conn"])

        f <- current[current$Var1 == targetConnections[i, "conn"], "Freq"]
        types <- c(type1, type2)
        bool <- current$type1 %in% types | current$type2 %in% types
        s <- sum(current[bool, "Freq"])
        currentPercent <- (f / s) * 100
        
        #calculate number to add
        numerator <- (100 * f) - (s * targetConnections[i, "target"])
        denominator <- targetConnections[i, "target"] - 100
        add <- ceiling(numerator / denominator)
        mess <- paste(
            "Adjusting connections percentages by adding ",
            add,
            " multiplets.",
            sep = ""
        )
        print(mess)
        
        #synthesize and add connections to achieve target percentages
        for(y in 1:add) {
            tmp <- .makeMultuplet(
                2,
                types,
                multuplets,
                colnames(multuplets),
                .syntheticSinglets(ngenes, 1, cellTypes)[, types]
            )
            multuplets <- tmp[[1]]
            
            #remove self connections
            n <- ncol(multuplets)
            multuplets <- multuplets[, -c(n, (n - 2))]
            n <- ncol(multuplets)
            colnames(multuplets)[n] <- paste(types, collapse = "")
        }
        
        #check frequency
        after <- .quantifyConnections(colnames(multuplets))
        f <- after[after$Var1 == targetConnections[i, "conn"], "Freq"]
        bool <- after$type1 %in% types | after$type2 %in% types
        s <- sum(after[bool, "Freq"])
        currentPercent <- (f / s) * 100
        targetConnections[i, "current"] <- round(currentPercent)
    }
    
    if(!all.equal(targetConnections$target, targetConnections$current)) {
        stop("target connection percentage was not met. check code")
    }
    return(multuplets)
}

#Uses the multuplet naming convention to quantify the connections.
.quantifyConnections <- function(
    multupletNames
){
    split <- trimws(gsub("(.{2})", "\\1 ", multupletNames))
    suffixRemove <- gsub("^([A-Z0-9]* [A-Z0-9]*) \\..*$", "\\1", split)
    ss <- strsplit(suffixRemove, " ")
    l <- lapply(ss, function  (x) combn(x, 2))
    srt <- lapply(l, function(x) apply(x, 2, sort))
    p <- lapply(srt, function(x) paste(x[1, ], x[2, ], sep = "-"))
    d <- as.data.frame(table(unlist(p)), stringsAsFactors = FALSE)
    d$type1 <- gsub("^(..)...$", "\\1", d$Var1)
    d$type2 <- gsub("^...(..)$", "\\1", d$Var1)
    return(d)
}

#decide which connections should exist.
#This randomly assigns connections between all cell types and sets a target
#for the percentage of the total connections that the assigned connection should
#have.

#target: specifies the target percentage of connections for cell types that
#should be connected.
decideConnections <- function(
    cellTypes,
    target
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

.adjustSelf <- function(
    multuplets,
    singlets
){
    
    #count number of current connections
    conn <- .quantifyConnections(colnames(multuplets))
    currConn <- sum(conn$Freq)
    
    #calculate connections to add per cell type
    #currently dilutes the perfered connections by 3
    cellNames <- unique(c(conn$type1, conn$type2))
    perCellType <- (currConn * 2) / length(cellNames)
    mess <- paste(
        "Adding ",
        perCellType * length(cellNames),
        " self connections",
        sep = ""
    )
    print(mess)
    
    #add self connections
    #doublets = 1 connection
    #triplets = 3 connections
    #quadruplets = 6 connections
    connPerCellComb <- perCellType / 3
    dToAdd <- ceiling(connPerCellComb)
    tToAdd <- ceiling(connPerCellComb / 3)
    qToAdd <- ceiling(connPerCellComb / 6)
    
    for(i in cellNames) {
        #set valid self connections
        valid <- c(
            paste(rep(i, 2), collapse = ""),
            paste(rep(i, 3), collapse = ""),
            paste(rep(i, 4), collapse = "")
        )
        valid <- paste("^", valid, "$", sep="")
        
        #generate multuplets
        #doublets
        tmp <- data.frame()
        for(u in 1:dToAdd) {
            if(u == 1) {
                tmp <- .makeMultuplet(
                    2,
                    c(i, i),
                    tmp,
                    tmp,
                    singlets
                )
            } else {
                tmp <- .makeMultuplet(
                    2,
                    c(i, i),
                    tmp[[1]],
                    tmp[[2]],
                    singlets
                )
            }
        }
        
        bool <- grepl(valid[1], tmp[[2]])
        tmp[[1]] <- tmp[[1]][, bool]
        tmp[[2]] <- tmp[[2]][grepl(valid[1], tmp[[2]])]
        
        #triplets
        for(u in 1:tToAdd) {
            tmp <- .makeMultuplet(
                3,
                c(i, i),
                tmp[[1]],
                tmp[[2]],
                singlets
            )
        }
        
        bool1 <- grepl(valid[1], tmp[[2]])
        bool2 <- grepl(valid[2], tmp[[2]])
        bool <- bool1 | bool2
        tmp[[1]] <- tmp[[1]][, bool]
        tmp[[2]] <- tmp[[2]][bool]
        
        #quadruplets
        for(u in 1:qToAdd) {
            tmp <- .makeMultuplet(
                4,
                c(i, i),
                tmp[[1]],
                tmp[[2]],
                singlets
            )
        }
        
        bool1 <- grepl(valid[1], tmp[[2]])
        bool2 <- grepl(valid[2], tmp[[2]])
        bool3 <- grepl(valid[3], tmp[[2]])
        bool <- bool1 | bool2 | bool3
        tmp[[1]] <- tmp[[1]][, bool]
        tmp[[2]] <- tmp[[2]][bool]
        
        colnames(tmp[[1]]) <- tmp[[2]]
        
        #add
        multuplets <- cbind(multuplets, tmp[[1]])
        
    }
    
    return(multuplets)
}


################################################################################
#                                                                              #
# Plot distribution of connections per cell type.                              #
#                                                                              #
################################################################################

.plotConnectionDist <- function(
    multupletNames
){
    freq <- .quantifyConnections(multupletNames)
    cellTypes <- unique(c(freq$type1, freq$type2))
    percents <- sapply(cellTypes, function(x) {
        curr <- subset(freq, type1 == x | type2 == x)
        sum <- sum(curr$Freq)
        curr$percent <- 100 * (curr$Freq / sum)
    })
    rownames(percents) <- cellTypes
    
    p <- percents %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        rename(from = rowname) %>%
        gather(to, percent, -from) %>%
        add_column(connection = paste(.$from, .$to, sep="-")) %>%
        ggplot(., aes(connection, percent)) +
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







