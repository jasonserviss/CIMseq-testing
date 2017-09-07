
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
    n = 100,
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
    idx <- sample(1:ncol(multuplets), n, replace = FALSE)
    testMultuplets <- multuplets[, idx]
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
            nrow = ngenes
        )
        if( i == 1 ) {
            singlets <- counts
        } else {
            singlets <- cbind(singlets, counts)
        }
    }
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
    
    multuplets <- .adjustMultuplets(singlets, multuplets, ngenes, cellTypes)
    
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
        
        set.seed(2918834 + u)
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

.adjustMultuplets <- function(singlets, multuplets, ngenes, cellTypes) {
    
    #decide which connection preferences should exist
    targetConnections <- decideConnections(unique(colnames(singlets)), 30)
    
    #quantify current connections
    current <- .quantifyConnections(colnames(multuplets))
    
    #adjust connection frequency
    
    for(i in 1:nrow(targetConnections)) {
        
        #calculate current percent
        type1 <- gsub("^(..)...$", "\\1", targetConnections[i, "conn"])
        type2 <- gsub("^...(..)$", "\\1", targetConnections[i, "conn"])

        f <- current[current$Var1 == targetConnections[i, "conn"], "Freq"]
        s <- sum(current[current$type1 %in% c(type1, type2) | current$type2 %in% c(type1, type2), "Freq"])
        currentPercent <- (f / s) * 100
        
        #calculate number to add
        numerator <- (100 * f) - (s * targetConnections[i, "target"])
        denominator <- targetConnections[i, "target"] - 100
        add <- numerator / denominator
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
                c(type1, type2),
                multuplets,
                colnames(multuplets),
                .syntheticSinglets(ngenes, 1, cellTypes)[,c(type1, type2)]
            )
            multuplets <- tmp[[1]]
            
            #remove self connections
            n <- ncol(multuplets)
            multuplets <- multuplets[, -c(n, (n - 2))]
            colnames(multuplets)[ncol(multuplets)] <- paste(type1, type2, sep = "")
        }
        
        #check frequency
        after <- .quantifyConnections(colnames(multuplets))
        f <- after[after$Var1 == targetConnections[i, "conn"], "Freq"]
        s <- sum(after[after$type1 %in% c(type1, type2) | after$type2 %in% c(type1, type2), "Freq"])
        currentPercent <- (f / s) * 100
        targetConnections[i, "current"] <- currentPercent
    }
    
    if(!all.equal(targetConnections$target, targetConnections$current)) {
        stop("target connection percentage was not met. check code")
    }
    return(multuplets)
}

#Uses the multuplet naming convention to quantify the connections.
.quantifyConnections <- function(multupletNames) {
    split <- trimws(gsub("(.{2})", "\\1 ", multupletNames))
    suffixRemove <- gsub("^([A-Z0-9]* [A-Z0-9]*) \\..*$", "\\1", split)
    ss <- strsplit(suffixRemove, " ")
    l <- lapply(ss, function  (x) combn(x, 2))
    srt <- lapply(l, function(x) apply(x, 2, sort))
    p <- lapply(srt, function(x) paste(x[1, ], x[2, ], sep = "-"))
    d <- as.data.frame(table(unlist(p)), stringsAsFactors = FALSE)
    d$type1 <- gsub("^(..)...$", "\\1", d$Var1)
    d$type2 <- gsub("^...(..)$", "\\1", d$Var1)
    d
}

#decide which connections should exist
decideConnections <- function(cellTypes, target) {
    
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

################################################################################
#                                                                              #
# Plot distribution of connections per cell type.                              #
#                                                                              #
################################################################################

.plotConnectionDist <- function(multupletNames) {
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
        geom_bar(aes(fill = to), stat = "identity", position = position_dodge(width = 1)) +
        facet_grid(from~to, scales = "free", space = "free") +
        ggthemes::scale_fill_ptol() +
        ggthemes::theme_few() +
        theme(
            axis.text.x = element_text(angle = 90)
        ) +
        guides(fill = FALSE)
}








