
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
            nrow = ngenes
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
    
    #multuplets <- .adjustMultuplets(singlets, multuplets)
    
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


.adjustMultuplets <- function(singlets, multuplets) {
    
    combos <- combn(unique(colnames(singlets)), 2)
    
    #find connections in multuplets
    register <- sapply(
        1:ncol(combos),
        function(x)
            grepl(combos[1, x], colnames(multuplets)) & grepl(combos[2, x], colnames(multuplets))
    )
    rownames(register) <- colnames(multuplets)
    colnames(register) <- paste(combos[1, ], combos[2,], sep = "-")
    
    #identify multuplets with 1 and >1 connections
    one <- which(rowSums(register) == 1)
    gtOne <- which(rowSums(register) > 1)
    
    #decide which connection preferences should exist
    targetConnections <- decideConnections(unique(colnames(singlets)))
    
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
        
        #synthesize and add
        #the problem with 100 cells each cell type is that it is not enough to
        #meet the target so either synthesize more cells or drop the target
        for(y in 1:add) {
            tmp <- .makeMultuplet(
                2,
                c(type1, type2),
                multuplets,
                colnames(multuplets),
                singlets
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
decideConnections <- function(cellTypes) {
    
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
        target = 50,
        current = 0,
        stringsAsFactors = FALSE
    ))
}


.optimizeFun <- function(x) {
    x <- round(x)
    for(i in 1:length(x)) {
        #subset multiplets that include this connection
    }
}



#ratio: the ratio of cell types to have fixed connection frequencies
.pickFixedFreq <- function(singlets, multuplets, cellTypes, ratio = 0.5) {
    nFixed <- round(cellTypes * ratio)
    
    #pick cell types to fix
    set.seed(22389)
    fixed <- sample(unique(colnames(singlets)), size = nFixed, replace = FALSE)
    combos <- combn(fixed, 2) #add self connections?
    
    #find multuplets with fixed connections
    register <- sapply(
        1:ncol(combos),
        function(x)
            grepl(combos[1, x], colnames(multuplets)) & grepl(combos[2, x], colnames(multuplets))
    )
    rownames(register) <- colnames(multuplets)
    colnames(register) <- paste(combos[1, ], combos[2,], sep = "-")
    
    #identify multuplets with 1 and >1 fixed connections
    one <- which(rowSums(register) == 1)
    gtOne <- which(rowSums(register) > 1)
    
    #quantify connections per fixed connection type in multiplets with >1
    #fixed connections
    freq <- .quantifyConnections(rownames(register)[gtOne], fixed)
    
    #set up target frequencies for existing connections
    #each cell type should have one other cell type that it interacts with 30%
    #of the time
    targetFreq <- .targetFreq(freq)
    
    #compare target frequencies to quantified connections and adjust
    
}

.targetFreq <- function(freq) {

    freq$target <- NA
    freq$self <- ifelse(freq$type1 == freq$type2, TRUE, FALSE)
    done <- c()
    
    for(i in fixed) {
        if(!all(is.na(freq[freq$type1 == i | freq$type2 == i, 4]))) {
            next
        }
        
        idxs <- which(
            (freq$type1 == i | freq$type2 == i ) &
            freq$self == FALSE &
            !freq$type1 %in% done &
            !freq$type2 %in% done
        )
        
        set.seed(23089)
        idx <- sample(idxs, size = 1)
        freq[idx, "target"] <- 50
        done <- c(done, i)
    }
    
    return(conTypes[is.na(conTypes$target) == FALSE, ])
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
