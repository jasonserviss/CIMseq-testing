
#' noiseDataTest
#'
#' Subtitle
#'
#'
#' @name noiseDataTest
#' @rdname noiseDataTest
#' @aliases noiseDataTest
#' @param dataType Data type to use for test. Can be "synthetic" or "expression"
#' @param cores Number of cores to run spSwarm on.
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords noiseDataTest
#' @examples
#'
#' #use demo data
#' noiseDataTest()
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq
#' @import foreach
#' @import doMC

noiseDataTest <- function(innerCores=1, outerCores=10, outPath='data', n=10, dataType, ...) {
    
    #load test data
    tmp <- .generateNoiseData(n, dataType)
    noise <- tmp[[1]]
    
    #extract data info
    names <- tmp[[2]]
    noisePercent <- dimnames(noise)[[3]]
    multCellNr <- dimnames(noise)[[4]]
    
    #load spUnsupervised for dataset
    if(dataType == "synthetic") {
        data(syntheticDataTest)
        uObj <- syntheticDataUnsupervised
        rm(list=c(
            "syntheticDataTable",
            "syntheticDataTest",
            "syntheticDataSwarm",
            "syntheticDataUnsupervised"
        ))
    } else {
        data(expressionTest)
        uObj <- expressionTestUnsupervised
    }
    
    if(outerCores > 1) {registerDoMC(outerCores)}
    results <- .outer(noise, names, uObj, innerCores)
    if(outerCores > 1) {registerDoSEQ()}
    
    #add names
    names(results) <- multCellNr
    names(results[[1]]) <- noisePercent
    names(results[[2]]) <- noisePercent
    names(results[[3]]) <- noisePercent
    
    #save
    if(dataType == "synthetic") {
        noiseDataSyn <- tmp
        noiseDataSwarmSyn <- results
        noiseDataUnsupervisedSyn <- uObj
        
        save(noiseDataSyn, file=paste(outPath, 'noiseDataSyn.rda', sep='/'), compress="bzip2")
        save(noiseDataUnsupervisedSyn, file=paste(outPath, 'noiseDataUnsupervisedSyn.rda', sep='/'), compress="bzip2")
        save(noiseDataSwarmSyn, file=paste(outPath, 'noiseDataSwarmSyn.rda', sep='/'), compress="bzip2")
        return(noiseDataSwarmSyn)
    } else {
        noiseDataExp <- tmp
        noiseDataSwarmExp <- results
        noiseDataUnsupervisedExp <- uObj
        
        save(noiseDataExp, file=paste(outPath, 'noiseDataExp.rda', sep='/'), compress="bzip2")
        save(noiseDataUnsupervisedExp, file=paste(outPath, 'noiseDataUnsupervisedExp.rda', sep='/'), compress="bzip2")
        save(noiseDataSwarmExp, file=paste(outPath, 'noiseDataSwarmExp.rda', sep='/'), compress="bzip2")
        return(noiseDataSwarmExp)
    }
}

#iterates over cell Number dimension in noise variable
.outer <- function(noise, names, uObj, innerCores) {
    results <- list()
    for(i in 1:dim(noise)[4]) {
        currMultNumb <- noise[,,,i]
        currNames <- names[,,,i]
        inner <- .inner(currMultNumb, currNames, uObj, innerCores)
        results[[i]] <- inner
    }
    return(results)
}

#iterates over amount of noise dimension in noise variable
.inner <- function(currMultNumb, currNames, uObj, innerCores) {
    inner <- foreach(u = 1:dim(currMultNumb)[3]) %dopar% {
        currMatrix <- currMultNumb[,,u]
        colnames(currMatrix) <- currNames[,u]
        
        #incorporate multuplet data into uObj
        sampleType <- getData(uObj, "sampleType")
        counts <- getData(uObj, "counts")
        counts <- counts[ ,sampleType=="Singlet"]
        sampleType <- sampleType[sampleType=="Singlet"]
        
        counts(uObj) <- cbind(counts, currMatrix)
        sampleType(uObj) <- c(sampleType, rep("Multuplet", ncol(currMatrix)))
        counts.log(uObj) <- sp.scRNAseq:::.norm.log.counts(cbind(counts, currMatrix))
        
        #run pySwarm
        noiseDataSwarm <- spSwarm(uObj, limit="none", maxiter=10, swarmsize=250, cores=innerCores)
    }
    return(inner)
}

#master function for noise data generation
.generateNoiseData <- function(n, dataType) {
    set.seed(11)
    
    if(dataType == "synthetic") {
        data <- syntheticData
    } else {
        data <- expressionTestData
    }
    
    #subset multuplets
    multuplets <- data[ ,grep("m.", colnames(data))]
    
    #declare array to hold noise data
    noise <- array(NA, dim=c(nrow(multuplets), n, 11, 3))
    dimnames(noise)[[3]] <- as.character(seq(0.00,1,0.1))
    dimnames(noise)[[4]] <- c("two", "three", "four")
    
    percent <- seq(0.1,1,0.1)
    
    #pick multuplets to be used and store the names in the names array
    names <- array(NA, dim=c(1, n, 11, 3))
    names[,,,1] <- .pickMultuplets(multuplets, 'two', n)
    names[,,,2] <- .pickMultuplets(multuplets, 'three', n)
    names[,,,3] <- .pickMultuplets(multuplets, 'four', n)
    
    #add multuplets with no noise according to the names in the names array
    noise[,,1, 'two'] <- data[ ,names[,,,1][,1]]
    noise[,,1, 'three'] <- data[ ,names[,,,2][,1]]
    noise[,,1, 'four'] <- data[ ,names[,,,3][,1]]
    
    #inject noise
    noise <- .noise(noise, percent, 'two')
    noise <- .noise(noise, percent, 'three')
    noise <- .noise(noise, percent, 'four')
    
    noiseData <- list(noise, names)
    return(noiseData)
}

#handles picking the multuplets and returning their names
.pickMultuplets <- function(multuplets, nCells, nMultuplets) {
    nChar <- c(
        two=6,
        three=8,
        four=10
    )
    
    allTheseMult <- colnames(multuplets[ ,nchar(colnames(multuplets)) == nChar[nCells]])
    justTheseMult <- allTheseMult[sample(1:length(allTheseMult), nMultuplets, replace=FALSE)]
    return(justTheseMult)
}

#handles generating the noise and adding it to the noise array
.noise <- function(noise, percent, nCells) {
    justTheseMult <- noise[,,1, nCells]
    
    for( i in 1:length(percent) ) {
        tmp <- justTheseMult
        
        #pick rows to be substituted
        rows <- .subsetRows(justTheseMult, i, percent)
        picked <- rows[[1]]
        sub <- rows[[2]]
        
        #add noise
        tmp <- .addNoise(sub, tmp, picked)
        
        #insert in array
        noise[,,i+1, nCells] <- tmp
    }
    
    return(noise)
}

#subsets the rows where noise should be added dependant on the percentage of desired noise
.subsetRows <- function(justTheseMult, i, percent) {
    picked <- sample(
        1:nrow(justTheseMult),
        nrow(justTheseMult) * percent[i],
        replace = FALSE
    )
    sub <- justTheseMult[picked, ]
    return(list(picked, sub))
}

#adds noise to the previously subsetted rows and
#inserts it back into the tmp variable for subsequent
#insertion into the noise arrray
.addNoise <- function(sub, tmp, picked) {
    #reordered <- sub[sample.int(nrow(sub)),]
    reordered <- matrix(
        sample(
            seq(
                min(sub),
                max(sub),
                0.0005
            ),
            length(sub)
        ),
        ncol=ncol(sub)
    )
    
    #insert
    tmp[picked, ] <- reordered
    return(tmp)
}

#' noiseDataTestPlot
#'
#' Subtitle
#'
#'
#' @name noiseDataTestPlot
#' @rdname noiseDataTestPlot
#' @aliases noiseDataTestPlot
#' @param cores Number of cores to run spSwarm on.
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords noiseDataTestPlot
#' @examples
#'
#' #use demo data
#' noiseDataTestPlot()
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq
#' @importFrom plyr rbind.fill ddply
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_fill_economist
#' @importFrom dplyr combine

noiseDataTestPlot <- function(
    dataType="synthetic",
    edge.cutoff = 1/10.5
){
    
    if(dataType == "synthetic") {
        data <- noiseDataSwarmSyn
    } else {
        data <- noiseDataSwarmExp
    }
    
    res <- lapply(
        data,
        function(e)
        lapply(
            e,
            function(p)
                calculateConnections(getData(p, "spSwarm"))[[1]]
        )
    )
    
    pc <- rbind.fill(combine(res))
    pc$cells <- sort(rep(names(res), nrow(pc)/3), decreasing=TRUE)
    pc$percentNoise <- sort(rep(names(res[[1]]), 10))
    
    ggplot(pc, aes(x=percentNoise, y=percentFoundConnections, fill=cells))+
        geom_boxplot()+
        theme_few()+
        scale_fill_economist()
    
}


