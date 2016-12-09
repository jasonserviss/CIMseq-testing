#' cellTypesIncludedTest
#'
#' Subtitle
#'
#'
#' @name cellTypesIncludedTest
#' @rdname cellTypesIncludedTest
#' @aliases cellTypesIncludedTest
#' @param dataType Data type to use for test. Can be "synthetic" or "expression"
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords cellTypesIncludedTest
#' @examples
#'
#' #use demo data
#' cellTypesIncludedTest()
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq

cellTypesIncludedTest <- function(dataType, ...) {
    
    if(dataType == "synthetic") {
        data <- syntheticData
    } else {
        data <- expressionTestData
    }
    number <- seq(10, 4, -1)
    cellTypes <- unique(colnames(data[ ,grepl("s.", colnames(data))]))
    
    cellTypesInclcObj <- list()
    cellTypesIncluObj <- list()
    cellTypesInclsObj <- list()
    
    for(o in 1:length(number)) {
        
        idxSng <- .cellTypesIncSubsetSng(data, cellTypes, number[o])
        sng <- data[ ,idxSng]
        idxMul <- .cellTypesIncSubsetMul(data, sng, 6)
        mult <- data[ ,idxMul]
        counts <- cbind(sng, mult)
        
        cObj <- spCounts(counts, matrix(), "m.")
        uObj <- spUnsupervised(cObj, max=1000, max_iter=1000)
        sObj <- spSwarm(uObj, cores=6, limit="none", swarmsize=250, cutoff=0.01)
        
        cellTypesInclcObj[[o]] <- cObj
        cellTypesIncluObj[[o]] <- uObj
        cellTypesInclsObj[[o]] <- sObj
    }
    
    names(cellTypesInclcObj) <- as.character(number)
    names(cellTypesIncluObj) <- as.character(number)
    names(cellTypesInclsObj) <- as.character(number)
    
    if(dataType == "synthetic") {
        save(
            cellTypesInclcObj,
            cellTypesIncluObj,
            cellTypesInclsObj,
            file="data/cellTypesIncTestSyn.rda",
            compress="bzip2"
        )
    } else {
        save(
            cellTypesInclcObj,
            cellTypesIncluObj,
            cellTypesInclsObj,
            file="data/cellTypesIncTestExp.rda",
            compress="bzip2"
        )
    }
    
}

.cellTypesIncSubsetSng <- function(data, cellTypes, number) {
    idx <- which(colnames(data) %in% cellTypes[1:number])
    return(idx)
}

.cellTypesIncSubsetMul <- function(data, sng, n) {
    
    #get sng and mult types
    sngTypes <- gsub("s\\.(.*)", "\\1", unique(colnames(sng)))
    
    multTypes <- strsplit(
        gsub(
            "(.{2})",
            "\\1 ",
            gsub(
                "m\\.(.*)",
                "\\1",
                colnames(
                    data[ ,grepl("m.", colnames(data))]
                )
            )
        ), " "
    )
    
    #get indexes
    idx1 <- which(unlist(lapply(multTypes, function(r) all(r %in% sngTypes) == TRUE)))
    idx2 <- which(grepl("m.", colnames(data)) == TRUE )[idx1]
    
    #use only n doublets
    idx3 <- idx2[which(nchar(colnames(data[ ,idx2])) == 6)]
    
    if(n == "all") {
        return(idx3)
    } else {
       return(sample(idx3, n, replace=FALSE))
    }
    
}

#' plotCellTypesIncludedTest
#'
#' Subtitle
#'
#'
#' @name plotCellTypesIncludedTest
#' @rdname plotCellTypesIncludedTest
#' @aliases plotCellTypesIncludedTest
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords plotCellTypesIncludedTest
#' @examples
#'
#' #use demo data
#' plotCellTypesIncludedTest()
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq
#' @import ggplot2
#' @importFrom plyr rbind.fill
#' @importFrom reshape2 melt

plotCellTypesIncludedTest <- function(...) {
    
    l <- lapply(cellTypesInclsObj, function(t) getData(t, "codedSwarm"))
    res <- lapply(l, function(d) calculateConnections(d))
    df <- rbind.fill(res)
    df$cellTypesInc <- names(res)
    df$cellTypesInc <- sort(df$cellTypesInc)
    
    p <- .pCellTypesIncludedTest(df)
    p
    return(p)
}


.pCellTypesIncludedTest <- function(data) {
    p <- ggplot(data, aes(x=cellTypesInc, y=percentFoundConnections))+
        geom_boxplot()
    return(p)
}








