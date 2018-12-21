#' groupSizeTest
#'
#' Subtitle
#'
#'
#' @name groupSizeTest
#' @rdname groupSizeTest
#' @aliases groupSizeTest
#' @param dataType Data type to use for test. Can be "synthetic" or "expression"
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords groupSizeTest
#' @examples
#'
#' #use demo data
#' groupSizeTest()
#'
#'
#'
NULL
#' @export
#' @import CIMseq

groupSizeTest <- function(dataType, ...) {
    
    if(dataType == "synthetic") {
        data <- syntheticData
        multNames <- c("m.A1B1", "m.I1J1", "m.A1B1C1", "m.H1I1J1", "m.A1B1C1D1", "m.G1H1I1J1")
        s <- c(seq(100, 10, -10), 5, 2)
    } else {
        #note that the expression data isn't really working due to the fact that the H1 group in the uObj only coonsists of 5 cells. If the analysis starts using 5 cells per group, it already is unable to correctly deconvolute the multuplets on the first iteration.
        data <- expressionTestData
        multNames <- c("m.A1B1", "m.F1H1", "m.A1B1C1", "m.F1G1H1", "m.A1B1C1D1", "m.E1F1G1H1")
        s <- seq(5, 2, -1)
    }
    
    sngNames <- paste("s.", LETTERS[1:10], 1, sep="")
    
    sng <- data[,grepl("s.", colnames(data))]
    mult <- data[ ,colnames(data) %in% multNames]
    
    groupSizeTestcObj <- list()
    groupSizeTestuObj <- list()
    groupSizeTestsObj <- list()
    

    for(i in 1:length(s)) {
        if(i==1) {
            testCounts <- cbind(sng, mult)
        } else {
            testCounts <- .makeCounts(sngNames, sng, mult, i, s)
        }
        #run optimization
        optOutput <- .runOptimization(testCounts, groupSizeTestcObj, groupSizeTestuObj, groupSizeTestsObj, i)
        groupSizeTestcObj <- optOutput[[1]]
        groupSizeTestuObj <- optOutput[[2]]
        groupSizeTestsObj <- optOutput[[3]]
    }
    
    names(groupSizeTestcObj) <- as.character(s)
    names(groupSizeTestuObj) <- as.character(s)
    names(groupSizeTestsObj) <- as.character(s)
    
    if(dataType == "synthetic") {
        
        groupSizeTestcObjSyn <- groupSizeTestcObj
        groupSizeTestuObjSyn <- groupSizeTestuObj
        groupSizeTestsObjSyn <- groupSizeTestsObj
        
        save(
            groupSizeTestcObjSyn,
            groupSizeTestuObjSyn,
            groupSizeTestsObjSyn,
            file="data/groupSizeTestSyn.rda",
            compress="bzip2"
        )
        return(groupSizeTestsObj)
    } else {
        
        groupSizeTestcObjExp <- groupSizeTestcObj
        groupSizeTestuObjExp <- groupSizeTestuObj
        groupSizeTestsObjExp <- groupSizeTestsObj
        
        save(
            groupSizeTestcObjExp,
            groupSizeTestuObjExp,
            groupSizeTestsObjExp,
            file="data/groupSizeTestExp.rda",
            compress="bzip2"
        )
        return(groupSizeTestsObj)
    }
}

.makeCounts <- function(sngNames, sng, mult, i, s) {
    for(o in 1:length(sngNames)) {
        currName <- sngNames[o]
        currData <- sng[ ,grepl(currName, colnames(sng))][ ,1:s[i]]
        
        if(o==1) {
            testCounts <- currData
        } else {
            testCounts <- cbind(testCounts, currData)
        }
    }
    testCounts <- cbind(testCounts, mult)
    return(testCounts)
}

.runOptimization <- function(testCounts, groupSizeTestcObj, groupSizeTestuObj, groupSizeTestsObj, i) {
    cObj <- spCounts(testCounts, matrix(), "m.")
    
    if(i == 12) {
        p <- 6
    } else {
        p <- 10
    }
    
    uObj <- spUnsupervised(cObj, max=1000, max_iter=1000, perplexity=p)
    sObj <- spSwarm(uObj, limit="none", maxiter=10, swarmsize=250, cores=6)
    
    groupSizeTestcObj[[i]] <- cObj
    groupSizeTestuObj[[i]] <- uObj
    groupSizeTestsObj[[i]] <- sObj
    
    return(list(groupSizeTestcObj, groupSizeTestuObj, groupSizeTestsObj))
}

#' plotGroupSizeTest
#'
#' Subtitle
#'
#'
#' @name plotGroupSizeTest
#' @rdname plotGroupSizeTest
#' @aliases plotGroupSizeTest
#' @param ... additional arguments to pass on
#' @return Saves results in .rda file.
#' @author Jason T. Serviss
#' @keywords plotGroupSizeTest
#' @examples
#'
#' #use demo data
#' plotGroupSizeTest()
#'
#'
#'
NULL
#' @export
#' @importFrom plyr rbind.fill
#' @import ggplot2
#' @importFrom ggthemes theme_few scale_fill_economist

plotGroupSizeTest <- function(dataType = "synthetic", ...) {
    
    if(dataType == "synthetic") {
        sObj <- groupSizeTestsObjSyn
    } else {
        sObj <- groupSizeTestsObjExp
    }
    
    plot <- .processGroupSizeTest(sObj)
    p <- .plotGroupSizeProcessedRes(plot)
    p
    return(p)
}

.plotGroupSizeProcessedRes <- function(data) {
    p <- ggplot(data, aes(x=nrPerGrp, y=percentFoundConnections, fill=cells))+
    geom_boxplot()+
    theme_few()+
    scale_fill_economist()+
    labs(
        x="Points per cell type",
        y="% connections detected",
        title="Group Size Test"
        )+
        theme(
            legend.position="top",
            legend.title=element_text(size=15),
            legend.text=element_text(size=15),
            axis.title=element_text(size=17),
            axis.text=element_text(size=15),
            plot.title=element_text(
                hjust=0.5,
                family="Arial",
                face="bold",
                size=24,
                margin=margin(b=15)
            )
        )+
        guides(
            fill=guide_legend(
                title="Cells in multuplet",
                title.position = "top",
                title.hjust =0.5
            )
        )

    return(p)
}

.processGroupSizeTest <- function(sObj) {
    r <- lapply(
        sObj,
        function(i)
            calculateConnections(
                getData(i, "codedSwarm")
            )[[1]]
    )
    
    plot <- rbind.fill(r)
    
    plot$nrPerGrp <- sort(
        as.numeric(
            rep(
                names(r),
                nrow(r[[1]])
            )
        ),
        decreasing=TRUE
    )
    
    plot$multuplet <- as.character(plot$multuplet)
    plot$cells <- ifelse(nchar(plot$multuplet) == 6, "two",
        ifelse(nchar(plot$multuplet) == 8, "three", "four")
    )
    
    plot$nrPerGrp <- factor(plot$nrPerGrp, levels=c(seq(100, 10, -10), 5, 2))
    return(plot)
}