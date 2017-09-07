
#' calculateConnections
#'
#' Subtitle
#'
#'
#' @name calculateConnections
#' @rdname calculateConnections
#' @aliases calculateConnections
#' @param codedSwarm Output in spSwarm object.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss
#' @keywords calculateConnections
#' @examples
#'
#' #use demo data
#' calculateConnections()
#'
#'
#'
NULL
#' @export
#' @import sp.scRNAseq
#' @importFrom plyr ddply

calculateConnections <- function(
    swarm,
    type = "codedSwarm",
    edge.cutoff = 0
){
    
    if(type == "multuplets") {
        conns <- .processSampleNames(colnames(swarm))
        tab <- table(conns[ ,c("from", "to")])
        return(tab)
    }
    
    #detectedConnections <- .swarmNetworkDF(swarm)
    detectedConnections <- .swarmNetworkDF(swarm, edge.cutoff)
    realConnections <- .processSampleNames(rownames(swarm))
    tab <- table(
        rbind(
            realConnections[,c("from", "to", "connType")],
            detectedConnections[,c("from", "to", "connType")]
        )
    )
    
    data <- rbind(realConnections, detectedConnections)
    
    pfc <- ddply(data, c("multuplet"), function(x)
        .percentFoundConnections(
            x[x$connType=="detected", ],
            x[x$connType=="real", ]
        )
    )
    colnames(pfc) <- c("multuplet", "percentFoundConnections")
    return(list(pfc, tab))
}

.percentFoundConnections <- function(dC, rC) {
    detected <- paste(dC$from, dC$to, sep="-")
    real <- paste(rC$from, rC$to, sep="-")
    found <-  real %in% detected
    
    nrFound <- length(which(found == TRUE))
    total <- length(found)

    if(nrFound == 0) {
        percent <- 0
    } else {
        percent <- nrFound / total
    }
        
    return(percent*100)
}

.processSampleNames <- function(data) {
    realConn1 <- gsub("m.([A-Z]1[A-Z]1)", "\\1", data)
    realConn2 <- trimws(gsub("(.{2})", "\\1 ", realConn1), which = "both")
    suffixRemove <- gsub("^([A-Z0-9]* [A-Z0-9]*) \\..*$", "\\1", realConn2)
    realConn3 <- strsplit(suffixRemove, " ")
    realConn4 <- lapply(realConn3, function(x) combn(x, 2))
    realConn5 <- lapply(realConn4, function(o) apply(o, 2, sort))
    nConnections <- lapply(realConn5, function(y) ncol(y))
    
    from <- unlist(lapply(realConn5, function(u) u[1,]))
    to <- unlist(lapply(realConn5, function(u) u[2,]))
    
    output <- data.frame(
        from = from,
        to = to,
        multuplet = as.character(
            rep(
                data,
                unlist(nConnections)
            )
        ),
        connType="real",
        stringsAsFactors=FALSE
    )
    return(output)
}

.swarmNetworkDF <- function(x, edge.cutoff) {
    for(o in 1:nrow(x)) {
        ind <- which(x[o,] > edge.cutoff)
        
        if(length(ind) == 0) {
            combs <- data.frame(V1=NA, V2=NA)
        } else if(length(ind) == 1) {
            combs <-  data.frame(V1=names(x)[ind], V2=names(x)[ind], stringsAsFactors=FALSE)
        } else {
            combs <- as.data.frame(t(combn(names(x)[ind],2)), stringsAsFactors=FALSE)
        }
        
        
        combs$multuplet <- rep(rownames(x)[o], nrow(combs))
        
        if( o == 1 ) {
            connections <- combs
        } else {
            connections <- rbind(connections, combs)
        }
    }
    
    connections$connType <- "detected"
    colnames(connections) <- c("from", "to", "multuplet", "connType")
    return(connections)
}

