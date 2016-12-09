#context("generalFunctions")

##run test .swarmNetworkDF
test_that("check that the .swarmNetworkDF function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    
    res <- as.data.frame(matrix(
        c(
            0,1,1,0,
            1,1,0,0,
            0,0,1,1,
            1,0,0,1
        ),
        nrow=4,
        byrow=TRUE
    ))
    
    colnames(res) <- LETTERS[1:4]
    
    input <- data.frame(
        fopt = "holder",
        sampleName = as.character(1:4),
        as.data.frame(res)
    )
    
    #setup expected data
    expected <- data.frame(
        from = c("B", "A", "C", "A"),
        to = c("C", "B", "D", "D"),
        multuplet = as.character(1:4),
        connType = rep("detected", 4),
        stringsAsFactors = FALSE
    )
    
    #run function
    output <- .swarmNetworkDF(input)
    
    #test
    expect_equivalent(expected, output)
    
})

##run test .processSampleNames
test_that("check that the .processSampleNames function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    a <- paste("m.", LETTERS[1:4], 1, sep="")
    b <- paste(rev(LETTERS)[1:4], 1, sep="")
    input <- data.frame(sampleNames = paste(a, b, sep=""))
    
    #setup expected data
    expected <- data.frame(
        from = c("A1", "B1", "C1", "D1"),
        to = c("Z1", "Y1", "X1", "W1"),
        multuplet = paste(a, b, sep=""),
        connType = rep("real", 4),
        stringsAsFactors = FALSE
    )
    
    #run function
    output <- .processSampleNames(input)
    
    #test
    expect_equivalent(expected, output)
    
})

##run test .percentFoundConnections
test_that("check that the .percentFoundConnections function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    inputDC <- data.frame(from = LETTERS[1:10], to=rev(LETTERS)[1:10])[1:9,]
    inputRC <- data.frame(from = LETTERS[1:10], to=rev(LETTERS)[1:10])
    
    #setup expected data
    expected <- 90
    
    #run function
    output <- .percentFoundConnections(inputDC, inputRC)
    
    #test
    expect_equivalent(expected, output)
    
})

##run test calculateConnections
test_that("check that the calculateConnections function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    res <- as.data.frame(matrix(
        c(
            0,1,1,0,
            1,1,0,0,
            0,0,1,1,
            1,0,0,1
        ),
        nrow=4,
        byrow=TRUE
    ))
    
    colnames(res) <- c("A1", "B1", "C1", "D1")
    
    input <- data.frame(
        fopt = "holder",
        sampleName = paste("m.", c("B1C1", "A1C1", "A1B1", "A1D1"), sep=""),
        as.data.frame(res),
        stringsAsFactors=FALSE
    )
    
    #setup expected data
    expected <- data.frame(
        multuplet = sort(paste("m.", c("B1C1", "A1C1", "A1B1", "A1D1"), sep="")),
        percentFoundConnections = c(0,0,100,100),
        stringsAsFactors=FALSE
    )
    
    #run function
    output <- calculateConnections(input)[[1]]
    
    #test
    expect_equivalent(expected, output)
    
    
    ###TEST2####
    #prepare normal input data
    input <- as.data.frame(matrix(rep(NA,4), ncol=4))
    colnames(input) <- c("m.A1B1", "m.A1C1", "m.A1D1", "m.A1E1")
    
    #setup expected data
    expected <- table(data.frame(
        from=rep("A1", 4),
        to=c("B1", "C1", "D1", "E1")
    ))
    
    #run function
    output <- calculateConnections(input, type="multuplets")
    
    #test
    expect_equivalent(expected, output)

})

