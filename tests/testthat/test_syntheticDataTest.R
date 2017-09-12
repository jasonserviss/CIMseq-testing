#context("syntheticDataTest")

################################################################################
#                                                                              #
# Singlets.                                                                    #
#                                                                              #
################################################################################

##run test syntheticSinglets
test_that("check syntheticSinglets function outputs the expected result", {
    
    ###TEST1####
    #setup input
    nGenes <- 100
    nCells <- 10
    nCellTypes <- 3
    singletExpansion <- 2
    
    #setup expected data
    expectedDim <- c(100, 30)
    expectedDimFull <- c(100, 1500)
    
    prefix <- sort(rep(LETTERS[1:nCellTypes], (nCells * singletExpansion)))
    suffix <- 1:(nCells * singletExpansion)
    expectedNames <- paste(prefix, suffix, sep = "")
    expectedClass <- "data.frame"
    
    #run function
    output <- .syntheticSinglets(nGenes, nCells, nCellTypes, singletExpansion)
    outputSingletsFull <- output[[1]]
    outputSinglets <- output[[2]]
    singletNames <- c(colnames(outputSingletsFull), colnames(outputSinglets))
    
    #test
    expect_equivalent(dim(outputSinglets), expectedDim)
    expect_equivalent(dim(outputSingletsFull), expectedDim)

    expect_true(all(singletNames %in% expectedNames))
    
    expect_equivalent(class(outputSinglets), expectedClass)
    expect_equivalent(class(outputSingletsFull), expectedClass)
})

################################################################################
#                                                                              #
# Multiplets units.                                                            #
#                                                                              #
################################################################################

##run test makeMultuplet
test_that("check makeMultuplet with doublets", {
    
    #set up input data
    ##arguments to syntheticSinglets
    cellTypes <- paste(LETTERS[1:3], 1, sep = "")
    nGenes <- 100
    nCells <- 3
    nCellTypes <- length(cellTypes)
    singletExpansion <- 2
    

    ##arguments to makeMultuplet
    nCellsInMultiplet <- 2
    singlets <- .syntheticSinglets(
        nGenes,
        nCells,
        nCellTypes,
        singletExpansion
    )[[2]]
    colnames(singlets) <- sort(rep(cellTypes, 3))
    multuplets <- data.frame(row.names = 1:nrow(singlets))

    #setup expected data
    expectedDim <- c(nGenes, sum(length(cellTypes), ncol(combn(cellTypes, 2))))
    expectedMultuplets <- c("A1A1", "A1B1", "A1C1", "B1B1", "B1C1", "C1C1")
    
    #run function
    output <- makeMultuplet(
        nCellsInMultiplet,
        cellTypes,
        multuplets,
        singlets
    )
    
    #test
    expect_equivalent(dim(output), expectedDim)
    expect_equivalent(colnames(output), expectedMultuplets)
})

test_that("check makeMultuplet with triplets", {
    
    #set up input data
    ##arguments to syntheticSinglets
    cellTypes <- paste(LETTERS[1:3], 1, sep = "")
    nGenes <- 100
    nCells <- 3
    nCellTypes <- length(cellTypes)
    singletExpansion <- 2
    
    ##arguments to makeMultuplet
    nCellsInMultiplet <- 3
    singlets <- syntheticSinglets(
        nGenes,
        nCells,
        nCellTypes,
        singletExpansion
    )[[2]]
    colnames(singlets) <- sort(rep(cellTypes, 3))
    multuplets <- data.frame(row.names = 1:nrow(singlets))

    #setup expected data
    expectedDim <- c(nGenes, sum(length(cellTypes), 7))
    expectedMultuplets <- c(
        "A1A1A1", "A1A1B1", "A1A1C1", "A1B1B1", "A1B1C1",
        "A1C1C1", "B1B1B1", "B1B1C1", "B1C1C1", "C1C1C1"
    )
    
    #run function
    output <- makeMultuplet(
        nCellsInMultiplet,
        cellTypes,
        multuplets,
        singlets
    )
    
    #test
    expect_equivalent(dim(output), expectedDim)
    expect_equivalent(colnames(output), expectedMultuplets)
})

test_that("check makeMultuplet with quadruplets", {
    
    #set up input data
    ##arguments to syntheticSinglets
    cellTypes <- paste(LETTERS[1:3], 1, sep = "")
    nGenes <- 100
    nCells <- 3
    nCellTypes <- length(cellTypes)
    singletExpansion <- 2
    
    
    ##arguments to makeMultuplet
    nCellsInMultiplet <- 4
    singlets <- syntheticSinglets(
        nGenes,
        nCells,
        nCellTypes,
        singletExpansion
    )[[2]]
    colnames(singlets) <- sort(rep(cellTypes, 3))
    multuplets <- data.frame(row.names = 1:nrow(singlets))

    #setup expected data
    expectedDim <- c(nGenes, sum(length(cellTypes), 12))
    expectedMultuplets <- c(
        "A1A1A1A1", "A1A1A1B1", "A1A1A1C1", "A1A1B1B1", "A1A1B1C1",
        "A1A1C1C1", "A1B1B1B1", "A1B1B1C1", "A1B1C1C1", "A1C1C1C1",
        "B1B1B1B1", "B1B1B1C1", "B1B1C1C1", "B1C1C1C1", "C1C1C1C1"
    )
    
    #run function
    output <- makeMultuplet(
        nCellsInMultiplet,
        cellTypes,
        multuplets,
        singlets
    )
    
    #test
    expect_equivalent(dim(output), expectedDim)
    expect_equivalent(colnames(output), expectedMultuplets)
    
})

################################################################################
#                                                                              #
# Adjust multiplets.                                                           #
#                                                                              #
################################################################################

##run test decideConnections
test_that("check decideConnections function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    cellTypes <- paste(LETTERS[1:4], 1, sep = "")
    target <- 20
    
    #setup expected data
    expectedConn <- c("A1-D1", "B1-C1")
    expectedTarget <- rep(target, 2)
    expectedCurrent <- rep(0, 2)
    
    #run function
    output <- decideConnections(cellTypes, target)
    
    #test
    expect_equivalent(expectedConn, output$conn)
    expect_equivalent(expectedTarget, output$target)
    expect_equivalent(expectedCurrent, output$current)
})

##run test quantifyConnections
test_that("check quantifyConnections function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    multupletNames <- c("A1A1A1", "A1A1A1", "A1A1B1")
    
    #setup expected data
    expectedVar1 <- c("A1-A1", "A1-B1")
    expectedFreq <- c(7, 2)
    expectedType1 <- rep("A1", 2)
    expectedType2 <- c("A1", "B1")
    
    #run function
    output <- quantifyConnections(multupletNames)
    
    #test
    expect_equivalent(expectedVar1, output$Var1)
    expect_equivalent(expectedFreq, output$Freq)
    expect_equivalent(expectedType1, output$type1)
    expect_equivalent(expectedType2, output$type2)
})

##run test calculateFreqAndSum
test_that("check calculateFreqAndSum function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    i <- 1
    targetConnections <- data.frame(
        conn = c("A1-A1", "A1-B1"),
        target = rep(20, 2),
        current = rep(0, 2)
    )
    current <- quantifyConnections(c("A1A1A1", "A1A1A1", "A1A1B1"))
    
    #setup expected data
    expectedTypes <- c("A1", "A1")
    expectedFreq <- 7
    expectedSum <- 9
    
    #run function
    output <- calculateFreqAndSum(i, targetConnections, current)
    outputTypes <- output[[1]]
    outputFreq <- output[[2]]
    outputSum <- output[[3]]

    #test
    expect_equivalent(expectedTypes, outputTypes)
    expect_equivalent(expectedFreq, outputFreq)
    expect_equivalent(expectedSum, outputSum)
})

##run test calculateNumToAdd
test_that("check calculateNumToAdd function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    i <- 1
    currentFreq <- 0
    currentSum <- 0
    targetConnections <- data.frame(
        conn = c("A1-A1", "A1-B1"),
        target = rep(20, 2),
        current = rep(0, 2)
    )
    
    #setup expected data
    expected <- 0
    
    #run function
    output <- calculateNumToAdd(i, currentFreq, currentSum, targetConnections)
    
    #test
    expect_equivalent(expected, output)
    
    ###TEST2####
    #set up input data
    i <- 1
    currentFreq <- 0
    currentSum <- 100
    targetConnections <- data.frame(
        conn = c("A1-A1", "A1-B1"),
        target = rep(20, 2),
        current = rep(0, 2)
    )
    
    #setup expected data
    expected <- 25
    
    #run function
    output <- calculateNumToAdd(i, currentFreq, currentSum, targetConnections)
    
    #test
    expect_equivalent(expected, output)
})

##run test synthesizeAndAdd
test_that("check synthesizeAndAdd function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    add <- 13
    currentTypes <- c("A1", "B1")
    nGenes <- 1000
    multuplets <- data.frame(row.names = 1:nGenes)
    singlets <- syntheticSinglets(
        nGenes = nGenes,
        nCells = 5,
        nCellTypes = 2,
        singletExpansion = 2
    )[[2]]
    colnames(singlets) <- sort(rep(paste(LETTERS[1:2], 1, sep = ""), 5))
    
    #setup expected data
    expectedTypes <- c("A1", "B1")
    expectedFreq <- 7
    expectedSum <- 9
    
    #run function
    output <- synthesizeAndAdd(
        add,
        currentTypes,
        multuplets,
        singlets
    )
    
    #test
    outputTypes <- colnames(output) %>%
        gsub("(.{2})", "\\1 ", .) %>%
        trimws %>% strsplit(., "\\s+") %>%
        unlist %>%
        unique
    expect_equivalent(expectedTypes, outputTypes)
    
    freq <- quantifyConnections(colnames(output))
    outputFreq <- freq[freq$Var1 == "A1-B1", "Freq"]
    expect_true(all.equal(expectedFreq, outputFreq, tolerance = 5))
    
})

##run test checkFrequency
test_that("check checkFrequency function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    multuplets <- data.frame(1,2,3,4,5,6)
    colnames(multuplets) <- c(rep("A1B1", 3), rep("C1D1", 2), "C1E1")
    targetConnections <- data.frame(
        conn = c("A1-B1", "C1-D1"),
        target = 20,
        current = 0,
        stringsAsFactors = FALSE
    )
    
    #setup expected data
    expectedA1B1 <- 100
    expectedC1D1 <- round(100 * (2 / 3))
    
    #run function
    tmp <- checkFrequency(
        1,
        multuplets,
        targetConnections,
        currentTypes = c("A1", "B1")
    )
    output <- checkFrequency(
        2,
        multuplets,
        tmp,
        currentTypes = c("C1", "D1")
    )
    
    #test
    expect_equivalent(expectedA1B1, output[output$conn == "A1-B1", "current"])
    expect_equivalent(expectedC1D1, output[output$conn == "C1-D1", "current"])
})


################################################################################
#                                                                              #
# Adjust self.                                                                 #
#                                                                              #
################################################################################

test_that("check syntheticMultuplets function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    multuplets <- data.frame(lapply(1:20, function(x) 1:10))
    colnames(multuplets) <- c(
        rep("A1A1", 10),
        rep("B1B1", 10)
    )
    singlets <- .syntheticSinglets(
        nGenes = 10,
        nCells = 5,
        nCellTypes = 2,
        singletExpansion = 2
    )[[2]]
    colnames(singlets) <- c(rep("A1", 5), rep("B1", 5))
    
    #setup expected data
    expectedTotal <- 20
    
    #run function
    output <- adjustSelf(multuplets = multuplets, singlets = singlets)
    outputNewCols <- output[, (ncol(multuplets) - 1):ncol(output)]
    colnames(outputNewCols) <- gsub("(.*)\\.[0-9]", "\\1", colnames(outputNewCols))
    addedConns <- quantifyConnections(colnames(outputNewCols))
    
    #test
    expect_true(all.equal(addedConns[1, "Freq"], expectedTotal, tolerance = 5))
    expect_true(all.equal(addedConns[2, "Freq"], expectedTotal, tolerance = 5))
})

################################################################################
#                                                                              #
# Multuplets.                                                                  #
#                                                                              #
################################################################################

##run test syntheticMultuplets
test_that("check syntheticMultuplets function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    nGenes <- 1000
    nCells <- 100
    nCellTypes <- 4
    perplexity <- 5
    target <- 20
    singletExpansion <- 2
    
    #setup expected data
    expectedClassification <- c("A1", "B1", "C1", "D1")
    
    #run function
    output <- syntheticMultuplets(
        nGenes,
        nCells,
        nCellTypes,
        perplexity,
        target,
        singletExpansion
    )
    sng <- output[[1]]
    mult <- output[[2]]
    uObj <- output[[3]]
    
    #test that classification is correct
    expect_equivalent(expectedClassification, unique(getData(uObj, "classification")))
    #something else I can test here?
})

################################################################################
#                                                                              #
# Test data.                                                                   #
#                                                                              #
################################################################################

##run test syntheticTestData
test_that("check syntheticTestData function outputs the expected result", {
    
    ###TEST1####
    #set up input datanMultiplets,
    nMultiplets <- 100
    nGenes <- 1000
    nCells <- 20
    nCellTypes <- 4
    perplexity <- 5
    target <- 20
    singletExpansion <- 2
    
    #setup expected data
    expectedSngMultTotal <- table(c(rep(FALSE, 15), rep(TRUE, 40)))
    expectedSngMultTest <- table(c(rep(FALSE, 2), rep(TRUE, 40)))

    expectedClassification <- c("A1", "B1", "C1", "D1")
    
    #run function
    output <- syntheticTestData(
        nMultiplets,
        nGenes,
        nCells,
        nCellTypes,
        perplexity,
        target,
        singletExpansion
    )
    syntheticData <- output[[1]]
    testSyntheticData <- output[[2]]
    uObj <- output[[3]]
    table <- output[[4]]
    
    #test
    expect_equivalent(nrow(syntheticData), nGenes)
    expect_equivalent(nrow(testSyntheticData), nGenes)
    #more tests here?
})

