#context("syntheticDataTest")

################################################################################
#                                                                              #
# Singlets.                                                                    #
#                                                                              #
################################################################################

##run test .syntheticSinglets
test_that("check .syntheticSinglets function outputs the expected result", {
    
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

##run test .makeMultuplet
test_that("check .makeMultuplet with doublets", {
    
    #set up input data
    ##arguments to .syntheticSinglets
    cellTypes <- paste(LETTERS[1:3], 1, sep = "")
    nGenes <- 100
    nCells <- 3
    nCellTypes <- length(cellTypes)
    singletExpansion <- 2
    

    ##arguments to .makeMultuplet
    nCellsInMultiplet <- 2
    multuplets <- data.frame(row.names = 1:nrow(singlets))
    singlets <- .syntheticSinglets(
        nGenes,
        nCells,
        nCellTypes,
        singletExpansion
    )[[2]]
    colnames(singlets) <- sort(rep(cellTypes, 3))
    
    #setup expected data
    expectedDim <- c(nGenes, sum(length(cellTypes), ncol(combn(cellTypes, 2))))
    expectedMultuplets <- c("A1A1", "A1B1", "A1C1", "B1B1", "B1C1", "C1C1")
    
    #run function
    output <- .makeMultuplet(
        nCellsInMultiplet,
        cellTypes,
        multuplets,
        singlets
    )
    
    #test
    expect_equivalent(dim(output), expectedDim)
    expect_equivalent(colnames(output), expectedMultuplets)
})

test_that("check .makeMultuplet with triplets", {
    
    #set up input data
    ##arguments to .syntheticSinglets
    cellTypes <- paste(LETTERS[1:3], 1, sep = "")
    nGenes <- 100
    nCells <- 3
    nCellTypes <- length(cellTypes)
    singletExpansion <- 2
    
    ##arguments to .makeMultuplet
    nCellsInMultiplet <- 3
    multuplets <- data.frame(row.names = 1:nrow(singlets))
    singlets <- .syntheticSinglets(
        nGenes,
        nCells,
        nCellTypes,
        singletExpansion
    )[[2]]
    colnames(singlets) <- sort(rep(cellTypes, 3))
    
    #setup expected data
    expectedDim <- c(nGenes, sum(length(cellTypes), 7))
    expectedMultuplets <- c(
        "A1A1A1", "A1A1B1", "A1A1C1", "A1B1B1", "A1B1C1",
        "A1C1C1", "B1B1B1", "B1B1C1", "B1C1C1", "C1C1C1"
    )
    
    #run function
    output <- .makeMultuplet(
        nCellsInMultiplet,
        cellTypes,
        multuplets,
        singlets
    )
    
    #test
    expect_equivalent(dim(output), expectedDim)
    expect_equivalent(colnames(output), expectedMultuplets)
})

test_that("check .makeMultuplet with quadruplets", {
    
    #set up input data
    ##arguments to .syntheticSinglets
    cellTypes <- paste(LETTERS[1:3], 1, sep = "")
    nGenes <- 100
    nCells <- 3
    nCellTypes <- length(cellTypes)
    singletExpansion <- 2
    
    
    ##arguments to .makeMultuplet
    nCellsInMultiplet <- 4
    multuplets <- data.frame(row.names = 1:nrow(singlets))
    singlets <- .syntheticSinglets(
        nGenes,
        nCells,
        nCellTypes,
        singletExpansion
    )[[2]]
    colnames(singlets) <- sort(rep(cellTypes, 3))
    
    #setup expected data
    expectedDim <- c(nGenes, sum(length(cellTypes), 12))
    expectedMultuplets <- c(
        "A1A1A1A1", "A1A1A1B1", "A1A1A1C1", "A1A1B1B1", "A1A1B1C1",
        "A1A1C1C1", "A1B1B1B1", "A1B1B1C1", "A1B1C1C1", "A1C1C1C1",
        "B1B1B1B1", "B1B1B1C1", "B1B1C1C1", "B1C1C1C1", "C1C1C1C1"
    )
    
    #run function
    output <- .makeMultuplet(
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

##run test .decideConnections
test_that("check .decideConnections function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    cellTypes <- paste(LETTERS[1:4], 1, sep = "")
    target <- 20
    
    #setup expected data
    expectedConn <- c("A1-D1", "B1-C1")
    expectedTarget <- rep(target, 2)
    expectedCurrent <- rep(0, 2)
    
    #run function
    output <- .decideConnections(cellTypes, target)
    
    #test
    expect_equivalent(expectedConn, output$conn)
    expect_equivalent(expectedTarget, output$target)
    expect_equivalent(expectedCurrent, output$current)
})

##run test .quantifyConnections
test_that("check .quantifyConnections function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    multupletNames <- c("A1A1A1", "A1A1A1", "A1A1B1")
    
    #setup expected data
    expectedVar1 <- c("A1-A1", "A1-B1")
    expectedFreq <- c(7, 2)
    expectedType1 <- rep("A1", 2)
    expectedType2 <- c("A1", "B1")
    
    #run function
    output <- .quantifyConnections(multupletNames)
    
    #test
    expect_equivalent(expectedVar1, output$Var1)
    expect_equivalent(expectedFreq, output$Freq)
    expect_equivalent(expectedType1, output$type1)
    expect_equivalent(expectedType2, output$type2)
})

##run test .calculateFreqAndSum
test_that("check .calculateFreqAndSum function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    i <- 1
    targetConnections <- data.frame(
        conn = c("A1-A1", "A1-B1"),
        target = rep(20, 2),
        current = rep(0, 2)
    )
    current <- .quantifyConnections(c("A1A1A1", "A1A1A1", "A1A1B1"))
    
    #setup expected data
    expectedTypes <- c("A1", "A1")
    expectedFreq <- 7
    expectedSum <- 9
    
    #run function
    output <- .calculateFreqAndSum(i, targetConnections, current)
    outputTypes <- output[[1]]
    outputFreq <- output[[2]]
    outputSum <- output[[3]]

    #test
    expect_equivalent(expectedTypes, outputTypes)
    expect_equivalent(expectedFreq, outputFreq)
    expect_equivalent(expectedSum, outputSum)
})

##run test .calculateNumToAdd
test_that("check .calculateNumToAdd function outputs the expected result", {
    
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
    output <- .calculateNumToAdd(i, currentFreq, currentSum, targetConnections)
    
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
    output <- .calculateNumToAdd(i, currentFreq, currentSum, targetConnections)
    
    #test
    expect_equivalent(expected, output)
})

##run test .synthesizeAndAdd
test_that("check .synthesizeAndAdd function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    add <- 2
    currentTypes <- c("A1", "B1")
    nGenes <- 1000
    multuplets <- data.frame(row.names = 1:nGenes)
    nCellTypes <- 3
    
    #setup expected data
    expectedTypes <- c("A1", "A1")
    expectedFreq <- 7
    expectedSum <- 9
    
    #run function
    output <- .synthesizeAndAdd(add, currentTypes, nGenes, multuplets, nCellTypes)
    
    #test
    expect_equivalent(expectedTypes, outputTypes)
    expect_equivalent(expectedFreq, outputFreq)
    expect_equivalent(expectedSum, outputSum)
})

################################################################################
#                                                                              #
# Adjust self.                                                                 #
#                                                                              #
################################################################################

test_that("check .syntheticMultuplets function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    
    #setup expected data
    
    #run function
    
    #test
    expect_equivalent(expectedClassification, unique(getData(uObj, "classification")))
    
})

################################################################################
#                                                                              #
# Multuplets.                                                                  #
#                                                                              #
################################################################################

##run test .syntheticMultuplets
test_that("check .syntheticMultuplets function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    ngenes <- 1000
    ncells <- 100
    cellTypes <- 4
    perplexity <- 5
    target <- 20
    
    #setup expected data
    expectedClassification <- c("A1", "B1", "C1", "D1")
    
    #run function
    output <- .syntheticMultuplets(
        ngenes,
        ncells,
        cellTypes,
        perplexity,
        target
    )
    sng <- output[[1]]
    mult <- output[[2]]
    uObj <- output[[3]]
    
    #test that classification is correct
    expect_equivalent(expectedClassification, unique(getData(uObj, "classification")))
    
    #
    expect_equivalent(rowMeans(sng[, which(colnames(sng) %in% 'A1')]), mult[ ,'A1A1'])
    expect_equivalent(rowMeans(sng[, which(colnames(sng) %in% c('A1', 'B1'))]), mult[ ,'A1B1'])
    expect_equivalent(rowMeans(sng[, which(colnames(sng) %in% c('A1', 'B1', 'C1'))]), mult[ ,'A1B1C1'])
    expect_equivalent(rowMeans(sng[, which(colnames(sng) %in% c('A1', 'B1', 'C1', 'D1'))]), mult[ ,'A1B1C1D1'])
    
})

################################################################################
#                                                                              #
# Test data.                                                                   #
#                                                                              #
################################################################################

##run test .syntheticTestData
test_that("check .syntheticTestData function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    n=2
    ngenes=1000
    ncells=10
    cellTypes=4
    perplexity=5
    
    #setup expected data
    expectedSngMultTotal <- table(c(rep(FALSE, 15), rep(TRUE, 40)))
    expectedSngMultTest <- table(c(rep(FALSE, 2), rep(TRUE, 40)))

    expectedClassification <- c("A1", "B1", "C1", "D1")
    
    #run function
    output <- .syntheticTestData(n, ngenes, ncells, cellTypes, perplexity)
    syntheticData <- output[[1]]
    testSyntheticData <- output[[2]]
    uObj <- output[[3]]
    table <- output[[4]]
    
    #test
    expect_equivalent(ncol(syntheticData), 55)
    expect_equivalent(ncol(testSyntheticData), (ncells*cellTypes)+n)
    expect_equivalent(nrow(syntheticData), ngenes)
    expect_equivalent(nrow(testSyntheticData), ngenes)
    expect_equivalent(table(nchar(colnames(syntheticData)) == 4), expectedSngMultTotal)
    expect_equivalent(table(nchar(colnames(testSyntheticData)) == 4), expectedSngMultTest)
    expect_equivalent(unique(getData(uObj, "classification")), expectedClassification)
    expect_false(max(table) > 9)
    
    
})

