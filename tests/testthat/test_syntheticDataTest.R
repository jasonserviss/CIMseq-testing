#context("syntheticDataTest")

##run test .syntheticSinglets
test_that("check that the .syntheticSinglets function outputs the expected result", {
    
    ###TEST1####
    #setup input
    ngenes=100
    ncells=10
    cellTypes=3
    
    #setup expected data
    expectedDim <- c(100, 30)
    expectedNames <- paste(
        sort(
            rep(
                letters,
                ncells
            )
        )[1:(cellTypes*ncells)],
        1:ncells,
        sep=""
    )
    expectedClass <- "data.frame"
    expectedMin <- 0
    expectedMax <- 109
    
    #run function
    output <- .syntheticSinglets(ngenes, ncells, cellTypes)
    
    #test
    expect_equivalent(dim(output), expectedDim)
    expect_equivalent(colnames(output), expectedNames)
    expect_equivalent(class(output), expectedClass)
    expect_equivalent(round(min(output)), expectedMin)
    expect_equivalent(round(max(output)), expectedMax)
})

##run test .syntheticMultuplets
test_that("check that the .syntheticMultuplets function outputs the expected result", {
    
    ###TEST1####
    #set up input data
    ngenes=1000
    ncells=10
    cellTypes=4
    perplexity=5
    
    #setup expected data
    expectedClassification <- c("A1", "B1", "C1", "D1")
    
    #run function
    output <- .syntheticMultuplets(ngenes, ncells, cellTypes, perplexity)
    sng <- output[[1]]
    mult <- output[[2]]
    uObj <- output[[3]]
    
    #test
    expect_equivalent(expectedClassification, unique(getData(uObj, "classification")))
    expect_equivalent(rowMeans(sng[,which(colnames(sng) %in% 'A1')]), mult[ ,'A1A1'])
    expect_equivalent(rowMeans(sng[,which(colnames(sng) %in% c('A1', 'B1'))]), mult[ ,'A1B1'])
    expect_equivalent(rowMeans(sng[,which(colnames(sng) %in% c('A1', 'B1', 'C1'))]), mult[ ,'A1B1C1'])
    expect_equivalent(rowMeans(sng[,which(colnames(sng) %in% c('A1', 'B1', 'C1', 'D1'))]), mult[ ,'A1B1C1D1'])

})

##run test .syntheticTestData
test_that("check that the .syntheticTestData function outputs the expected result", {
    
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

