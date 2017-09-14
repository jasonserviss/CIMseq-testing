#context("expressionDataTest")

##run test .assemble
test_that("check that .assemble outputs the expected result for hetero type arg", {
    
    ###TEST1####
    #prepare normal input data
    A <- c(1,2,3)
    B <- c(10,20,30)
    C <- c(100,200,300)
    
    counts <- data.frame(
        A = A,
        A = A,
        B = B,
        B = B,
        C = C,
        C = C
    )
    colnames(counts) <- sort(rep(LETTERS[1:3], 2))
    classification <- sort(rep(LETTERS[1:3], 2))
    dataset <- data.frame(row.names = 1:3)
    names <- c()
    x = 2
    type = "hetero"

    #setup expected data
    expected1 <- data.frame(
        `1` = rowMeans(counts[, 2:3]),
        `2` = rowMeans(counts[, c(1, 6)]),
        `3` = rowMeans(counts[, 4:5])
    )
    colnames(expected1) <- as.character(1:3)
    expected2 <- c("AB", "AC", "BC")
    
    #run function
    tmp <- .assemble(counts, classification, dataset, names, x, type)
    output1 <- tmp[[1]]
    output2 <- tmp[[2]]
    
    #test
    expect_equivalent(expected1, output1)
    expect_equivalent(expected2, output2)

})

##run test .assemble
test_that("check that .assemble outputs the expected result for homo type arg", {
    
    ###TEST1####
    #prepare normal input data
    A <- c(1,2,3)
    B <- c(10,20,30)
    C <- c(100,200,300)
    
    counts <- data.frame(
        A = A,
        A = A,
        B = B,
        B = B,
        C = C,
        C = C
    )
    colnames(counts) <- sort(rep(LETTERS[1:3], 2))
    classification <- sort(rep(LETTERS[1:3], 2))
    dataset <- data.frame(row.names = 1:3)
    names <- c()
    x = 2
    type = "homo"
    
    #setup expected data
    expected1 <- data.frame(
        `1` = rowMeans(counts[, 1:2]),
        `2` = rowMeans(counts[, 3:4]),
        `3` = rowMeans(counts[, 5:6])
    )
    colnames(expected1) <- as.character(1:3)
    expected2 <- c("AA", "BB", "CC")
    
    #run function
    tmp <- .assemble(counts, classification, dataset, names, x, type)
    output1 <- tmp[[1]]
    output2 <- tmp[[2]]
    
    #test
    expect_equivalent(expected1, output1)
    expect_equivalent(expected2, output2)
})
