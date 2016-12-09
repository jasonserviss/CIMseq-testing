#context("generalFunctions")

##run test .noNoise
test_that("check that the .pickMultuplets function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    multuplets <- matrix(
        rep(NA, 12),
        ncol=4
    )
    names <- c("m.A1B1", "m.B1C1", "m.C1D1", "m.A1B1C1D1")
    colnames(multuplets) <- names
    
    nCells <- 'two'
    nMultuplets <- 3
    
    #setup expected data
    expected <- names[1:3]
    
    #run function
    output <- .pickMultuplets(multuplets, nCells, nMultuplets)
    
    #test
    expect_true(all(expected %in% output == TRUE))
    
})

##run test .noise
test_that("check that the .noise function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    noise <- array(NA, dim=c(10, 3, 2, 2))
    noise[,,1,1] <- matrix(c(rep(1,15), rep(2, 15)), ncol=3)
    dimnames(noise)[[3]] <- as.character(c(0,0.2))
    dimnames(noise)[[4]] <- c('two', 'three')
    
    percent <- 0.2
    nCells <- 'two'
    
    #setup expected data
    expected <- table(c(rep(TRUE, 8), rep(FALSE, 2)))
    
    #run function
    output <- .noise(noise, percent, nCells)
    out1 <- output[,,1,1]
    out2 <- output[,,2,1]
    
    out3 <- table(sapply(
        1:nrow(out1),
        function(p)
            identical(out1[p,], out2[p,])
    ))
    
    #test
    expect_equivalent(out1, noise[,,1,1])
    expect_equivalent(expected, out3)
    
})

##run test .subsetRows
test_that("check that the .subsetRows function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    justTheseMult <- matrix(
        rep(NA, 400),
        ncol=4
    )

    percent <- 0.1
    i <- 1
    
    #setup expected data
    expected1 <- 10
    
    #run function
    output <- .subsetRows(justTheseMult, i, percent)
    out1 <- output[[1]]
    out2 <- output[[2]]
    
    #test
    expect_equivalent(length(out1), expected1)
    expect_equivalent(nrow(out2), expected1)
    expect_equivalent(length(out1), nrow(out2))
    
})

##run test .addNoise
test_that("check that the .addNoise function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    tmp <- matrix(
        rep(NA, 40),
        ncol=4
    )
    
    sub <- matrix(
        sample(1:2, 8, replace=TRUE),
        ncol=4
    )
    
    picked <- c(1,2)
    
    #setup expected data
    expected1 <- table(c(rep(TRUE, 32), rep(FALSE, 8)))
    expected2 <- table(rep(FALSE, 8))
    
    #run function
    output <- .addNoise(sub, tmp, picked)
    
    #test
    expect_equivalent(table(is.na(output)), expected1)
    expect_equivalent(table(is.na(output[1:2,])), expected2)
    expect_equivalent(nrow(output), nrow(tmp))
    
})

##run test .generateNoiseData
test_that("check that the .generateNoiseData function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    n <- 15
    
    #setup expected data
    expected <- as.integer(c(22673, 15, 11, 3))
    
    #run function
    output <- .generateNoiseData(n, "x")[[1]]
    
    #test
    expect_equivalent(dim(output), expected)
    
    
    ####test that the names, noise, and syntheticData correspond for the noise level 0  
})