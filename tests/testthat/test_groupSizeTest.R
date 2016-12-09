#context("groupSizeTest")

##run test .makeCounts
test_that("check that the .makeCounts function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    sng <- matrix(
        rep(NA, 32),
        ncol=8
    )
    colnames(sng) <- rep(LETTERS[1:4],2)
    sngNames <- unique(colnames(sng))

    mult <- matrix(rep(1, 4), ncol=1)
    i <- 1
    s <- 2
    
    #run function
    output <- .makeCounts(sngNames, sng, mult, i, s)
    
    #test
    expect_equivalent(colnames(output), c(sort(rep(LETTERS[1:4], 2)), ""))
    expect_equivalent(output[ ,9], rep(1,4))
    expect_equivalent(table(is.na(output)), table(c(rep(FALSE, 4), rep(TRUE, 32))))

})
