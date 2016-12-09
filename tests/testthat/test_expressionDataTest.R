#context("expressionDataTest")

##run test .assemble
test_that("check that the .assemble function outputs the expected result", {
    
    ###TEST1####
    #prepare normal input data
    A <- c(1,2,3)
    B <- c(10,20,30)
    C <- c(100,200,300)
    
    means <- data.frame(A=A, B=B, C=C)
    classification <- c("A", "B", "C")
    dataset <- data.frame(row.names=1:3)
    names <- c()
    x=2
    type="hetero"

    #setup expected data
    expected1 <- data.frame(
        currentMult = c(6, 11, 16),
        currentMult = c(50, 101, 152),
        currentMult = c(55, 110, 165)
    )
    
    expected2 <- c("AB", "AC", "BC")
    
    #run function
    tmp <- .assemble(means, classification, dataset, names, x, type)
    output1 <- tmp[[1]]
    output2 <- tmp[[2]]
    
    #test
    expect_equivalent(expected1, round(output1))
    expect_equivalent(expected2, output2)

})

