#context("multipletSynthMethods")

################################################################################
#                                                                              #
# Run test .makeMultCombos                                                     #
#                                                                              #
################################################################################

test_that("test unit .makeMultCombos with normal input and default args", {
    
    #prepare normal input data
    nCellsInMultiplet <- 2
    cellTypes <- c("A1", "B1", "C1")
    self <- TRUE
    repetitions <- NULL
    
    #setup expected outcome
    expected <- matrix(
        c(
            "A1", "A1", "A1", "B1", "B1", "C1",
            "A1", "B1", "C1", "B1", "C1", "C1"
        ),
        byrow = TRUE,
        nrow = 2
    )
    
    #run function
    output <- .makeMultCombos(
        nCellsInMultiplet,
        cellTypes,
        self,
        repetitions
    )
    
    #test
    expect_identical(output, expected)
})

test_that("test unit .makeMultCombos with normal input and self = FALSE", {
    
    #prepare normal input data
    nCellsInMultiplet <- 2
    cellTypes <- c("A1", "B1", "C1")
    self <- FALSE
    repetitions <- NULL
    
    #setup expected outcome
    expected <- matrix(
        c(
            "A1", "A1", "B1",
            "B1", "C1", "C1"
        ),
        byrow = TRUE,
        nrow = 2
    )
    
    #run function
    output <- .makeMultCombos(
        nCellsInMultiplet,
        cellTypes,
        self,
        repetitions
    )
    
    #test
    expect_identical(output, expected)
})

test_that("test unit .makeMultCombos with normal input and repetitions != NULL", {
    
    #prepare normal input data
    nCellsInMultiplet <- 2
    cellTypes <- c("A1", "B1", "C1")
    self <- FALSE
    repetitions <- 2
    
    #setup expected outcome
    expected <- matrix(
        c(
            "A1", "A1", "B1",
            "B1", "C1", "C1"
        ),
        byrow = TRUE,
        nrow = 2
    )
    expected <- cbind(expected, expected)
    
    #run function
    output <- .makeMultCombos(
        nCellsInMultiplet,
        cellTypes,
        self,
        repetitions
    )
    
    #test
    expect_identical(output, expected)
})

################################################################################
#                                                                              #
# Run test .readPool                                                           #
#                                                                              #
################################################################################

test_that("test unit .readPool with normal input", {
    
    #prepare normal input data
    cellTypes <- c("A1", "B1", "C1")
    rownames <- letters[1:5]
    colnames <- c("A1", "B1", "C1", "A1", "B1", "C1")
    singlets <- matrix(
        c(
            rep(1, 5),
            rep(2, 5),
            rep(3, 5),
            rep(1, 5),
            rep(2, 5),
            rep(3, 5)
        ),
       ncol = 6,
       dimnames = list(rownames, colnames)
    )
    
    #setup expected outcome
    expected <- list(
        A1 = rep(letters[1:5], 2),
        B1 = c(sort(rep(letters[1:5], 2)), sort(rep(letters[1:5], 2))),
        C1 = c(sort(rep(letters[1:5], 3)), sort(rep(letters[1:5], 3)))
    )
    names(expected) <- cellTypes
    
    #run function
    output <- .readPool(cellTypes, singlets)
    
    #test
    expect_identical(output, expected)
})

test_that("test unit .readPool with missing cell types", {
    
    #TEST 1
    #prepare input data
    cellTypes <- c("A1", "B1")
    rownames <- letters[1:5]
    colnames <- c("A1", "B1", "C1", "A1", "B1", "C1")
    singlets <- matrix(
        c(
            rep(1, 5),
            rep(2, 5),
            rep(3, 5),
            rep(1, 5),
            rep(2, 5),
            rep(3, 5)
        ),
       ncol = 6,
       dimnames = list(rownames, colnames)
    )
    
    #test
    expect_error(.readPool(cellTypes, singlets))
    
    #TEST 2
    #prepare input data
    cellTypes <- c("A1", "B1", "C1")
    rownames <- letters[1:5]
    colnames <- c("A1", "B1", "A1", "B1")
    singlets <- matrix(
        c(
            rep(1, 5),
            rep(2, 5),
            rep(1, 5),
            rep(2, 5)
        ),
       ncol = 4,
       dimnames = list(rownames, colnames)
    )
    
    #test
    expect_error(.readPool(cellTypes, singlets))

})

test_that("test unit .readPool with missing output", {
    
    #TEST 1
    #prepare input data
    cellTypes <- c("A1", "B1")
    rownames <- letters[1:5]
    colnames <- c("A1", "B1", "C1", "A1", "B1", "C1")
    singlets <- matrix(
        c(
            rep(0, 5),
            rep(2, 5),
            rep(3, 5),
            rep(0, 5),
            rep(2, 5),
            rep(3, 5)
        ),
       ncol = 6,
       dimnames = list(rownames, colnames)
    )
    
    #test
    expect_error(.readPool(cellTypes, singlets))
})

################################################################################
#                                                                              #
# Run test .meanCountsPerType                                                  #
#                                                                              #
################################################################################

test_that("test unit .meanCountsPerType with normal input", {
    
    #prepare normal input data
    cellTypes <- c("A1", "B1", "C1")
    rownames <- letters[1:5]
    colnames <- c("A1", "B1", "C1", "A1", "B1", "C1")
    singlets <- matrix(
        c(
            rep(1, 5),
            rep(2, 5),
            rep(3, 5),
            rep(1, 5),
            rep(2, 5),
            rep(3, 5)
        ),
       ncol = 6,
       dimnames = list(rownames, colnames)
    )
    
    #setup expected outcome
    expected <- c(5, 10, 15)
    names(expected) <- cellTypes
    
    #run function
    output <- .meanCountsPerType(cellTypes, singlets)
    
    #test
    expect_identical(output, expected)
})

test_that("test unit .meanCountsPerType with only one sample in one cell type",{
    
    #prepare normal input data
    cellTypes <- c("A1", "B1", "C1")
    colnames <- c("A1", "A1", "B1", "B1", "C1")
    rownames <- letters[1:5]
    singlets <- matrix(
        c(
            rep(1, 5),
            rep(1, 5),
            rep(2, 5),
            rep(2, 5),
            rep(3, 5)
        ),
       ncol = 5,
       dimnames = list(rownames, colnames)
    )
    
    #test
    expect_error(.meanCountsPerType(cellTypes, singlets))
})

################################################################################
#                                                                              #
# Run test .makeMultipletsList                                                 #
#                                                                              #
################################################################################

test_that("test .makeMultipletsList with self connection input and doublet", {
    
    #prepare normal input data
    combos <- matrix(c("A1", "A1"), nrow = 2)
    fractions <- tibble(cellTypes = c("A1", "B1"), fractions = rep(1, 2))
    pool <- list(
        A1 = rep(letters[1:5], 2),
        B1 = rep(letters[6:10], 3)
    )
    means <- c(5, 5)
    
    #setup expected outcome
    expected <- list(table(.=c("b", "c")))
    
    #run function
    output <- .makeMultipletsList(combos, fractions, pool, means)
    
    #test
    expect_identical(output, expected)
    expect_false(all(names(unlist(output)) %in% pool[[2]]))
    expect_true(all(names(unlist(output)) %in% pool[[1]]))
})

test_that("test .makeMultipletsList with non-self connection input and doublet", {
    
    #prepare normal input data
    combos <- matrix(c("A1", "B1"), nrow = 2)
    fractions <- tibble(cellTypes = c("A1", "B1"), fractions = rep(1, 2))
    pool <- list(
        A1 = rep(letters[1:5], 2),
        B1 = rep(letters[6:10], 3)
    )
    means <- c(5, 5)
    
    #setup expected outcome
    expected <- list(table(.=c("c", "h")))
    
    #run function
    output <- .makeMultipletsList(combos, fractions, pool, means)
    
    #test
    expect_identical(output, expected)
    expect_true(all(names(unlist(output)) %in% unlist(pool)))
})

test_that("test .makeMultipletsList with non-self connection input and triplet", {
    
    #prepare normal input data
    combos <- matrix(c("A1", "B1", "C1"), nrow = 3)
    fractions <- tibble(cellTypes = c("A1", "B1", "C1"), fractions = rep(1, 3))
    pool <- list(
        A1 = rep(letters[1:5], 2),
        B1 = rep(letters[6:10], 3),
        C1 = rep(letters[11:20], 2)
    )
    means <- c(10, 10, 10)
    
    #setup expected outcome
    expected <- list(table(.=c("c", "h", "n")))
    
    #run function
    output <- .makeMultipletsList(combos, fractions, pool, means)
    
    #test
    expect_identical(output, expected)
    expect_true(all(names(unlist(output)) %in% unlist(pool)))
})

################################################################################
#                                                                              #
# Run test .makeMultipletsList                                                 #
#                                                                              #
################################################################################

test_that("test .getCurrentFractions with normal input #1", {
    
    #prepare normal input data
    combos <- matrix(c("A1", "B1", "C1", "D1"), nrow = 2)
    fractions <- tibble(
        cellTypes = c("A1", "B1", "C1"),
        fractions = c(0.5, 0.5, 0.5)
    )
    x <- 1
    
    #setup expected outcome
    expected <- c(0.5, 0.5)
    
    #run function
    output <- .getCurrentFractions(combos, fractions, x)
    
    #test
    expect_identical(output, expected)
    expect_true(sum(output) == 1)
})

test_that("test .getCurrentFractions with normal input #2", {
    
    #prepare normal input data
    combos <- matrix(c("A1", "B1", "C1", "D1"), nrow = 2)
    fractions <- tibble(
        cellTypes = c("A1", "B1", "C1"),
        fractions = c(0.2, 0.8, 0.5)
    )
    x <- 1
    
    #setup expected outcome
    expected <- c(0.2, 0.8)
    
    #run function
    output <- .getCurrentFractions(combos, fractions, x)
    
    #test
    expect_identical(output, expected)
    expect_true(sum(output) == 1)
})

test_that("test .getCurrentFractions with normal input #3", {
    
    #prepare normal input data
    combos <- matrix(c("A1", "B1", "C1", "D1"), nrow = 2)
    fractions <- tibble(
        cellTypes = c("A1", "B1", "C1"),
        fractions = c(0.2, 0.2, 0.5)
    )
    x <- 1
    
    #setup expected outcome
    expected <- c(0.5, 0.5)
    
    #run function
    output <- .getCurrentFractions(combos, fractions, x)
    
    #test
    expect_identical(output, expected)
    expect_true(sum(output) == 1)
})

################################################################################
#                                                                              #
# Run test .convertMultipletListToMatrix                                       #
#                                                                              #
################################################################################

test_that("test .convertMultipletListToMatrix with normal input.", {
    
    #prepare normal input data
    combos <- matrix(c("A1", "B1", "C1", "D1"), nrow = 2)
    countsList <- list(
        table(rep(letters[1:3], 2)),
        table(rep(letters[4:6], 2))
    )
    
    #setup expected outcome
    colnames <- paste(combos[1, ], combos[2, ], sep = "-")
    rownames <- letters[1:6]
    expected <- matrix(c(
        c(rep(2, 3), rep(0, 3)),
        c(rep(0, 3), rep(2, 3))
    ), dimnames = list(rownames, colnames), ncol = 2)
    
    
    #run function
    output <- .convertMultipletListToMatrix(countsList, combos)
    
    #test
    expect_identical(output, expected)
})
