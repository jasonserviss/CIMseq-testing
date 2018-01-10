
#context("funToCheckKnownResults")

##run test .tp
test_that("check that .tp outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  detected <- tibble(
    multiplet = sort(c(rep(LETTERS[1], 2), rep(LETTERS[2], 5))),
    from = c("A", "B", rep("A", 3), rep("B", 2)),
    to = c("B", "C", "B", "C", "D", "C", "D")
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  known <- tibble(
    multiplet = sort(c(LETTERS[1], rep(LETTERS[2], 6))),
    from = c("A", rep("A", 3), rep("B", 2), "C"),
    to = c("B", "B", "C", "D", "C", "D", "D")
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  data <- full_join(detected, known, by = "multiplet")
  
  #setup expected data
  expected <- c(1, 5)
  
  #run function
  output <- .tp(data)
  
  #test
  expect_equivalent(expected, output)
})

##run test .fp
test_that("check that .fp outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  detected <- tibble(
    multiplet = sort(c(rep(LETTERS[1], 2), rep(LETTERS[2], 5))),
    from = c("A", "B", rep("A", 3), rep("B", 2)),
    to = c("B", "C", "B", "C", "D", "C", "D")
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  known <- tibble(
    multiplet = sort(c(LETTERS[1], rep(LETTERS[2], 6))),
    from = c("A", rep("A", 3), rep("B", 2), "C"),
    to = c("B", "B", "C", "D", "C", "D", "D")
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  data <- full_join(detected, known, by = "multiplet")
  
  #setup expected data
  expected <- c(1, 0)
  
  #run function
  output <- .fp(data)
  
  #test
  expect_equivalent(expected, output)
})

##run test .fn
test_that("check that .fn outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  detected <- tibble(
    multiplet = sort(c(rep(LETTERS[1], 2), rep(LETTERS[2], 5))),
    from = c("A", "B", rep("A", 3), rep("B", 2)),
    to = c("B", "C", "B", "C", "D", "C", "D")
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  known <- tibble(
    multiplet = sort(c(LETTERS[1], rep(LETTERS[2], 6))),
    from = c("A", rep("A", 3), rep("B", 2), "C"),
    to = c("B", "B", "C", "D", "C", "D", "D")
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  data <- full_join(detected, known, by = "multiplet")
  
  #setup expected data
  expected <- c(0, 1)
  
  #run function
  output <- .fn(data)
  
  #test
  expect_equivalent(expected, output)
})

##run test .tn
test_that("check that .tn outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  .process <- function(data) {
    data %>%
    getEdgesForMultiplet(., 0.2, rownames(getData(., "spSwarm"))) %>%
    unite(connections, from, to, sep = "-") %>%
    nest(-multiplet)
  }
  
  detected <- new("spSwarm",
    spSwarm = data.frame(
      A = c(0),
      B = c(0.5),
      C = c(0.5),
      row.names = c("A1")
    ),
    costs = vector(mode = "numeric"),
    convergence = vector(mode = "character"),
    stats = list(),
    arguments = list()
  )
  
  known <- new("spSwarm",
    spSwarm = data.frame(
      A = c(0),
      B = c(0.5),
      C = c(0.5),
      row.names = c("A1")
    ),
    costs = vector(mode = "numeric"),
    convergence = vector(mode = "character"),
    stats = list(),
    arguments = list()
  )
  
  data <- full_join(.process(detected), .process(known), by = "multiplet")
  
  #setup expected data
  expected <- c(2)
  
  #run function
  output <- .tn(data, detected, known)
  
  #test
  expect_equivalent(expected, output)
})
