
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
    multiplet = sort(c(rep(LETTERS[1], 2), rep(LETTERS[2], 6))),
    from = c("A", "B", rep("A", 3), rep("B", 2), "C"),
    to = c("B", "C", "B", "C", "D", "C", "D", "D")
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  data <- full_join(detected, known, by = "multiplet")
  
  #setup expected data
  expected <- c(2, 5)
  
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
    to = c("B", "D", "B", "C", "D", "C", "D")
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  known <- tibble(
    multiplet = sort(c(rep(LETTERS[1], 2), rep(LETTERS[2], 6))),
    from = c("A", "B", rep("A", 3), rep("B", 2), "C"),
    to = c("B", "C", "B", "C", "D", "C", "D", "D")
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
    multiplet = sort(c(rep(LETTERS[1], 2), rep(LETTERS[2], 6))),
    from = c("A", "B", rep("A", 3), rep("B", 2), "C"),
    to = c("B", "C", "B", "C", "D", "C", "D", "D")
  ) %>%
  unite(connections, from, to, sep = "-") %>%
  nest(-multiplet)
  
  data <- full_join(detected, known, by = "multiplet")
  
  #setup expected data
  expected <- c(1, 0)
  
  #run function
  output <- .fn(data)
  
  #test
  expect_equivalent(expected, output)
})

#WRITE CHECK FOR .tn
