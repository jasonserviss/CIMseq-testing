
edgeStats <- spSwarmPoisson(sObj, 0.01)
edgeStats <- edgeStats[, 1:3]
edgeStats <- filter(edgeStats, from != to) #removes self connections from total edges

# Calculating the total number of edges in the dataset
total.edges <- sum(edgeStats$weight) * 2

#calculate expected
#how can I calculate this without permuting
#multiply the total edges with the joint probability for cell type A and cell type B
#edgeStats$expected.edges <- total.edges * (table(classes) / length(classes))[edgeStats$to] * (table(classes) / length(classes))[edgeStats$from]

#If each cell has the ability to bind to each other cell, that also means that each other cell has the possibility to
#bind to that cell. Thus, each cell will have (ncells - 1) *2 connections. Note that the minus 1 removes the self connection.
#Since there are length(classes) = 1325 cells, we expect 2648 per cell.
#names <- paste(classes, 1:length(classes), sep = ".")
#grid <- expand.grid(names, names, stringsAsFactors = FALSE)
#colnames(grid) <- c("from", "to")
#grid <- filter(grid, from != to)
#filter(grid, from == "SI.Paneth.58") %>% nrow(.)
#filter(grid, to == "SI.Paneth.58") %>% nrow(.)
#filter(grid, from == "SI.Paneth.58" | to == "SI.Paneth.58") %>% nrow(.)

classes <- getData(cObjSng, "classification")
niter <- 10000

grid <- expand.grid(classes, classes, stringsAsFactors = FALSE)
colnames(grid) <- c("from", "to")
grid <- filter(grid, from != to) #removes self connections from possible edges

cmbs <- expand.grid(unique(classes), unique(classes))
cmbs <- filter(cmbs, Var1 != Var2)
cmbs <- paste(cmbs[, 1], cmbs[, 2], sep = "-")

calculateExpectedEdges <- function(
  grid, iter, total.edges, cmbs, seed = 2383, calcMeans = TRUE
){
  perms <- map_dfr(1:iter, function(i) {
    #randomly pick indices in the grid to assign edges to
    set.seed(seed + i)
    idx <- sample(1:nrow(grid), total.edges, replace = FALSE)

    #count the class edges
    t <- table(paste(grid[idx, 1], grid[idx, 2], sep = "-"))

    #add 0's
    out <- vector(mode = "logical", length = length(cmbs))
    names(out) <- cmbs
    out[match(names(t), names(out))] <- t
    out[!out] <- 0

    #reformat and return data
    data.frame(combination = names(out), Freq = out, stringsAsFactors = FALSE)
  }, .id = "iteration") %>%
    separate(combination, into = c("from", "to"), sep = "-")

  if(calcMeans) {
    perms <- perms %>%
      group_by(from, to) %>%
      summarize(expected.edges = mean(Freq)) %>%
      ungroup()
    return(perms)
  } else {
    return(perms)
  }
}

#check n iter
# test.iter <- map_dbl(c(500, 1000, 10000, 100000), function(i) {
#   res <- calculateExpectedEdges(grid, i, total.edges, cmbs)
#   filter(res, from == "SI.Paneth" & to == "SI.Stem")$expected.edges
# })

expected <- calculateExpectedEdges(grid, niter, total.edges, cmbs, calcMeans = FALSE) %>%
  as_tibble() %>%
  select(-iteration) %>%
  nest(Freq) %>%
  mutate(expected.edges = map_dbl(data, ~mean(.x[[1]])))
edgeStats <- inner_join(edgeStats, expected, by = c("from", "to"))

#calculate expected without permutation
#calculating total.edges * P(cell type A, cell type B) underestimates expected 
#edges compared to permuted. It would seem that the reason for this is that 
#P(Blood) * P(Blood) gives P(Blood, Blood) and thus, the probabilities "account"
#for self connections. Since self connections are not included in the 
#permutation, P(A, B) overestimates the desired P.
#I fix this by calculating the P(A, B) for all possible interactions, removing
#the probabilities for the self interactions, and restandardizing to [0, 1].
ct.freq <- as.numeric(table(classes) / length(classes))
names(ct.freq) <- names(table(classes))

allProbs <- expand.grid(names(ct.freq), names(ct.freq))
allProbs$p1 <- ct.freq[allProbs$Var1]
allProbs$p2 <- ct.freq[allProbs$Var2]
allProbs$jp <- allProbs$p1 * allProbs$p2
allProbs <- filter(allProbs, Var1 != Var2)
allProbs$jp <- allProbs$jp / sum(allProbs$jp)
edgeStats <- mutate(edgeStats, expected.edges.calculated = map2_dbl(from, to, function(f, t) {
    total.edges * allProbs[allProbs$Var2 == f & allProbs$Var1 == t, "jp"]
  }))


#calculate p-value based on permutations
edgeStats <- edgeStats %>%
  mutate(p.permute = map2_dbl(data, weight, function(p, w) {
    p <- sum(p[[1]] >= w) / length(p[[1]])
    if_else(
      p == 0, 10^-log10(niter), p)
  }))

#add CI 95 for expected.edges
edgeStats <- mutate(
  edgeStats,
  CI95h = map_dbl(expected.edges, ~qchisq(0.975, 2 * (.x + 1)) / 2)
)

# Calculate p-value based on expected/observed
edgeStats$p.pois <- ppois(q = edgeStats$weight, lambda = edgeStats$expected.edges.calculated, lower.tail = FALSE)

edgeStats %>%
  ggplot() + 
  geom_point(aes(-log10(p.permute), -log10(p.pois))) + 
  labs(x = "p-value (permuted)", y = "p-value (calculated)")

#calculate score
edgeStats$score <- edgeStats$weight / edgeStats$expected.edges

##Plots
edgeStats %>%
  filter(from == "Blood") %>%
  unnest() %>%
  ggplot() +
  geom_histogram(aes(Freq), binwidth = 1) +
  facet_wrap(~ to) +
  geom_density(
    data = . %>%
      select(from, to, expected.edges) %>%
      distinct() %>%
      mutate(rand.p = map(expected.edges, ~rpois(niter, .x))) %>%
      unnest(),
    aes(x = rand.p, y = stat(count)), adjust = 10, colour = "green"
  ) +
  geom_vline(
    data = . %>%
      group_by(from, to) %>%
      summarize(mean = mean(Freq)) %>%
      mutate(CI95h = map_dbl(mean, ~qchisq(0.975, 2 * (.x + 1)) / 2)),
    aes(xintercept = CI95h), colour = "blue"
  ) +
  geom_vline(
    data = . %>%
      group_by(from, to) %>%
      summarize(mean = mean(Freq)),
    aes(xintercept = mean), colour = "red"
  ) +
  geom_vline(
    data = . %>%
      select(from, to, weight) %>%
      distinct(),
    aes(xintercept = weight), colour = "purple"
  ) +
  scale_x_continuous(breaks = seq(0, 15, 2))
