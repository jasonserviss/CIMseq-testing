packages <- c("CIMseq", "CIMseq.data", "tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

#1) run the deconvolution using 5, 10, 100, 500, 1000, 5000, 10000 synthetic 
# multiplets for 6 different samples.
#2) for each number of synthetic multiplets run the deconvolution 5 times 
# providing a different set of synthetic multiplets each time.
#3) For each sample and each number of synthetic multiplets, calculate the
# intersect (%, the number of intersecting solutions of all the reported 
# solutions in 2 replicates * 100) of the identified edges when comparing all 
# replicates to each other.

load('~/Github/CIMseq.testing/inst/analysis/MGA.analysis_engeSimple/data/CIMseqData.rda')
s <- c(
  "m.NJA01201.E21", "m.NJA01202.B19", "m.NJA01301.H19", "m.NJA01301.H21", 
  "m.NJA01301.O23", "m.NJA01302.D21"
)
cObjMul.2 <- CIMseqMultiplets(
  getData(cObjMul, "counts")[, s], getData(cObjMul, "counts.ercc")[, s],
  getData(cObjMul, "features")
)

si <- map(1:8, ~swarmInit(cObjSng, 2)) %>% do.call(cbind, .)
maxiter <- 100
eps.stagnate <- 1
maxit.stagnate <- 5
tests <- c(5, 10, 100, 500, 1000, 5000, 10000)


run <- function(nsm, reps) {
  base.seed <- 76238945
  sidx <- map(reps, function(i) {
    set.seed(base.seed + i)
    purrr::map(1:nsm, ~sampleSinglets(getData(cObjSng, "classification")))
  })
  
  plan(multiprocess, workers = 6)
  map(reps, function(i) {
    CIMseqSwarm(
      cObjSng, cObjMul.2, maxiter = maxiter, swarmsize = ncol(si), 
      nSyntheticMultiplets = nsm, swarmInit = si, singletIdx = sidx[[i]],
      psoControl = list(
        eps.stagnate = eps.stagnate, maxit.stagnate = maxit.stagnate
      )
    )
  })
}

res <- map(tests, ~run(.x, 1:5))

p <- map_dfr(1:length(res), function(i) {
  t <- res[[i]]
  edges <- map_dfr(
    t, ~getEdgesForMultiplet(.x, cObjSng, cObjMul.2), 
    .id = "rep.nr"
  ) %>%
    mutate(
      conn = map2_chr(from, to, function(f, t) {
        paste(sort(c(f, t)), collapse = "_")
      })) %>%
    select(-from, -to) %>%
    distinct()
  combs <- combn(reps, 2)
  s <- unique(pull(edges, sample))
  m.int <- map(s, function(samp) {
    int <- map(1:ncol(combs), function(i) {
      a <- filter(edges, sample == samp & rep.nr == combs[1, i])$conn
      b <- filter(edges, sample == samp & rep.nr == combs[2, i])$conn
      length(intersect(a, b)) / length(unique(filter(edges, sample == samp)$conn))
    })
    as.numeric(int)
  })
  tibble(sample = s, intersect = m.int, test = tests[i])
})

plot <- p %>% 
  unnest() %>%
  ggplot() + 
  geom_jitter(
    aes(factor(test), intersect * 100, colour = sample), width = 0.2, height = 0
  ) +
  labs(x = "Synthetic multiplets (n)", y = "Intersect (%)") + 
  theme_bw() +
  facet_wrap(~sample) +
  guides(colour = FALSE)


save(res, file = "~/Desktop/nrSynthMulVsIntersect.rda")
ggsave(plot, file = "~/Desktop/nrSynthMulVsIntersect.pdf")