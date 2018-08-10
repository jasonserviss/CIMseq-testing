
#FUNCTIONS
determineEdgeCutoff <- function(cObjSng, uObj, sObj, reps = 10) {
  frac <- .determineEdgeCutoff_generateFractions_doublet(uObj, reps = 10)
  sm <- .determineEdgeCutoff_generateSM(cObjSng, uObj, sObj, frac, seed = 8923)
  swarm <- .determineEdgeCutoff_swarm(cObjSng, uObj, sm)
  .determineEdgeCutoff_process(swarm)
}

.determineEdgeCutoff_generateFractions_triplet <- function(uObj, reps = 10){
  classes <- getData(uObj, "classification")
  uClass <- sort(unique(classes))
  lc <- length(uClass)
  li <- length(c(1/lc, 10^-(1:5)))
  
  tibble(
    class_reduced = rep(uClass, each = li * reps),
    interval = rep(rep(round(c(1/lc, 10^-(1:5)), digits = 5), each = reps), lc),
    reps = rep(1:reps, lc * li)
  ) %>%
    mutate(fractions = map2(class_reduced, interval, function(c, i) {
      fractions <- rep(0, lc)
      names(fractions) <- uClass
      fractions[names(fractions) == c] <- i
      rest <- (1 - sum(fractions)) / (length(fractions) - 1)
      fractions[names(fractions) != c] <- rest
      fractions
    }))
}

.determineEdgeCutoff_generateFractions_doublet <- function(uObj, reps = 10){
  classes <- getData(uObj, "classification")
  uClass <- sort(unique(classes))
  key <- map(uClass, function(c) {
    combos <- combn(uClass, 2)
    x <- combos[, combos[1,] == c | combos[2,] == c]
    c(x)[c(x) != c]
  })
  names(key) <- uClass
  
  li <- length(c(0.5, 10^-(1:5)))
  lk <- length(key[[1]])
  lc <- length(uClass)
  
  tibble(
    class_reduced = rep(uClass, each = li * reps * lk),
    interval = rep(rep(round(c(0.5, 10^-(1:5)), digits = 5), each = reps * lk), lc),
    rep = rep(1:reps, lc * li * lk),
    cidx = rep(rep(1:lk, each = reps), (li * lk * lc) / 2)
  ) %>%
    mutate(combo = map2_chr(class_reduced, cidx, function(class, i) {
      key[[class]][i]
    })) %>%
    mutate(fractions = pmap(list(class_reduced, combo, interval), function(c1, c2, i) {
      fractions <- rep(0, length(uClass))
      names(fractions) <- uClass
      fractions[names(fractions) == c1] <- i
      fractions[names(fractions) == c2] <- 1 - i
      fractions[!names(fractions) %in% c(c1, c2)] <- 0
      fractions
    })) %>%
    select(-cidx)
}

.determineEdgeCutoff_generateSM <- function(cObjSng, uObj, sObj, fractions, seed = 8923) {
  selectInd <- getData(sObj, "arguments")$selectInd
  counts <- getData(cObjSng, "counts.cpm")[selectInd, ]
  classes <- getData(uObj, "classification")
  
  sm <- map(1:nrow(fractions), function(i) {
    eps <- sample(1:1000, 1)
    f <- as.numeric(pull(slice(fractions, i), fractions)[[1]])
    syntheticMultipletsFromCounts(counts, classes, f, seed + eps)
  }) %>% 
    bind_cols() %>% 
    as.matrix()
  
  cn <- select(fractions, -fractions) %>% unite(name, 1:ncol(.)) %>% pull(name)
  colnames(sm) <- cn
  rownames(sm) <- rownames(counts)
  spCounts(sm, matrix(NA))
}

.determineEdgeCutoff_swarm <- function(cObjSng, uObj, sm) {
  maxiter <- getData(sObj, "arguments")$maxiter
  swarmsize <- getData(sObj, "arguments")$swarmsize
  nSyntheticMultiplets <- getData(sObj, "arguments")$nSyntheticMultiplets
  selectInd <- getData(sObj, "arguments")$selectInd
  
  print(paste0("starting deconvolution: ", Sys.time()))
  plan(multiprocess)
  sobj <- spSwarm(
    cObjSng, sm, uObj, maxiter = maxiter, swarmsize = swarmsize, 
    nSyntheticMultiplets = nSyntheticMultiplets, selectInd = selectInd
  )
  print(paste0("finished deconvolution: ", Sys.time()))
  sobj
}

.determineEdgeCutoff_process <- function(uObj, sObj_out) {
  uClass <- unique(getData(uObj, "classification"))
  getData(sObj_out, "spSwarm") %>%
    rownames_to_column("info") %>%
    as_tibble() %>%
    separate(info, into = c("celltype_reduced", "interval", "repetition"), sep = "_") %>%
    nest(-(celltype_reduced:repetition)) %>%
    mutate(data = map(data, function(d) as.numeric(as.data.frame(d)[1, ]))) %>%
    mutate(mat = map(data, function(d) {
      matrix(d, nrow = 1, dimnames = list(NULL, uClass))
    })) %>%
    mutate(edges = map_chr(mat, function(m) {
      logic <- sp.scRNAseq:::.fractionCutoff(m, cutoff = 0)
      all.edges <- sp.scRNAseq:::.calculateWeight(m, logic)
      d <- filter(all.edges, weight > 0)
      paste(d$from, d$to, sep = "-", collapse = ", ")
    }))
}

syntheticMultipletsFromCounts <- function(
  counts, classes, fractions, seed = 87909023, ...
){
  set.seed(seed)
  out <- as.numeric(sampleSinglets(classes)) %>%
    subsetSinglets(counts, .) %>%
    adjustAccordingToFractions(fractions, .) %>%
    multipletSums()
  
  rownames(out) <- rownames(counts)
  out
}

#PACKAGES
packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse", "future", "future.apply")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

#DATA
s <- grepl("^s", colnames(countsSorted2))
cObjSng <- spCounts(countsSorted2[, s], countsSortedERCC2[, s])
load('../testingPoissonSorted/data/uObj.rda')
load('../testingPoissonSorted/data/sObj.rda')

#SCRIPT
frac <- .determineEdgeCutoff_generateFractions_doublet(uObj, reps = 10)
save(frac, file = "data/frac.rda", compress = "bzip2")

sm <- .determineEdgeCutoff_generateSM(cObjSng, uObj, sObj, frac, seed = 8923)
save(sm, file = "data/sm.rda", compress = "bzip2")

swarm <- .determineEdgeCutoff_swarm(cObjSng, uObj, sm)
save(swarm, file = "data/swarm.rda", compress = "bzip2")

edgeData <- .determineEdgeCutoff_process(swarm)

#SAVE
save(frac, sm, swarm, edgeData, file = "data/determineEdgeCutoff.rda", compress = "bzip2")
