costLandscape <- function(singlets, multiplets, swarm, samples, nCellsInMultiplet) {
  singletIdx <- getData(swarm, "singletIdx")
  syntheticMultiplets <- t(appropriateSinglets(
    singlets, singletIdx, getData(multiplets, "features")
  ))
  combs <- combn(unique(getData(cObjSng, "classification")), nCellsInMultiplet)
  fractions <- t(sapply(1:ncol(combs), function(i) {
    u.classes <- unique(getData(singlets, "classification"))
    frac <- vector(mode = "numeric", length = length(u.classes))
    names(frac) <- sort(u.classes)
    frac[names(frac) %in% combs[, i]] <- 1 / nrow(combs)
    frac
  }))
  #rownames(fractions) <- paste(combs[1, ], combs[2, ], sep = "_")
  rownames(fractions) <- apply(combs, 2, function(c) paste(sort(c), collapse = "_"))
  
  mulCPM <- getData(multiplets, "counts.cpm")
  selectInd <- getData(multiplets, "features")
  mul <- mulCPM[selectInd, , drop = FALSE]
  
  calc.costs <- map_dfr(samples, function(s) {
    oneMultiplet <- round(mul[, s])
    costs <- sapply(1:nrow(fractions), function(i) {
      calculateCost(
        oneMultiplet, syntheticMultiplets, fractions[i, ], length(singletIdx)
      )
    })
    tibble(sample = s, cost = costs, comb = rownames(fractions))
  })
}

#step = 0.001, range = 0.1
costLandscape <- function(singlets, multiplets, swarm, sample, step, range) {
  singletIdx <- getData(swarm, "singletIdx")
  syntheticMultiplets <- t(appropriateSinglets(
    singlets, singletIdx, getData(multiplets, "features")
  ))
  adj <- adjustFractions(singlets, multiplets, swarm, theoretical.max = 4)
  adj <- adj[sample, ]
  fractions <- getData(swarm, "fractions")[sample, ]
  new.fractions <- map(1:length(which(adj == 1)), function(i) {
    idx <- which(adj == 1)[i]
    f <- fractions[names(adj)[idx]]
    target <- seq((f - range), (f + range), step)
    out <- do.call("rbind", replicate(length(target), fractions, simplify = FALSE))
    out[, idx] <- target
    out * 1 / rowSums(out)
  }) %>%
    do.call("rbind", .)
  
  rownames(new.fractions) <- paste(
    rep(
      names(adj)[adj == 1], 
      each = length(seq(-range, range, step))
    ),
    seq(-range, range, step),
    sep = "_"
  )
  
  mulCPM <- getData(multiplets, "counts.cpm")
  selectInd <- getData(multiplets, "features")
  mul <- mulCPM[selectInd, , drop = FALSE]
  oneMultiplet <- round(mul[, sample])
  future::plan(future::multiprocess)
  options(future.globals.maxSize = Inf)
  costs <- future.apply::future_lapply(1:nrow(new.fractions), FUN = function(i) {
    calculateCost(
      oneMultiplet, syntheticMultiplets, new.fractions[i, ], length(singletIdx)
    )
  })
  costs <- as.numeric(costs)
  tibble(cost = costs, name = rownames(new.fractions)) %>%
    separate(name, into = c("class", "change"), sep = "_") %>%
    mutate(change = as.numeric(change)) %>%
    ggplot() +
    geom_point(aes(change, cost, colour = class), size = 0.25) +
    facet_wrap(~class, scales = "free_x") +
    theme_bw() +
    guides(colour = FALSE)
}
