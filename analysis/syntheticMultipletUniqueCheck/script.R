#PACKAGES
library(tidyverse)
library(CIMseq)
library(CIMseq.data)

#DATA
load('~/Github/CIMseq.testing/inst/analysis/SCM.analysis/data/CIMseqData.rda')
classes <- getData(cObjSng, "classification")

set.seed(8329)
d <- tibble(nSyntheticMultiplets = seq(2, 10000, 500)) %>%
  mutate(sm = map(nSyntheticMultiplets, function(n) {
    purrr::map(1:n, ~sampleSinglets(classes))
  })) %>%
  mutate(unique.idx = map_int(sm, ~length(unique(.x)))) %>%
  mutate(sng = map(sm, ~appropriateSinglets(cObjSng, .x, 1))) %>%
  mutate(unique.sng = map_int(sng, function(s) {
    s %>%
      apply(., 1, sort) %>%
      t() %>%
      matrix_to_tibble(drop = TRUE) %>%
      distinct() %>%
      nrow()
  }))

d %>%
  select(nSyntheticMultiplets, unique.idx, unique.sng) %>%
  gather(type, value, -nSyntheticMultiplets) %>%
  ggplot() +
  geom_point(aes(nSyntheticMultiplets, value, colour = type), size = 1) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Synthetic multiplets generated (n)", y = "n") +
  guides(colour = guide_legend(title = ""))
