library(Seurat)
library(tidyverse)
library(future.apply)
load('data/seuratObj.rda')
options("future.globals.maxSize" = Inf)
plan(multiprocess)
classes <- unique(mca@meta.data$res.1)
de <- future.apply::future_lapply(classes, FUN = function(c) {
  FindMarkers(mca, ident.1 = c, only.pos = TRUE, min.diff.pct = 0.2, logfc.threshold = log(2), test.use = "roc")
}) %>%
  map_dfr(function(x) rownames_to_column(x, "gene"), .id = "idx") %>%
  as_tibble() %>%
  mutate(idx = as.integer(idx)) %>%
  mutate(cluster = classes[idx]) %>%
  select(-idx)


save(de, file = "data/DE.rda")