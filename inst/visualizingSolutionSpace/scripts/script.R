packages <- c(
  "sp.scRNAseqData", "sp.scRNAseq", "seqTools",
  "printr", "ggthemes", "tidyverse", "future"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currentDir <- getwd()

s <- str_detect(colnames(countsSorted2), "^s")
sng <- countsSorted2[, s]
mul <- countsSorted2[, !s]
cObjSng <- spCounts(sng, countsSortedERCC2[, s])
bool <- colnames(mul) %in% c("m.NJB00204.G04", "m.NJB00204.D07", "m.NJB00204.F12")
cObjMul <- spCounts(mul[, bool], countsSortedERCC2[, !s][, bool])
uObj <- spUnsupervised(cObjSng)
plan(multiprocess)
sObj <- spSwarm(
  cObjSng, cObjMul, uObj, maxiter = 40, swarmsize = 200, 
  nSyntheticMultiplets = 500, report = TRUE, reportRate = 1
)

ss <- getData(sObj, "spSwarm")
colnames(ss) <- c("HOS", "HCT116", "A375")
sObj@spSwarm <- ss

save(cObjSng, cObjMul, uObj, sObj, file = file.path(currentDir, 'data/algoOutput.rda'))

processReports <- function(sObj) {
  tibble(
    sample = rownames(getData(sObj, "spSwarm")),
    iteration = map(getData(sObj, "stats"), function(x) x$it),
    error = map(getData(sObj, "stats"), function(x) x$error),
    fitness = map(getData(sObj, "stats"), function(x) x$f),
    position = map(getData(sObj, "stats"), function(x) {
      map(x$x, function(y) t(y) * 1/colSums(y))
    })
  ) %>% 
    unnest() %>%
    mutate(position = map(position, function(x) {
      x %>% 
        as.data.frame() %>% 
        setNames(colnames(getData(sObj, "spSwarm"))) %>% 
        as_tibble() %>%
        add_column(swarmMemberID = 1:nrow(.), .before = 1) %>%
        ungroup() 
    }))
}

report <- processReports(sObj)

save(report, file = file.path(currentDir, "data/report.rda"))
writeLines(capture.output(sessionInfo()), "logs/sessionInfo.txt")

