packages <- c("sp.scRNAseq", "sp.scRNAseqData", "ggthemes", "tidyverse", "seqTools")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

#setup spCounts
s <- str_detect(colnames(countsMgfp), "^s")
countsMgfpMeta <- countsMgfpMeta %>%
  mutate(mouseName = if_else(
  str_detect(plate, "NJA004"),
  "C57black/6J", "C57black/6J/Lgr5gfp"
  ))
m <- colnames(countsMgfp) %in% filter(countsMgfpMeta, mouseName == "C57black/6J")$sample

cObjSng <- spCounts(countsMgfp[, s & m], cbind(countsMgfpERCC[, s & m]))
cObjMul <- spCounts(countsMgfp[, !s & m], cbind(countsMgfpERCC[, !s & m]))

#spUnsupervised
uObj <- spUnsupervised(cObjSng, max_iter = 7000, initial_dims = sum(s), max = 250, perplexity = 7, seed = 73223)
save(uObj, file = "./uObj.rda")

#Rename the classes
classifications <- tibble(
  sample = rownames(getData(uObj, "tsne")),
  oldClass = getData(uObj, "classification")
) %>%
  mutate(newClass = case_when(
    oldClass == "C1" ~ "C.Stem",
    oldClass == "F1" ~ "C.Colonocyte",
    oldClass == "B1" ~ "C.Goblet",
    oldClass == "G1" ~ "Blood.Tufft.Endocrine",
    oldClass == "D1" ~ "SI.Stem",
    oldClass == "A1" ~ "SI.Enterocyte",
    oldClass == "E1" ~ "SI.Goblet",
    TRUE ~ "error"
  ))

classification(uObj) <- classifications$newClass[match(rownames(getData(uObj, "tsne")), classifications$sample)]
groupMeans(uObj) <- averageGroupExpression(cObjSng, getData(uObj, "classification"), FALSE)
tsneMeans(uObj) <- tsneGroupMeans(getData(uObj, "tsne"), getData(uObj, "classification"))

#KS test strategy 1
ks.res <- KStest(getData(cObjSng, "counts.cpm"), getData(uObj, "classification"), cores = 8)
p.ks.res <- processKStest(ks.res, getData(uObj, "classification"), 0.05) %>%
  filter(sigBool) %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(p.prod = map_dbl(p.values, prod)) %>%
  mutate(meanCpm = map2_dbl(id, gene, function(i, g) {
    mean(getData(cObjSng, "counts.cpm")[g, getData(uObj, "classification") == i])
  })) %>%
  mutate(log2.fc = map2_dbl(id, gene, function(i, g) {
    log2(
    mean(getData(cObjSng, "counts.cpm")[g, getData(uObj, "classification") == i]) /
    mean(getData(cObjSng, "counts.cpm")[g, getData(uObj, "classification") != i])
    )
  })) %>%
  arrange(desc(log2.fc), p.prod, desc(statSum), desc(meanCpm))

save(ks.res, file = './differentialExpressionStrategy.ks.res1.rda', compress = "bzip2")
save(p.ks.res, file = './differentialExpressionStrategy1.rda', compress = "bzip2")
openxlsx::write.xlsx(p.ks.res, file = './differentialExpressionStrategy1.xlsx')

#KS test strategy 2
idx <- which(getData(uObj, "classification") %in% c("C.Stem", "SI.Stem"))
counts <- getData(cObjSng, "counts.cpm")[, idx]
classes <- getData(uObj, "classification")[idx]

#find top 50 cells in each class with highest Lgr5
top50 <- tibble(Lgr5 = counts["Lgr5", ], class = classes, sample = colnames(counts)) %>%
group_by(class) %>%
top_n(50, Lgr5) %>%
ungroup()

ks.res2 <- KStest(counts[, match(top50$sample, colnames(counts))], top50$class, cores = 8)
p.ks.res2 <- ks.res2 %>%
filter(p.value < 0.001) %>%
separate(combination, c("Celltype1", "Celltype2"), sep = "-")  %>%
mutate(
meanCpmCelltype1 = map2_dbl(Celltype1, gene, function(i, g) {
  mean(getData(cObjSng, "counts.cpm")[g, getData(uObj, "classification") == i])
}),
meanCpmCelltype2 = map2_dbl(Celltype2, gene, function(i, g) {
  mean(getData(cObjSng, "counts.cpm")[g, getData(uObj, "classification") == i])
})
) %>%
mutate(log2.fc = log2(meanCpmCelltype1 / meanCpmCelltype2)) %>%
arrange(desc(log2.fc), p.value)

save(ks.res2, file = './differentialExpressionStrategy.ks.res2.rda', compress = "bzip2")
save(p.ks.res2, file = './differentialExpressionStrategy2.rda', compress = "bzip2")
openxlsx::write.xlsx(p.ks.res2, file = './differentialExpressionStrategy2.xlsx')
