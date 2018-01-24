
addComplexity <- function(uObj, k = 6, plot = TRUE) {

  #add complexity
  tsne <- getData(uObj, "tsne")
  km <- kmeans(tsne, k, iter.max = 100, nstart = 100)[[1]]
  km <- tibble(class = km, multiplet = names(km))
  
  #plot
  p <- c()
  if(plot) {
    p <- tsne %>%
    as.data.frame() %>%
    rownames_to_column(var = "multiplet") %>%
    as_tibble() %>%
    full_join(km, by = "multiplet") %>%
    mutate(class = parse_factor(class, levels = unique(class))) %>%
    ggplot() +
    geom_point(aes(x = V1, y = V2, colour = class))
  }
  
  #rename and add to spUnsupervised object
  newClass <- tibble(
    multiplet = rownames(tsne),
    classReal = getData(uObj, "classification")
  ) %>%
  full_join(km, by = "multiplet") %>%
  unite(combined, classReal, class) %>%
  pull(combined)
  
  classification(uObj) <- newClass
  return(list(uObj, p))
}
