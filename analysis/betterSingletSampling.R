averageGroupExpression <- function(singlets) {
  classes <- getData(singlets, "classification")
  u.classes <- unique(classes)
  exp <- getData(singlets, "counts.cpm")
  
  sapply(u.classes, function(x) {
    matrixStats::rowMeans2(exp[, classes == x])
  })
}

medianGroupExpression <- function(singlets) {
  classes <- getData(singlets, "classification")
  u.classes <- unique(classes)
  exp <- getData(singlets, "counts.cpm")
  
  sapply(u.classes, function(x) {
    matrixStats::rowMedians(exp[, classes == x])
  })
}

distToCenter <- function(singlets, method) {
  if(method == "mean") {
    center <- averageGroupExpression(singlets)
  } else if(method == "median") {
    center <- medianGroupExpression(singlets)
  } else {
    stop("Enter a valid method.")
  }
  
  cpm <- getData(singlets, "counts.cpm")
  classes <- getData(singlets, "classification")
  
  dists <- sapply(1:ncol(cpm), function(i) {
    d <- as.matrix(dist(t(cbind(cpm[, i], center))))[, 1]
    unname(d[classes[i]])
  })
  tibble(
    sample = colnames(cpm),
    class = classes,
    distToCenter = dists
  )
}

pickBestClassRepresentation <- function(
  singlets, method, topPercent, minPerGroup = NULL, safe = TRUE
){
  filtered <- distToCenter(singlets, method) %>%
    group_by(class) %>%
    arrange(class, distToCenter) %>%
    filter(distToCenter < quantile(distToCenter, topPercent)) 
  
  if(safe & !is.null(minPerGroup)) {
    n <- filtered %>% summarize(n = n()) %>% pull(n)
    if(any(n < minPerGroup)) stop("Too few in group. Increase top percent.")
  }
  
  pull(filtered, sample)
}  