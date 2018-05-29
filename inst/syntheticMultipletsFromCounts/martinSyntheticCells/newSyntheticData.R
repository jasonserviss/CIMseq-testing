#combos: a matrix the combinations of cell types to be included in the 
# multiplets. Typically generated with combn.
#adjustment: a numeric vector with length = the number of unique cell types
#indicating the fraction of contribution for each cell type. Names should 
#indicate the corresponding cell type.
#seed: seed for random number initiation. When generating many multiplets by 
# running the function multiple times in a loop, the seed should be set 
# dynamically in the function calling the makeSyntheticData function to avoid 
# all multiplets being identical.
makeSyntheticData <- function(
  counts,
  classes,
  combos,
  adjustment,
  seed = 87909023
){
  #generates one multiplet per combo
  set.seed(seed)
  sums <- t(rowsum(t(counts), classes))
  output <- data.frame(gene = rownames(counts), stringsAsFactors = FALSE)
  
  combos %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    map(., function(x) {
      adj <- adjustment[match(x, names(adjustment))]
      adjusted <- t(t(sums[, colnames(sums) %in% x]) * adj)
      mixed <- rowSums(adjusted)
      
      sample(
        x = names(mixed), replace = TRUE, prob = mixed,
        size = sum(mixed) / nrow(combos) / sum(classes %in% x)
      ) %>%
        table() %>%
        as.data.frame(., row.names = names(.), stringsAsFactors = FALSE) %>%
        setNames(c("gene", paste(sort(x), collapse = "-")))
    }) %>%
    reduce(full_join, by = "gene") %>%
    full_join(output, by = "gene") %>%
    column_to_rownames("gene") %>% 
    replace(is.na(.), 0)
}

makeSyntheticData <- function(
  counts,
  classes,
  combos,
  adjustment,
  seed = 87909023
){
  set.seed(seed)
  sums <- t(rowsum(t(counts), classes, reorder = FALSE))
  output <- data.frame(gene = rownames(counts), stringsAsFactors = FALSE)
  
  combos %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    map(., function(x) {
      adj <- adjustment[match(x, names(adjustment))]
      adjusted <- t(t(sums[, colnames(sums) %in% x]) * adj)
      rs <- rowSums(adjusted)
      
      sample(
        x = names(rs), replace = TRUE, prob = rs,
        size = sum(rs) / sum(classes %in% x)
      ) %>%
        table() %>%
        as.data.frame(., row.names = names(.), stringsAsFactors = FALSE) %>%
        setNames(c("gene", paste(sort(x), collapse = "-")))
    }) %>%
    reduce(full_join, by = "gene") %>%
    full_join(output, by = "gene") %>%
    replace(is.na(.), 0) %>%
    as_tibble()
}

makeSyntheticData2 <- function(
  counts,
  classes,
  combos,
  adjustment,
  seed = 87909023
){
  set.seed(seed)
  output <- data.frame(gene = rownames(counts), stringsAsFactors = FALSE)
  classCounts <- split(counts, classes)
  classCounts <- lapply(classCounts, function(x) matrix(x, nrow = nrow(counts)))
  
  combos %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    map(., function(x) {
      c <- classCounts[x]
      adj <- adjustment[match(names(c), names(adjustment))]
      adjusted <- map2(c, adj, function(x, y) x * y)
      
      map_dbl(1:nrow(counts), function(y) 
        map_dbl(adjusted, function(z) z[y, sample(1:ncol(z), 1)]) %>%
          mean()
      ) %>%
        as_tibble() %>%
        setNames(paste(sort(x), collapse = "-")) %>%
        add_column(gene = rownames(counts))
    }) %>%
    reduce(full_join, by = "gene") %>%
    full_join(output, by = "gene") %>%
    replace(is.na(.), 0) %>%
    as_tibble()
}

makeSyntheticData3 <- function(
  counts,
  classes,
  combos,
  adjustment,
  seed = 87909023
){
  set.seed(seed)
  sums <- t(rowsum(t(counts), classes, reorder = FALSE))
  output <- data.frame(gene = rownames(counts), stringsAsFactors = FALSE)
  
  combos %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    map(., function(x) {
      adj <- adjustment[match(x, names(adjustment))]
      adjusted <- t(t(sums[, colnames(sums) %in% x]) * adj)
      rs <- rowSums(adjusted)
      rm <- rowMeans(adjusted)
      
      pois <- rpois(length(rs), rm)
      #nbin <- rnbinom(length(rs), mu = runif(length(rs), 0, rs), size = 1)
      #nbin <- rnbinom(length(rmin), mu = runif(length(rmin), rmin, rmax), size = 1)
      
      sample(
        x = rownames(counts), replace = TRUE, prob = nbin, 
        size = sum(rs) / sum(classes %in% x)
      ) %>%
        table() %>%
        as.data.frame(., row.names = names(.), stringsAsFactors = FALSE) %>%
        setNames(c("gene", paste(sort(x), collapse = "-")))
    }) %>%
    reduce(full_join, by = "gene") %>%
    full_join(output, by = "gene") %>%
    replace(is.na(.), 0) %>%
    as_tibble()
}

makeSyntheticData2 <- function(
counts,
classes,
combos,
adjustment,
noise = FALSE,
seed = 87909023
){
  set.seed(seed)
  
  #adjust counts by adjustment vector
  adj <- rep(adjustment, length(classes) / length(adjustment))
  adjusted <- t(t(counts) * adj[match(classes, names(adj))])
  
  #calculate number of dropouts and non-dropouts per class per gene
  dos <- countDropOuts(counts, classes, FUN = function(y) {length(which(y == 0))})
  ndos <- countDropOuts(counts, classes, FUN = function(y) {length(which(y != 0))})
  
  #calculate mean and sd gene expression per gene per class
  ms <- means(adjusted, classes)
  ss <- sds(adjusted, classes)
  
  #setup output structure
  output <- data.frame(gene = rownames(counts), stringsAsFactors = FALSE)
  
  #generate synthetic multiplet
  combos %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  map(., function(x) {
    map_dbl(1:nrow(counts), function(y) {
      if(!is.nan(ms[y, x[1]]) & !is.na(ss[y, x[1]])){
        n1 <- rnorm(ndos[y, x[1]], ms[y, x[1]], ss[y, x[1]])
      } else {
        n1 <- NA
      }
      
      if(!is.nan(ms[y, x[2]]) & !is.na(ss[y, x[2]])) {
        n2 <- rnorm(ndos[y, x[2]], ms[y, x[2]], ss[y, x[2]])
      } else {
        n2 <- NA
      }
      
      dos1 <- rep(0, dos[y, x[1]])
      dos2 <- rep(0, dos[y, x[2]])
      
      s1 <- sample(c(n1, dos1)[!is.na(c(n1, dos1))], 1)
      s2 <- sample(c(n2, dos2)[!is.na(c(n2, dos2))], 1)
      
      if_else(mean(s1, s2) < 0, 0, mean(s1, s2))
    }) %>%
    data.frame(gene = rownames(counts), ., stringsAsFactors = FALSE) %>%
    setNames(c("gene", paste(sort(x), collapse = "-")))
  }) %>%
  reduce(full_join, by = "gene") %>%
  full_join(output, by = "gene") %>%
  replace(is.na(.), 0) %>%
  as_tibble()
}

countDropOuts <- function(counts, classes, FUN = function(y) {length(which(y == 0))}) {
  counts %>%
  t() %>%
  split(classes) %>%
  map(function(x) {
    m <- matrix(x, nrow = nrow(counts), byrow = TRUE)
    dos <- apply(m, 1, FUN)
    tibble(gene = rownames(counts), dos)
  }) %>%
  reduce(full_join, by = "gene") %>%
  as.data.frame() %>%
  column_to_rownames("gene") %>%
  setNames(sort(unique(classes)))
}

means <- function(counts, classes) {
  counts %>%
  t() %>%
  split(classes) %>%
  map(function(x) {
    m <- matrix(x, nrow = nrow(counts), byrow = TRUE)
    means <- apply(m, 1, function(z) mean(z[z != 0]))
    tibble(gene = rownames(counts), means)
  }) %>%
  reduce(full_join, by = "gene") %>%
  as.data.frame() %>%
  column_to_rownames("gene") %>%
  setNames(sort(unique(classes)))
}

sds <- function(counts, classes) {
  counts %>%
  t() %>%
  split(classes) %>%
  map(function(x) {
    m <- matrix(x, nrow = nrow(counts), byrow = TRUE)
    sds <- apply(m, 1, function(z) sd(z[z != 0])) / ncol(m)
    tibble(gene = rownames(counts), sds)
  }) %>%
  reduce(full_join, by = "gene") %>%
  as.data.frame() %>%
  column_to_rownames("gene") %>%
  setNames(sort(unique(classes)))
}


makeSyntheticData5 <- function(
counts,
classes,
combos,
adjustment,
seed = 87909023
){
  set.seed(seed)
  output <- data.frame(gene = rownames(counts), stringsAsFactors = FALSE)
  
  combos %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  map(., function(x) {
    exCounts <- counts[, map_int(x, ~sample(which(classes == .x), 1))]
    colnames(exCounts) <- x
    
    adj <- adjustment[match(x, names(adjustment))]
    adjusted <- t(t(exCounts) * adj)
    
    rs <- rowSums(adjusted)
    expanded <- rep(names(rs), rs)
    
    sample(
    x = expanded, replace = FALSE,
    size = length(expanded) / sum(classes %in% x)
    ) %>%
    table() %>%
    as.data.frame(., row.names = names(.), stringsAsFactors = FALSE) %>%
    setNames(c("gene", paste(sort(x), collapse = "-")))
  }) %>%
  reduce(full_join, by = "gene") %>%
  full_join(output, by = "gene") %>%
  replace(is.na(.), 0) %>%
  as_tibble()
}
