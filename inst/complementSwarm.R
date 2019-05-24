complementSwarm <- function(seedFile, outputDir, newSeedFile) {
  expected.samples <- read.table(seedFile, stringsAsFactors = FALSE)[[1]]
  detected.files <- list.files(outputDir, pattern = "*.rda")
  detected.samples <- gsub("sObj_(.*)_uppmax.rda", "\\1", detected.files)
  if(length(detected.samples) == length(expected.samples)) {
    message('All expected samples are present.')
    return(NULL)
  } else {
    toAdd <- expected.samples[!expected.samples %in% detected.samples]
    write.table(
      data.frame(X1 = toAdd, stringsAsFactors = FALSE), file = newSeedFile,
      col.names = FALSE, row.names = FALSE, quote = FALSE
    )
  }
}