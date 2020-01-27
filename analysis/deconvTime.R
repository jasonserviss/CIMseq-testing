#function that estimates uppmax CIMseq deconvolution times

deconvTime <- function(dir = file.path(getwd(), "logs")) {
  f <- list.files(dir, pattern = "CIMseqSwarm.out", recursive = TRUE, full.names = TRUE)
  lines <- readLines(f)
  idx <- purrr::map_chr(lines, function(l) {
    dplyr::case_when(
      !is.na(.start(l)) ~ "start",
      !is.na(.end(l)) ~ "end",
      TRUE ~ "empty"
    )
  })
  
  start <- unlist(purrr::map(lines[idx == "start"], .start))
  finish <- unlist(purrr::map(lines[idx == "end"], .end))
  diff <- finish - start
  tibble(
    min = lubridate::seconds_to_period(min(diff)),
    mean = lubridate::seconds_to_period(mean(diff)),
    median = lubridate::seconds_to_period(median(diff)),
    max = lubridate::seconds_to_period(max(diff))
  )
}

.start <- function(string) {
  if(stringr::str_detect(string, "Starting deconvolution at ")) {
    chr <- stringr::str_replace(string, ".*Starting deconvolution at (.*)", "\\1")
    lubridate::ymd_hms(chr)
  } else {
    NA
  }
}

.end <- function(string) {
  if(stringr::str_detect(string, "Finished deconvolution at ")) {
    chr <- stringr::str_replace(string, ".*Finished deconvolution at (.*)", "\\1")
    lubridate::ymd_hms(chr)
  } else {
    NA
  }
}

.sample <- function(string) {
  if(stringr::str_detect(string, "Running ")) {
    stringr::str_replace(string, ".*Running (.*).", "\\1")
  } else {
    NA
  }
}
