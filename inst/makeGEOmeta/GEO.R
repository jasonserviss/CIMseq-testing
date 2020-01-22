library(CIMseq.data)
library(tidyverse)

keep.plates <- c(
  "NJA01201", "NJA01202", "NJA01203", "NJA01205", "NJA01301", 
  "NJA01302", "NJA01303", "NJA01401", "NJA01501", "NJA01503", "NJA01504", 
  "NJA01801", "NJA01803", "NJA01901", "NJA02001", "NJD00101", "NJD00102", 
  "NJD00103", "NJD00104"
)


GEOmeta <- MGA.Meta %>%
  filter(!filtered) %>%
  filter(unique_key %in% keep.plates) %>%
  {tibble(
  `sample name` = str_replace(pull(., sample), "..(.*)", "\\1"),
  title = pull(., sample),
  `source name` = pull(., sub_tissue),
  organism = "Mus musculus",
  `characteristics: strain` = pull(., subject_strain),
  `characteristics: tissue` = pull(., sub_tissue),
  `characteristics: age` = pull(., subject_age),
  molecule = "RNA",
  description = pull(., cellNumber),
  `processed data file` = "MGA.Counts.txt",
  `raw file` = NA
  )}

counts <- MGA.Counts %>%
  .[, pull(GEOmeta, title)]

ercc <- MGA.CountsERCC %>%
  .[, pull(GEOmeta, title)]

GEOcounts <- rbind(counts, ercc)

write_tsv(GEOmeta, '~/Github/CIMseq.testing/inst/makeGEOmeta/meta.txt')
write_tsv(as.data.frame(GEOcounts), '~/Github/CIMseq.testing/inst/makeGEOmeta/MGA.Counts.txt')

# 180810_BEA18P131 - NJA01201(?), NJA01202, NJA01203, NJA01205
# 180822_BEA - NJA01301, NJA01303
# 180828_BEA - NJA01302, NJA01401
# 181011_BEA18P142d - NJA01501, NJA01503
# 181030_BEA18P142g - NJA01201(?), NJA01504
# 181105_BEA18P191 - NJA01602, NJA01701
# 181116_BEA18P191e - NJD00101, NJD00103
# 181210_BEA18P191g - NJA01801
# 190107_BEA18P191h - NJA01803
# 190111_BEA19P004b - NJD00102, NJD00104
# 190228_BEA19P004e - NJA01901(?), NJA02001(?)
# 190306_BEA19P004e - NJA01901(?), NJA02001(?)