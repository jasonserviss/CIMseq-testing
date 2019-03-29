library(CIMseq.data)
library(tidyverse)
library(RNAseqFunctions)

cpm <- RNAseqFunctions::cpm(MGAB.Counts)

#plot function
plotHeatSI <- function(data, nColors) {
  p <- data %>%
    inner_join(MGAB.Meta) %>%
    filter(sub_tissue == "small_intestine") %>%
    mutate(replicate = if_else(
      str_detect(sample, "NJE003"), 
      "Replicate 1", "Replicate 2"
    )) %>%
    mutate(xname = case_when(
      str_detect(sample, "S01") ~ "Section 1",
      str_detect(sample, "S02") ~ "Section 2",
      str_detect(sample, "S03") ~ "Section 3",
      str_detect(sample, "S04") ~ "Section 4",
      str_detect(sample, "S05") ~ "Section 5",
      str_detect(sample, "S06") ~ "Section 6",
      str_detect(sample, "S07") ~ "Section 7",
      str_detect(sample, "S08") ~ "Section 8"
    )) %>%
    mutate(xname = parse_factor(xname, levels = paste0("Section ", 1:8))) %>%
    group_by(HGN) %>%
    arrange(cpm, .by_group = TRUE) %>%
    mutate(colour = viridis::viridis(nColors)) %>%
    ungroup() %>%
    arrange(type) %>%
    mutate(HGN = parse_factor(HGN, levels = unique(HGN))) %>%
    ggplot() + 
    geom_tile(aes(xname, HGN, fill = colour)) + 
    facet_wrap(~replicate, scales = "free") +
    scale_fill_identity() +
    labs(
      y = paste0("Genes (", length(unique(pull(data, HGN))), ")"), 
      x = "proximal -> distal", title = "Small intestine"
    ) +
    theme_bw()
  
  return(p)
}

#find genes
pairs <- list(
  c("NJE003S01", "NJE004S01"), c("NJE003S02", "NJE004S02"),
  c("NJE003S03", "NJE004S03"), c("NJE003S04", "NJE004S04"),
  c("NJE003S05", "NJE004S05"), c("NJE003S06", "NJE004S06"),
  c("NJE003S07", "NJE004S07"), c("NJE003S08", "NJE004S08")
)

means <- map(pairs, ~matrixStats::rowMeans2(cpm[, .x]))
means <- matrix(unlist(means), ncol = length(pairs))
colnames(means) <- paste0("Section ", 1:length(pairs))
rownames(means) <- rownames(cpm)
corr <- cor(t(means[rowSums(means) != 0, ]), 1:length(pairs))
corr <- corr[order(corr[,1]), ]

#plot
n <- 1000
negative <- cpm[rownames(cpm) %in% head(names(corr), n = n), ] %>% 
  matrix_to_tibble("HGN") %>% 
  gather(sample, cpm, -HGN) %>%
  mutate(type = "negative")

positive <- cpm[rownames(cpm) %in% tail(names(corr), n = n), ] %>% 
  matrix_to_tibble("HGN") %>% 
  gather(sample, cpm, -HGN) %>%
  mutate(type = "positive")

plotHeatSI(bind_rows(negative, positive), 16) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

#save some data
selected <- corr[corr > 0.9 | corr < -0.9]
dataToSave <- tibble(gene = names(selected), correlation = selected)
#write.table(dataToSave, file = "", sep = "\t", quote = FALSE, row.names = FALSE)