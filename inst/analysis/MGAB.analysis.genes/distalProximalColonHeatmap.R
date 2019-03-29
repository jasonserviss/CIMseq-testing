library(tidyverse)
library(CIMseq)
library(CIMseq.data)

cpm <- CIMseq:::.norm.counts(MGAB.Counts)

#plot function
plotHeatColon <- function(data) {
  p <- data %>%
    inner_join(MGAB.Meta) %>%
    filter(sub_tissue == "colon") %>%
    mutate(replicate = if_else(
      str_detect(sample, "NJE003"), 
      "Replicate 1", "Replicate 2"
    )) %>%
    mutate(xname = case_when(
      str_detect(sample, "C01") ~ "Section 1",
      str_detect(sample, "C02") ~ "Section 2",
      str_detect(sample, "C03") ~ "Section 3",
      str_detect(sample, "C04") ~ "Section 4"
    )) %>%
    mutate(xname = parse_factor(xname, levels = paste0("Section ", 1:4))) %>%
    group_by(HGN) %>%
    arrange(cpm, .by_group = TRUE) %>%
    mutate(colour = viridis::viridis(8)) %>%
    ungroup() %>%
    arrange(type) %>%
    mutate(HGN = parse_factor(HGN, levels = unique(HGN))) %>%
    ggplot() + 
    geom_tile(aes(xname, HGN, fill = colour)) + 
    facet_wrap(~replicate, scales = "free") +
    scale_fill_identity() +
    labs(
      y = paste0("Genes (", length(unique(pull(data, HGN))), ")"), 
      x = "proximal -> distal", title = "Colon"
    ) +
    theme_bw()
  
  return(p)
}

plotLineColon <- function(data) {
  p <- data %>%
    inner_join(select(MGAB.Meta, sample, sub_tissue)) %>%
    filter(sub_tissue == "colon") %>%
    mutate(replicate = if_else(
      str_detect(sample, "NJE003"), 
      "Replicate 1", "Replicate 2"
    )) %>%
    mutate(xname = case_when(
      str_detect(sample, "C01") ~ "Section 1",
      str_detect(sample, "C02") ~ "Section 2",
      str_detect(sample, "C03") ~ "Section 3",
      str_detect(sample, "C04") ~ "Section 4"
    )) %>%
    mutate(xname = parse_factor(xname, levels = paste0("Section ", 1:4))) %>%
    group_by(HGN) %>%
    arrange(cpm, .by_group = TRUE) %>%
    mutate(colour = viridis::viridis(8)) %>%
    ungroup() %>%
    arrange(type) %>%
    mutate(HGN = parse_factor(HGN, levels = unique(HGN))) %>%
    ggplot() + 
    geom_line(aes(xname, cpm, colour = HGN, group = HGN)) + 
    facet_wrap(~replicate) +
    labs(y = "CPM", x = "proximal -> distal", title = "Colon"
    ) +
    theme_bw()
  
  return(p)
}

##Show "known" DE genes
#gene list from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3882682/
genes <- c(
  "Hoxb13", "Slc37a2", "Hoxd10", "Hoxd13", "Sval1", "Chst5", "Fxyd4", "B3gnt7",
  "Scgb1a1", "Tgm3", "Atp12a", "Lbp", "Best2", "Pla2g2f", "Sfrp2", "Ceacam10", 
  "Klk1c9", "Vsig1", "Ccdc153"
)

data <- cpm[rownames(cpm) %in% genes, ] %>% 
  matrix_to_tibble("HGN") %>% 
  gather(sample, cpm, -HGN) %>%
  mutate(type = NA)

plotHeatColon(data)

##Identify genes
pairs <- list(
  c("NJE003C01", "NJE004C01"), c("NJE003C02", "NJE004C02"),
  c("NJE003C03", "NJE004C03"), c("NJE003C04", "NJE004C04")
)

means <- map(pairs, ~matrixStats::rowMeans2(cpm[, .x]))
means <- matrix(unlist(means), ncol = 4)
colnames(means) <- paste0("Section ", 1:length(pairs))
rownames(means) <- rownames(cpm)
corr <- cor(t(means[rowSums(means) != 0, ]), 1:4)
corr <- corr[order(corr[,1]), ]

#Plot with replicates
n <- 1000
n.cor <- head(names(corr), n = n)
p.cor <- tail(names(corr), n = n)
colonSamples <- str_detect(colnames(cpm), "^NJE003C") | str_detect(colnames(cpm), "^NJE004C")
negative <- cpm[rownames(cpm) %in% n.cor, colonSamples] %>% 
  matrix_to_tibble("HGN") %>% 
  gather(sample, cpm, -HGN) %>%
  mutate(type = "negative")

positive <- cpm[rownames(cpm) %in% p.cor, colonSamples] %>% 
  matrix_to_tibble("HGN") %>% 
  gather(sample, cpm, -HGN) %>%
  mutate(type = "positive")

plotHeatColon(bind_rows(positive, negative)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

#"effect size". Haven't really finished this part...
bind_rows(positive, negative) %>%
  inner_join(select(MGAB.Meta, sample, section)) %>% 
  group_by(section, HGN, type) %>%
  summarize(mean = mean(cpm)) %>%
  ungroup() %>%
  group_by(HGN) %>%
  mutate(abs.diff = abs(mean[section == 1] - mean[section == 4])) %>%
  ungroup() %>%
  arrange(desc(abs.diff), HGN)

plotLineColon(bind_rows(positive, negative)) +
  scale_colour_manual(values = col40())

selected <- corr[corr > 0.95 | corr < -0.95]
dataToSave <- tibble(gene = names(selected), correlation = selected)
#write.table(dataToSave, file = "", sep = "\t", quote = FALSE, row.names = FALSE)

#highlights
low_genes <- c(
  "Lig4", "Kit", "Ccpg1", "Nod2", "Casp6", "Foxd2", "Fbxo10", "Ccl5", 
  "Il10rb", "Cxxc5", "Cdx2", "Pdgfrl", "Stat6", "Cdkn1c", "Mapk4", "Braf", 
  "Etv5", "Casp3", "Foxa1", "Prom1", "Ccl9", "Cdc34", "Sox6", "Il2rb", "Cdk2",
  "Gata4", "Cdk17", "Cd7", "Runx3", "Cd226", "Cd276", "Plet1os"
)  
high_genes <- c(
  "Map2k7", "Skil", "Lgr4", "Klf2", "Cdk8", "Il22ra1", "Jund", "Cd81", "Smad4", 
  "Cd82", "Cdh1", "Casp9", "Cdx1", "Foxp1", "Hoxb13", "Max", "Il17ra", "Myb", 
  "Sox13", "Amotl1", "Jak3", "Mapk15", "Brd1", "Mapk15", "Cbx3", "Bmp2", 
  "Foxp2", "Map3k13", "Il1rn"
)

negative <- cpm[rownames(cpm) %in% low_genes, ] %>% 
  matrix_to_tibble("HGN") %>% 
  gather(sample, cpm, -HGN) %>%
  mutate(type = "negative")

positive <- cpm[rownames(cpm) %in% high_genes, ] %>% 
  matrix_to_tibble("HGN") %>% 
  gather(sample, cpm, -HGN) %>%
  mutate(type = "positive")

plotHeatColon(bind_rows(positive, negative))

#identify genes not linear
nlg <- tibble(
  direction = 1:2,
  condition1 = list(
    c("NJE003C01", "NJE004C01"), c("NJE003C04", "NJE004C04")
  ),
  condition2 = list(
    c("NJE003C02", "NJE004C02", "NJE003C03", "NJE004C03", "NJE003C04", "NJE004C04"),
    c("NJE003C01", "NJE004C01", "NJE003C02", "NJE004C02", "NJE003C03", "NJE004C03")
  )
) %>%
  mutate(
    sample1 = map(condition1, function(c) {
      cpm[, c]
    }),
    sample2 = map(condition2, function(c) {
      cpm[, c]
    })
  ) %>%
  mutate(genes = map2(sample1, sample2, function(s1, s2) {
    lgl <- map_lgl(1:nrow(s1), function(i) {
      val <- s2[i, ]
      sd <- sd(val) * 1.25
      m <- mean(val)
      if(!all((val < m + sd) & (val > m - sd))) {
        FALSE
      } else {
        all(val > (mean(s1[i, ]) + (sd(s1[i, ]) * 2)))
      }
    })
    names(lgl) <- rownames(s1)
    names(lgl[lgl])
  }))

n.ends <- unlist(filter(nlg, direction == 1)$genes)
p.ends <- unlist(filter(nlg, direction == 2)$genes)

negative <- cpm[rownames(cpm) %in% n.ends, ] %>% 
  matrix_to_tibble("HGN") %>% 
  gather(sample, cpm, -HGN) %>%
  mutate(type = "negative")

positive <- cpm[rownames(cpm) %in% p.ends, ] %>% 
  matrix_to_tibble("HGN") %>% 
  gather(sample, cpm, -HGN) %>%
  mutate(type = "positive")

plotHeatColon(bind_rows(positive, negative)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())



#####
nlg <- tibble(
  condition1 = list(
    c("NJE003C01", "NJE004C01", "NJE003C02", "NJE004C02")
  ),
  condition2 = list(
    c("NJE003C03", "NJE004C03", "NJE003C04", "NJE004C04")
  )
) %>%
  mutate(
    sample1 = map(condition1, function(c) {
      cpm[, c]
    }),
    sample2 = map(condition2, function(c) {
      cpm[, c]
    })
  ) %>%
  mutate(genes = map2(sample1, sample2, function(s1, s2) {
    lgl <- map_lgl(1:nrow(s1), function(i) {
      val1 <- s2[i, ]
      sd1 <- sd(val1) * 1.25
      m1 <- mean(val1)
      if(!all((val < m + sd) & (val > m - sd))) {
        FALSE
      } else {
        all(val > (mean(s1[i, ]) + (sd(s1[i, ]) * 2)))
      }
    })
    names(lgl) <- rownames(s1)
    names(lgl[lgl])
  }))

