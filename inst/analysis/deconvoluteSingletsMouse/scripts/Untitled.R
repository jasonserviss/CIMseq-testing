library(tidyverse)
library(CIMseq)

load('~/Github/CIMseq.testing/inst/analysis/MGA.analysis_enge/data/CIMseqData.rda')
load('~/Github/CIMseq.testing/inst/analysis/deconvoluteSingletsMouse/data/sObj.rda')

cObjMul.2 <- CIMseqMultiplets(
  getData(cObjSng, "counts")[, rownames(getData(sObj, "fractions"))],
  getData(cObjSng, "counts.ercc")[, rownames(getData(sObj, "fractions"))],
  getData(cObjMul, "features")
)

renameClasses <- function(class) {
  case_when(
    class == "0" ~ "C.Goblet.d",
    class == "1" ~ "C.Lgr5+.p.1",
    class == "2" ~ "C.Lgr5+.d",
    class == "3" ~ "SI.Lgr5+.KI67h",
    class == "4" ~ "C.TA.d",
    class == "5" ~ "C.Lgr5+.p.2",
    class == "6" ~ "SI.Lgr5+.KI67l",
    class == "7" ~ "C.Lgr5+.KI67+.p",
    class == "8" ~ "C.Colonocytes.p",
    class == "9" ~ "SI.Entero.early",
    class == "10" ~ "SI.Lgr5-.KI67l",
    class == "11" ~ "SI.Lgr5+.KI67-",
    class == "12" ~ "SI.Goblet",
    class == "13" ~ "C.Goblet.p",
    class == "14" ~ "C.Colonocyte.d",
    class == "15" ~ "SI.Entero",
    class == "16" ~ "Enteroendo",
    class == "17" ~ "SI.Lgr5-.KI67+",
    class == "18" ~ "Tufft",
    class == "19" ~ "C.Lgr5+.KI67+.d",
    class == "20" ~ "SI.Paneth",
    class == "21" ~ "SI.Goblet.KI67+",
    class == "22" ~ "C.Goblet.KI67+",
    class == "23" ~ "Blood",
    TRUE ~ "error"
  )
}

getData(cObjSng, "classification") <- renameClasses(getData(cObjSng, "classification"))
fractions <- getData(sObj, "fractions")
colnames(fractions) <- renameClasses(colnames(fractions))
sObj@fractions <- fractions

data <- getCellsForMultiplet(sObj, cObjSng, cObjMul.2, theoretical.max = 4, drop = FALSE) %>%
  inner_join(tibble(sample = colnames(getData(cObjSng, "counts")), class = getData(cObjSng, "classification")))

data %>%
  mutate(correct = cells == class) %>%
  group_by(class) %>%
  summarize(precision = sum(correct) / n()) %>%
  arrange(precision) %>%
  print(n = nrow(.))

data %>% 
  count(class, cells) %>%
  group_by(class) %>%
  mutate(n = n / sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_tile(aes(class, cells, fill = n)) +
  scale_fill_viridis_c() +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))