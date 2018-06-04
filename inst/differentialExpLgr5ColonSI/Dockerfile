FROM rocker/tidyverse:3.5.0

MAINTAINER Jason Serviss <jason.serviss@ki.se>

# System dependencies for required R packages
RUN  rm -f /var/lib/dpkg/available \
  && rm -rf  /var/cache/apt/* \
  && apt-get update -qq \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libudunits2-dev \
    emacs \
    git 

# Install CRAN and Bioconductor packages
RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('mclust/5.4', 'Rtsne/0.13', 'pso/1.0.3', 'matrixStats/0.53.1', 'ggthemes/3.5.0', 'igraph/1.2.1', 'viridis/0.5.1', 'ggraph/1.0.1', 'tidygraph/1.1.0', 'testthat/2.0.0', 'print2/0.1', 'covr/3.1.0', 'openxlsx/4.0.17', 'googledrive/0.1.1', 'rlang/0.2.0', 'dplyr/0.7.5', 'tidyr/0.8.1', 'e1071/1.6-8', 'doMC/1.3.5'))"

RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('S4Vectors', 'BiocStyle', 'BioStrings'))"

# Clone and install remote R packages
RUN mkdir /home/Github
RUN git clone https://github.com/jasonserviss/sp.scRNAseqData.git /home/Github/sp.scRNAseqData
WORKDIR /home/Github/sp.scRNAseqData
RUN Rscript -e "devtools::install('/home/Github/sp.scRNAseqData')"
RUN Rscript -e "source('/home/Github/sp.scRNAseqData/inst/rawData/processRaw.R')"
RUN Rscript -e "devtools::install('/home/Github/sp.scRNAseqData')"

RUN git clone https://github.com/jasonserviss/seqTools.git /home/Github/seqTools

WORKDIR /home/