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
    git

# Install CRAN and Bioconductor packages
RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

# Clone and install remote R packages
RUN mkdir /home/Github

RUN git clone https://github.com/jasonserviss/sp.scRNAseq.git /home/Github/sp.scRNAseq

RUN git clone https://github.com/jasonserviss/sp.scRNAseqTesting.git /home/Github/sp.scRNAseqTesting

RUN Rscript -e "devtools::install('/home/Github/sp.scRNAseq')"
RUN Rscript -e "devtools::install('/home/Github/sp.scRNAseqTesting')"

WORKDIR /home/