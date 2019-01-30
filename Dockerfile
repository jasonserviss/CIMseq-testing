FROM rocker/tidyverse:3.5.1

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
    libhdf5-dev \
    libcairo2-dev \
    emacs \
    git \
    python-dev \
    python-pip

# Install app dependencies
RUN pip install --upgrade pip
RUN pip install "umap-learn==0.2.3"

# Install CRAN and Bioconductor packages
RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl', 'BiocManager'), repos = 'https://cran.rstudio.com')"

##CRAN Package Imports
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('googledrive/0.1.1')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('openxlsx/4.0.17')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('Rtsne/0.13')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('pso/1.0.3')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('matrixStats/0.53.1')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('ggthemes/3.5.0')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('viridis/0.5.1')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('ggraph/1.0.1')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('tidygraph/1.1.0')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('Rcpp/0.12.19')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('future.apply/0.2.0')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('RANN/2.6')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('gmodels/2.18.1')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('circlize/0.4.4')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('gridBase/0.4-7')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('Seurat/2.3.4')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('RcppArmadillo/0.9.200.5.0')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('rmarkdown/1.11')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('printr/0.1')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('knitr/1.20')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('remotes/2.0.2')"
RUN Rscript -e "remotes::install_github('thomasp85/patchwork@fd7958bae3e7a1e30237c751952e412a0a1d1242')"

#Bioconductor package imports
RUN Rscript -e "BiocManager::install('S4Vectors')"

# Clone and install EngeMetadata
RUN mkdir ~/Github
RUN git clone https://github.com/EngeLab/EngeMetadata.git ~/Github/EngeMetadata
RUN Rscript -e "devtools::install('~/Github/EngeMetadata', dependencies = FALSE)"

# Clone and install CIMseq-data
RUN git clone https://github.com/jasonserviss/CIMseq.data.git ~/Github/CIMseq.data
RUN Rscript -e "devtools::install('~/Github/CIMseq.data', dependencies = FALSE)"
RUN cd ~/Github/CIMseq.data && Rscript -e "source('~/Github/CIMseq.data/inst/rawData/processRaw.R'); processRaw(ignore = ignore, upload = FALSE, save = TRUE)"
RUN Rscript -e "devtools::install('~/Github/CIMseq.data', dependencies = FALSE)"

# Clone and install CIMseq
RUN git clone https://github.com/jasonserviss/CIMseq.git --branch devel ~/Github/CIMseq
RUN Rscript -e "devtools::install('~/Github/CIMseq', dependencies = FALSE)"

# Clone and install CIMseq-testing
RUN touch /tmp33.txt
RUN git clone https://github.com/jasonserviss/CIMseq.testing.git ~/Github/CIMseq.testing
RUN Rscript -e "devtools::install('~/Github/CIMseq.testing', dependencies = FALSE)"
