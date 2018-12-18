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
    emacs \
    git

# Install CRAN and Bioconductor packages
RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl', 'BiocManager'), repos = 'https://cran.rstudio.com')"

##CRAN Package Imports
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('googledrive/0.1.1')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('openxlsx/4.0.17')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('Rtsne/0.13')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('pso/1.0.3')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('matrixStats/0.53.1')"
RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('ggthemes/3.5.0')"
#RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R'); install_cran('igraph/1.2.1')"
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

#Bioconductor package imports
RUN Rscript -e "BiocManager::install('S4Vectors')"

# Clone and install EngeMetadata
RUN mkdir /home/Github

RUN git clone https://github.com/EngeLab/EngeMetadata.git /home/Github/EngeMetadata
RUN Rscript -e "devtools::install('/home/Github/EngeMetadata', dependencies = FALSE)"

# Clone and install sp.scRNAseqData
RUN git clone https://github.com/jasonserviss/sp.scRNAseqData.git /home/Github/sp.scRNAseqData
RUN Rscript -e "source('/home/Github/sp.scRNAseqData/inst/rawData/processRaw.R')"
RUN Rscript -e "devtools::install('/home/Github/sp.scRNAseqData', dependencies = FALSE)"

WORKDIR /home/
