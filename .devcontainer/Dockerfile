FROM rocker/r-devel:latest

# RUN \
#   echo 'options(repos=c(CRAN="https://cloud.r-project.org"))' >> ~/.Rprofile && \
#   Rscript --vanilla -e 'getOption("repos")'

# Adding Git
RUN apt-get update && apt-get install -y --no-install-recommends git

# Adding R packages
RUN \
  wget https://github.com/jgm/pandoc/releases/download/3.2.1/pandoc-3.2.1-1-amd64.deb && \
  dpkg -i pandoc-3.2.1-1-amd64.deb

RUN install2.r cpp11 rmarkdown roxygen2 tinytest data.table netplot

RUN apt-get install -y --no-install-recommends && \
  install2.r languageserver

CMD ["bash"]
