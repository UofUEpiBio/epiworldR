FROM rocker/r-base:latest

RUN \
  echo 'options(repos=c(CRAN="https://cloud.r-project.org"))' >> ~/.Rprofile && \
  Rscript --vanilla -e 'getOption("repos")'


RUN Rscript --vanilla -e "install.packages('cpp11', dependencies = TRUE, repos = 'https://cloud.r-project.org')"



RUN mkdir ~/.R && \
  echo "CXXFLAGS=-g -O0" > ~/.R/Makevars && \
  echo "CXX14FLAGS=-g -O0" >> ~/.R/Makevars && \
  echo "CXX17FLAGS=-g -O0" >> ~/.R/Makevars && \
  echo "CXX11FLAGS=-g -O0" >> ~/.R/Makevars && \
  echo "SAFE_FLAGS=-g -O0" >> ~/.R/Makevars

CMD ["bash"]