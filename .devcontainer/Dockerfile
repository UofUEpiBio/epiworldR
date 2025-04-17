FROM rocker/tidyverse:4.4.0

# RUN \
#   echo 'options(repos=c(CRAN="https://packagemanager.posit.co/cran/__linux__/bookworm/latest"))' >> ~/.Rprofile && \
#   Rscript --vanilla -e 'getOption("repos")'

# Adding R packages
RUN install2.r cpp11 roxygen2 tinytest data.table netplot \
  devtools decor

RUN install2.r languageserver httpgd

RUN apt-get update && apt-get install --no-install-recommends -y valgrind gdb

CMD ["bash"]
