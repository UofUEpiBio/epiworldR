# Capture the current value of the version of the package in DESCRIPTION
VERSION := $(shell grep Version DESCRIPTION | sed -e 's/Version: //')

build: docs clean
	cd .. && rm -f epiworldR_*.tar.gz & R CMD build epiworldR

debug: clean
	docker run --rm -ti -w/mnt -v $(PWD):/mnt uofuepibio/epiworldr:debug make docker-debug 

docker-debug:
	EPI_CONFIG="-DEPI_DEBUG -Wall -pedantic -g" R CMD INSTALL \
		--no-docs --build .

install: build 
	which R
	cd ../ && R CMD INSTALL epiworldR_*.tar.gz

README.md: README.Rmd
	Rscript --vanilla -e 'rmarkdown::render("README.Rmd")'

# update:
# 	wget https://raw.githubusercontent.com/UofUEpiBio/epiworld/master/epiworld.hpp && \
# 		mv epiworld.hpp inst/include/epiworld.hpp 
local-update:
	rsync -avz ../epiworld/include/epiworld inst/include/.

check: build
	R CMD check epiworldR_*.tar.gz

clean: 
	Rscript --vanilla -e 'devtools::clean_dll()'

docs:
	Rscript --vanilla -e 'devtools::document()'

.PHONY: build update check clean docs docker-debug

checkv: build
	R CMD check --as-cran --use-valgrind epiworldR*.tar.gz
	
