# Capture the current value of the version of the package in DESCRIPTION
VERSION := $(shell grep Version DESCRIPTION | sed -e 's/Version: //')


build:
	cd .. && R CMD build epiworldR

debug: clean
	docker run --rm -ti -w/mnt -v $(PWD):/mnt uofuepibio/epiworldr:debug make docker-debug

docker-debug:
	EPI_CONFIG="-DEPI_DEBUG -Wall -pedantic -g" R CMD INSTALL \
		--no-docs --build .

install-dev: clean
	sed -i -E 's/@useDynLib\s+[a-zA-Z]+/@useDynLib epiworldRdev/g' R/epiworldR-package.R
	sed -i -E 's/useDynLib\(+[a-zA-Z]+/useDynLib(epiworldRdev/g' NAMESPACE
	sed -i -E 's/^Package:.+/Package: epiworldRdev/g' DESCRIPTION
	sed -i -E 's/^library\(epiworldR\)/library(epiworldRdev)/g' README.*
	Rscript --vanilla -e 'roxygen2::roxygenize()'
	EPI_DEV=yes R CMD INSTALL .& $(MAKE) clean

install:
	cd .. && \
		R CMD INSTALL epiworldR_$(VERSION).tar.gz


README.md: README.Rmd
	Rscript --vanilla -e 'rmarkdown::render("README.Rmd")'

# update:
# 	wget https://raw.githubusercontent.com/UofUEpiBio/epiworld/master/epiworld.hpp && \
# 		mv epiworld.hpp inst/include/epiworld.hpp
local-update:
	rsync -avz ../epiworld/include/epiworld inst/include/.

check: build
	cd .. && R CMD check epiworldR_*.tar.gz

clean:
	rm -f src/*.dll src/*.so src/*.o
	sed -i -E 's/@useDynLib\s+[a-zA-Z]+/@useDynLib epiworldR/g' R/epiworldR-package.R
	sed -i -E 's/useDynLib\(+[a-zA-Z]+/useDynLib(epiworldR/g' NAMESPACE
	sed -i -E 's/^Package:.+/Package: epiworldR/g' DESCRIPTION
	# sed -i -E 's/^\\(name|alias|title)\{[a-zA-Z]+/\\\1{epiworldR/g' man/epiworldR-package.Rd
	sed -i -E 's/^library\(epiworldRdev\)/library(epiworldR)/g' README.*

docs:
	Rscript --vanilla -e 'roxygen2::roxigenize()'

.PHONY: build update check clean docs docker-debug

checkv: build
	R CMD check --as-cran --use-valgrind epiworldR*.tar.gz
