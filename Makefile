# Capture the current value of the version of the package in DESCRIPTION
VERSION := $(shell grep Version DESCRIPTION | sed -e 's/Version: //')

help:
	@echo "Makefile commands:"
	@echo "  build               Build the package tarball"
	@echo "  check               Check the built package"
	@echo "  checkv              Check the built package with valgrind"
	@echo "  install             Install the built package"
	@echo "  install-dev         Install the development version of the package"
	@echo "  debug               Build and run a debug Docker container"
	@echo "  docs                Generate documentation with roxygen2"
	@echo "  local-update        Update epiworld.hpp from local epiworld include folder"
	@echo "  local-update-diagrams  Update diagrams from local epiworld docs_src folder"
	@echo "  website             Build the pkgdown website"
	@echo "  clean               Clean compiled files and reset DESCRIPTION and NAMESPACE"

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


README.md: README.qmd
	quarto render README.qmd

# update:
# 	wget https://raw.githubusercontent.com/UofUEpiBio/epiworld/master/epiworld.hpp && \
# 		mv epiworld.hpp inst/include/epiworld.hpp
local-update:
	rsync -avz --delete ../epiworld/include/epiworld inst/include/. && \
	rm inst/include/epiworld/models/*.mmd

local-update-diagrams:
	rsync -avz --delete ../epiworld/docs_src/assets/img/* man/figures/

update:
	git clone --depth=1 https://github.com/UofUEpiBio/epiworld tmp_epiworld && \
	rsync -avz --delete tmp_epiworld/include/epiworld inst/include/. && \
	rm inst/include/epiworld/models/*.mmd && \
	rsync -avz --delete tmp_epiworld/docs_src/assets/img/* man/figures && \
	rm -rf tmp_epiworld


check:
	Rscript --vanilla -e 'devtools::check()'

clean:
	rm -f src/*.dll src/*.so src/*.o
	sed -i -E 's/@useDynLib\s+[a-zA-Z]+/@useDynLib epiworldR/g' R/epiworldR-package.R
	sed -i -E 's/useDynLib\(+[a-zA-Z]+/useDynLib(epiworldR/g' NAMESPACE
	sed -i -E 's/^Package:.+/Package: epiworldR/g' DESCRIPTION
	# sed -i -E 's/^\\(name|alias|title)\{[a-zA-Z]+/\\\1{epiworldR/g' man/epiworldR-package.Rd
	sed -i -E 's/^library\(epiworldRdev\)/library(epiworldR)/g' README.*

docs:
	Rscript -e 'devtools::document()'

checkv: build
	R CMD check --as-cran --use-valgrind epiworldR*.tar.gz

# Builds and installs without vignettes
dev: clean
	R CMD build --no-build-vignettes .
	R CMD INSTALL epiworldR_$(VERSION).tar.gz

website:
	Rscript -e 'pkgdown::build_site()'


test:
	Rscript --vanilla -e 'devtools::load_all(); tinytest::test_all()'

.PHONY: build update check clean docs docker-debug dev website test
