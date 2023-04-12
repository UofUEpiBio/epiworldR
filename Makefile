build: 
	Rscript -e 'Rcpp::compileAttributes()' && \
		Rscript -e 'roxygen2::roxygenize()' && \
		R CMD build .

debug: clean
	EPI_CONFIG=-DEPI_DEBUG R CMD INSTALL .

install: build
	R CMD INSTALL epiworldR_*tar.gz

README.md: README.Rmd
	Rscript --vanilla -e 'rmarkdown::render("README.Rmd")'

update:
	wget https://raw.githubusercontent.com/UofUEpi/epiworld/master/epiworld.hpp && \
		mv epiworld.hpp inst/include/epiworld.hpp 
local-update:
	rsync -avz ../epiworld/epiworld.hpp inst/include/epiworld.hpp

check: build
	R CMD check epiworldR_*.tar.gz

clean: 
	Rscript --vanilla -e 'devtools::clean_dll()'

.PHONY: build update check clean
