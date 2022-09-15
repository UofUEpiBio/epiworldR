build: update
	Rscript -e 'Rcpp::compileAttributes()' && \
		Rscript -e 'roxygen2::roxygenize()' && \
		R CMD INSTALL .

README.md: README.Rmd
	Rscript --vanilla -e 'rmarkdown::render("README.Rmd")'
update:
	wget https://raw.githubusercontent.com/UofUEpi/epiworld/master/epiworld.hpp && \
		mv epiworld.hpp inst/include/epiworld.hpp 

.PHONY: build update
