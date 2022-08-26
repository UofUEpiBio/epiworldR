build: update
	Rscript -e 'Rcpp::compileAttributes()' && \
		Rscript -e 'roxygen2::roxygenize()' && \
		R CMD INSTALL .

README.md: README.Rmd
	Rscript --vanilla -e 'rmarkdown::render("README.Rmd")'
update:
	rsync -avz ../../research/world-epi/epiworld.hpp inst/include/epiworld.hpp 

.PHONY: build update
