build: update
	Rscript -e 'Rcpp::compileAttributes()' && \
		Rscript -e 'roxygen2::roxygenize()' && \
		R CMD INSTALL .
update:
	rsync -avz ../../research/world-epi/epiworld.hpp inst/include/epiworld.hpp 

.PHONY: build update
