.PHONY: format check all install-deps

install-deps:
	Rscript -e "if (!requireNamespace('styler', quietly = TRUE)) install.packages('styler')"

format: 
	Rscript -e "styler::style_dir('R')"
	Rscript -e "styler::style_dir('tests')"


all: install-deps format
