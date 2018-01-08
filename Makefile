###############################################################################
## Setup for running model
output/nov_model_input.Rds : src/prep_model_input.R data/jamboree_data.csv data/jamboree_pop.csv
	@mkdir -p $(@D)
	./$<

output/nov_model.Rds : src/runmodel.R src/model.stan output/nov_model_input.Rds
	@mkdir -p $(@D)
	./$<

###############################################################################
## Post-processing
output/figures/daily_avg_r.pdf : src/avg_r_plot.R output/nov_model.Rds
	@mkdir -p $(@D)
	./$<

output/scalar_pars.csv : src/parameter_ci.R output/nov_model.Rds
	@mkdir -p $(@D)
	./$<
###############################################################################
## Below here is for generating output documents

UNAME_S := $(shell uname -s)
FONT :=
TEXENGINE := xelatex
BIBDIR := /home/jzelner/repos/bibtex-library/papers.bib
## Generate PDFs from TeX and Rmd files

FONT += --variable mainfont=Arial -V fontsize=11pt 
ifeq ($(UNAME_S), Darwin)
   TEXENGINE = xelatex
   BIBDIR = /Users/jzelner/Dropbox/bibtex/papers.bib
endif

CSL := presentations/config/american-journal-of-epidemiology.csl
############################################################
## Generic rule for translating Rmd to markdown
output/%.md : presentations/%.Rmd
	@echo ----Translating abstract text from RMD to Markdown----
	@mkdir -p $(@D)
	Rscript \
		-e "require(knitr)" \
	  -e "knitr::render_markdown()"\
		-e "knitr::knit('$<','$@')"

## Generic rule for translating markdown to tex
output/%.tex : presentations/%.md
	@echo ----Generating $(@F)----
	@mkdir -p $(@D)
	pandoc $< -V geometry:margin=0.75in --from=markdown -t latex -s -o $@ --bibliography=$(BIBDIR) --filter pandoc-citeproc --csl $(CSL) --latex-engine=xelatex $(FONT) --template=presentations/config/template.latex 

## Generic rule for translating markdown to tex
output/%.tex : output/%.md
	@echo ----Generating $(@F)----
	@mkdir -p $(@D)
	pandoc $< -V geometry:margin=0.75in --from=markdown -t latex -s -o $@ --bibliography=$(BIBDIR) --filter pandoc-citeproc --csl $(CSL) --latex-engine=xelatex $(FONT) --template=presentations/config/template.latex


## Generic rule for translating tex to pdf
output/%.pdf : output/%.tex $(FIGURES)
	@echo ----Generating $(@F)----
	latexmk -$(TEXENGINE) $< > texout.txt
	mv $(@F) $(@D)/.

## Generic rule for translating tex to pdf
output/%.pdf : presentations/%.tex
	@echo ----Generating $(@F)----
	latexmk -$(TEXENGINE) $< > texout.txt
	mv $(@F) $(@D)/.

.PHONY: pdf
pdf: output/methods.pdf

.PHONY: clean
clean:
	rm -rf output
