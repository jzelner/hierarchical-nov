###############################################################################
## Data figures
output/figures/all_camp_cases%pdf output/figures/cases_by_camp%pdf : src/data_figure%R data/jamboree_data.csv
	@mkdir -p $(@D)
	./$<

DATAFIGURES := output/figures/all_camp_cases.pdf output/figures/cases_by_camp.pdf
.PHONY: datafigures
datafigures : $(DATAFIGURES)

###############################################################################
## Setup for running model
output/nov_model_input.Rds : src/prep_model_input.R data/jamboree_data.csv data/jamboree_pop.csv
	@mkdir -p $(@D)
	./$<

output/nov_model.Rds : src/runmodel.R src/contact_model.stan output/nov_model_input.Rds
	@mkdir -p $(@D)
	./$<


###############################################################################
## Setup for running multiple models for sensitivity
output/random_model_runs.Rds : src/random_model_runs.R data/jamboree_data.csv data/jamboree_pop.csv src/contact_model.stan
	@mkdir -p $(@D)
	./$<

output/combined_values.Rds : src/combine_outputs.R output/random_model_runs.Rds
	@mkdir -p $(@D)
	./$<

###############################################################################
## Post-processing
output/figures/daily_avg_r.pdf : src/avg_r_plot.R output/combined_values.Rds
	@mkdir -p $(@D)
	./$<

output/scalar_pars.csv : src/parameter_ci.R output/combined_values.Rds
	@mkdir -p $(@D)
	./$<

output/infection_source.Rds : src/inf_source.R output/nov_model.Rds output/nov_model_input.Rds
	@mkdir -p $(@D)
	./$<

output/figures/p_by_camp.pdf : src/p_by_camp.R output/infection_source.Rds
	@mkdir -p $(@D)
	./$<

output/figures/marginal_inf_dist.pdf : src/marginal_r.R output/nov_model.Rds
	@mkdir -p $(@D)
	./$<
	
	
FIGURES := output/figures/daily_avg_r.pdf output/figures/p_by_camp.pdf output/figures/daily_camp_r.pdf
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
output/%.md : presentations/%.Rmd output/scalar_pars.csv output/infection_source.Rds 
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
output/methods.pdf : output/methods.tex $(FIGURES)
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
