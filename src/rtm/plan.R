#! /usr/bin/env Rscript
require(here)
require(drake)

rmdpath <- here("presentations/rtm.Rmd")
outpath <- here("output/rtm.pdf")

print(rmdpath)
plan <- drake_plan(
  rmarkdown::render(
    knitr_in("/home/jzelner/Dropbox/repos/projects/hierarchical-nov/presentations/rtm.Rmd"),
    output_file = file_out("/home/jzelner/Dropbox/repos/projects/hierarchical-nov/output/rtm.pdf")
  )
)

make(plan)