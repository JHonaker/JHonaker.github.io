#!/usr/bin/env Rscript

input <- commandArgs(trailingOnly = TRUE)
KnitPost <- function(input, base.url = "/") {
    require(knitr)
    opts_knit$set(base.url = base.url)
    fig.path <- paste0("../figs/", sub(".Rmd$", "", basename(input)), "/")
    opts_chunk$set(fig.path = fig.path)
    opts_chunk$set(fig.cap = "center")
    render_jekyll()
    outfile <- paste0("../_drafts/", Sys.Date(), '-', sub(".Rmd$", "", basename(input)), ".md")
    print(outfile)
    knit(input, output = outfile, envir = parent.frame())
}

KnitPost(input)
