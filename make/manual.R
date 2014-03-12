#!/usr/bin/env Rscript

## Emulate Rscript by not modifying the global environment.
Rscript <- function(filename) {
  source(filename, local=TRUE)
}

## 'data-raw'
if (file.exists(".downloaded_data.tar.gz")) {
  # Use the internal untar rather than system tar, because playforms
  # that don't have make are unlikely to have a good tar, and because
  # the defaults for external tar don't work on OSX.
  untar(".downloaded_data.tar.gz", tar="internal")
}

Rscript("make/data-zae.R")
Rscript("make/data-theplantlist.R")

## 'data-processed'
Rscript("make/output-genus_order_lookup.csv.R")
Rscript("make/output-woodiness.rds.R")
Rscript("make/output-dat.g.rds.R")
Rscript("make/output-dat.g.w.rds.R")
Rscript("make/output-dat.g.h.rds.R")
Rscript("make/output-phy.o.rds.R")

## wood.R -> wood.Rmd
library(sowsear)
sowsear('wood.R', 'Rmd')

## wood.Rmd -> wood.md
library(knitr)
knit('wood.Rmd', envir=new.env())

## wood.md -> wood.html
library(markdown)
markdownToHTML('wood.md', 'wood.html')
