#!/usr/bin/env Rscript
source("wood-functions.R")
dir.create("output", FALSE)
dat <- load.woodiness.data(TRUE)
