#!/usr/bin/env Rscript
source("R/load.R")
source("R/build.R")
saveRDS(build.woodiness.genus(extreme="herbaceous"),
        filename.woodiness.genus("herbaceous"))
