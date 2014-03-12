#!/usr/bin/env Rscript
source("R/load.R")
source("R/build.R")
saveRDS(build.woodiness.genus(extreme=FALSE),
        filename.woodiness.genus(FALSE))
