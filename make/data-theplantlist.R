#!/usr/bin/env Rscript
source("R/theplantlist.R")

path <- "data/theplantlist"
dir.create(file.path(path, "acceptedNames1.1"), FALSE)

plant.list.get.group("gymnosperm")
plant.list.get.group("pteridophyte")
plant.list.get.group("angiosperm")

files <- dir(file.path(path, "acceptedNames1.1"), full.names=TRUE)
dat <- do.call(rbind, lapply(files, load.names))
dat$gs <- paste(dat$genus, dat$species, sep="_")

dat.unique <- dat[!duplicated(dat),]
dat.accepted <- dat.unique[dat.unique$status == "Accepted",]

write.csv.pl <- function(obj, filename)
  write.csv(obj, file.path(path, filename), row.names=FALSE)

write.csv.pl(dat.accepted, "names_accepted.csv")
