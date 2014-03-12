#!/usr/bin/env Rscript
library(diversitree, quietly=TRUE)
suppressMessages(library(dplyr))
source("R/load.R")

drop.tip <- diversitree:::drop.tip.fixed

lookup <- load.genus.order.lookup()
phy <- read.tree("data/zae/Vascular_Plants_rooted.dated.tre")

## Try this in the simplest possible way.  We'll take one
## representative per order.
phy.order <- lookup$order[match(sub("_.+$", "", phy$tip.label),
                                lookup$genus)]
phy.order[grepl("^UnknownOrder-", phy.order)] <- NA

keep <- data.frame(tip=phy$tip.label, order=phy.order,
                   stringsAsFactors=FALSE) %.%
  filter(!is.na(order)) %.%
  group_by(order) %.%
  summarise(tip=tip[[1]])

## Update this for a few cases:
keep.change <- c("Cyatheales"="Dicksonia_antarctica",
                 "Ophioglossales"="Ophioglossum_lusitanicum")
keep$tip[match(names(keep.change), keep$order)] <- unname(keep.change)

phy.o <- drop.tip(phy, setdiff(phy$tip.label, keep$tip))
phy.o$tip.label <- keep$order[match(phy.o$tip.label, keep$tip)]
phy.o <- ladderize(phy.o)

## Get counts by order:
dat.g <- load.woodiness.genus()
n <- dat.g %.% group_by(order) %.% summarise(n=sum(N))

phy.o$n.taxa <- n$n[match(phy.o$tip.label, n$order)]
saveRDS(phy.o, "output/phy.o.rds")
