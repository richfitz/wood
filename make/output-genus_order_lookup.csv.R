#!/usr/bin/env Rscript
source("R/load.R")

## Build the genus -> family -> order lookup table.

## We do this in two steps.
##   1. Get the mapping of genus -> family from The Plant List
##   2. Get the mapping of family -> order from Zanne et al.

tpl <- load.theplantlist()
tpl <- unique(tpl[c("genus", "family", "major.clade")])
rownames(tpl) <- NULL

## Genus -> Family -> Order lookup table from Zanne et al., augmented
## by Dave Tank:
lookup1 <- read.csv("data/zae/genus_order_lookup.csv",
                    stringsAsFactors=FALSE)
lookup2 <- read.csv("data/genus_order_lookup_extra.csv",
                    stringsAsFactors=FALSE)
lookup <- rbind(lookup1, lookup2)
rownames(lookup) <- NULL

## Assign order to the Plant List families using this lookup:
tpl$order <- lookup$order[match(tpl$family, lookup$family)]

## TODO: We miss a few:
missing <- unique(tpl$family[is.na(tpl$order)])

## For now just nail them down so that we can keep going.
tpl$order[is.na(tpl$order)] <- ""

## Continuing:

## There are a handful of essentially unplaced families.  For now,
## these get their own pseudo-family
i <- tpl$order == ""
tpl$order[i] <- paste0("UnknownOrder-", tpl$family[i])

tpl <- tpl[c("genus", "family", "order", "major.clade")]

write.csv(tpl, "output/genus_order_lookup.csv", row.names=FALSE)
