#!/usr/bin/env Rscript

source("R/load.R")
source("R/build.R") # parse.count / summarise.count

## Aggressive warnings: turn into errors
options(warn=2)

## Where we'll store data:
dir.create("output", FALSE)

## Start by getting the woodiness information from the database
dat <- read.csv("data/zae/GlobalWoodinessDatabase.csv",
                stringsAsFactors=FALSE)
dat$woodiness[dat$woodiness == "variable"] <- "V"
names(dat)[names(dat) == "gs"] <- "species"
dat$species <- sub(" ", "_", dat$species)

## Check that we do recover the ZAE classes:
if (!identical(dat$woodiness,
               summarise.count(parse.count(dat$woodiness.count))))
  stop("Database classification failure")

## Map some species names to synonyms:
syn <- read.csv("data/synonyms.csv", stringsAsFactors=FALSE)

idx <- match(dat$species, syn$synonym)
message(sprintf("Resolving synonomy for %d species",
                sum(!is.na(idx))))
i <- !is.na(idx)
dat$species[i] <- syn$species[idx[i]]

## Read in The Plant List
tpl <- read.csv("data/theplantlist/names_accepted.csv",
                stringsAsFactors=FALSE)
keep <- dat$species %in% tpl$gs
message(sprintf("Dropping %d species not in Plant List",
                sum(!keep)))
dat <- dat[keep,]

## And look to see which species have now got duplicated records due
## to synonomy resolution:
dups <- unique(sort(dat$species[duplicated(dat$species)]))
message(sprintf("After synonym correction, %d duplicated entries",
                length(dups)))

## Merge the counts across the different instances of these names:
merge.counts <- function(x)
  paste(colSums(parse.count(x)), collapse=";")
dups.i <- which(dat$species %in% dups)
dups.count <- tapply(dat$woodiness.count[dups.i],
                     dat$species[dups.i], merge.counts)
dups.woodiness <- summarise.count(parse.count(dups.count))

dups.merged <- data.frame(species=names(dups.woodiness),
                          woodiness=unname(dups.woodiness),
                          woodiness.count=unname(dups.count),
                          stringsAsFactors=FALSE, row.names=NULL)

## Drop the duplicated records from the original vector, and add in
## the resolved entries here:
dat <- rbind(dat[-dups.i,], dups.merged)

## Tidy up, and we're done
dat <- dat[order(dat$species),
           c("species", "woodiness", "woodiness.count")]
names(dat)[1] <- "gs"
rownames(dat) <- NULL

saveRDS(dat, "output/woodiness.rds")
