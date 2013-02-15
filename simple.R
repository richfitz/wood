## % How much of the world is woody, first analysis.
library(multicore)
library(diversitree)
library(RColorBrewer)
path.forest <- readLines("~/.forest_path")

read.forest.csv <- function(filename)
  read.csv(file.path(path.forest, filename), stringsAsFactors=FALSE)

## Start by getting the woodiness information from the database
dat <- read.forest.csv("export/speciesTraitData.csv")

## Score the 633 species with no known information as NA
dat$woodiness[!(dat$woodiness %in% c("H", "W")) &
              !is.na(dat$woodiness)] <- NA

## Only the columns we care about:
dat <- data.frame(species=sub(" ", "_", dat$gs),
                  woodiness=dat$woodiness,
                  stringsAsFactors=FALSE)
## Filtered by whether or not they have woodiness information
dat <- dat[!is.na(dat$woodiness),]

## Next, normalise the species names.
spp <-
  read.forest.csv("taxonomic/spermatophyta_synonyms_PLANTLIST.csv")
spp$genus.synonym <- sub("_.+", "", spp$synonym)

## Locate the species names within the Plant List lookup table
idx <- match(dat$species, spp$synonym)
## 23% of 
mean(is.na(idx))

## Ouch -- here are all the names that are not found in the Plant
## list.
## 
## Will: This seems too high given the synonomy resolution that you
## did when exporting data.  Can you shed any light on this?
dat$species[is.na(idx)]

## So, for example we don't have
##   Ziziphus_itacaiunensis
## in the plant list, but we have species information for it.
##
## Using the Plant List online, this name comes back as "unresolved".
## It is not found in here at all (e.g., by misspelling).

## TODO: Decide what to do with this huge number of misses.  In the
## mean time I am deleting them all.
dat <- dat[!is.na(idx),]
idx <- match(dat$species, spp$synonym) # rematch

## About 90% of names are valid:
mean(spp$valid[idx])

## Assign the correct/current species name and adding the genus
dat$species <- spp$species[idx]
dat$genus   <- spp$genus[idx]

## And look to see which species have now got duplicated records due
## to synonomy resolution:
dups <- unique(sort(dat$species[duplicated(dat$species)]))

## 1436 species with more than one entry now.
length(dups)

## Find the indices of all entries corresponding to duplicated
## species:
dups.i <- which(dat$species %in% dups)

## Most duplicated species are of a single type (good) but a few
## aren't:
f <- function(x)
  if ( length(unique(x)) == 1 ) x[1] else NA
dups.fixed <- tapply(dat$woodiness[dups.i], dat$species[dups.i], f)
dups.fixed <- data.frame(species=names(dups.fixed),
                         woodiness=dups.fixed,
                         stringsAsFactors=FALSE, row.names=NULL)
dups.fixed$genus <- spp$genus[match(dups.fixed$species, spp$species)]

## These species had conflicting records and will be dropped:
dups.fixed[is.na(dups.fixed$woodiness),]
dups.fixed <- dups.fixed[!is.na(dups.fixed$woodiness),]

## Drop the duplicated records from the original vector, and add in
## the resolved entries here:
dat <- rbind(dat[-dups.i,], dups.fixed)

## So now 'dat' has species names sanitised to the same list that the
## plant list uses, and includes genus information.

## Next, compare the genus-> order lookup with the list of species
## that we have in the plant list.
lookup <- read.forest.csv("taxonomic/genus_order_lookup.csv")
lookup <- lookup[c("genus", "family", "order")]

## There are a handful of essentially unplaced families.  For now,
## these get their own pseudo-family
i <- lookup$order == ""
lookup$order[i] <- paste0("UnknownOrder-", lookup$family[i])

## First, check that the Plant List and the lookup table agree on the
## number of genera:

## TODO:
## There are three genera in the species list that are not found in
## the lookup table.
setdiff(spp$genus[spp$valid], lookup$genus)
## However, we don't have any data for these.

## More concerning -- there are a large number of genera in the lookup
## that we don't know about.
setdiff(lookup$genus, spp$genus[spp$valid])

## Here are the 1700 genera that are in the lookup table that are not
## in the valid species list.
tmp <- lookup[lookup$genus %in% setdiff(lookup$genus, spp$genus[spp$valid]),]

## If we allow for invalid species (i.e., the genus-family-order
## relationships are known for genera that are only known as synonyms,
## which is reasonable) then there are 934 genera that aren't in the
## species list but that we have order information for.  These include
## Lamiales, etc so the problem is not non-flowering plants again.
tmp <- lookup[lookup$genus %in% setdiff(lookup$genus, spp$genus.synonym),]

## Pressing on anyway.  All of the *data* genera are present in all
## cases.
all(dat$genus %in% spp$genus)
all(dat$genus %in% lookup$genus)

## Collapse these to get counts for all the genera that we know about.
dat.g <- table(factor(dat$genus, sort(unique(spp$genus[spp$valid]))),
               factor(dat$woodiness))
dat.g <- data.frame(genus=rownames(dat.g),
                    W=dat.g[,"W"],
                    H=dat.g[,"H"],
                    K=rowSums(dat.g),
                    stringsAsFactors=FALSE)
## Include the counts of known species:
spp.known <- table(spp$genus[spp$valid])
dat.g$N <- spp.known[dat.g$genus]

## There are no cases for which we know about more species than exist,
## which is a good check.
dat.g[which(dat.g$K > dat.g$N),]

## 1. Simple minded test.

## There are three different types of clades:
## u: Unknown -- no species known
## w: Woody -- all known species woody
## h: herbaceous -- all known species herbaceous
## p: polymorphic -- both woody and herbaceous species known
dat.g.u <- dat.g[dat.g$W == 0 & dat.g$H == 0,]
dat.g.w <- dat.g[dat.g$W > 0  & dat.g$H == 0,]
dat.g.h <- dat.g[dat.g$W == 0 & dat.g$H >  0,]
dat.g.p <- dat.g[dat.g$W > 0  & dat.g$H >  0,]

c(nrow(dat.g.u), sum(dat.g.u$N)) # unknown     -- 6950 genera,  35581 spp
c(nrow(dat.g.w), sum(dat.g.w$N)) # woody       -- 3890 genera,  86735 spp
c(nrow(dat.g.h), sum(dat.g.h$N)) # herbaceous  -- 3013 genera, 101904 spp
c(nrow(dat.g.p), sum(dat.g.p$N)) # polymorphic --  261 genera,  49964 spp

## So if we fill in the monomorphic genera we get
## woody, herbacious, still unknown
sum(dat.g.w$N) + sum(dat.g.p$W) # 88850
sum(dat.g.h$N) + sum(dat.g.p$H) # 103774
sum(dat.g.u$N) + sum(dat.g.p$N - dat.g.p$K) # 81560 species

## So we are in a range of 32% to 62% woody based on this approach.
c(88850, 88850+81560) / sum(dat.g$N)

## Interestingly, if we split the difference here (randomly assign
## species as woody/herbaceous for the 81560 unknown species) we get
(88850+81560/2) / sum(dat.g$N)
## 47.3%, which is rather similar to the 48.75 that we were getting
## before.  Could be coincidence.

## Looking at the polymorphic cases, these are pretty evenly spread
## over the distribution, so these are going to come out 50:50
## woody/herbacious which will drag the estimate towards the middle of
## its range.
dat.g.p$p <- dat.g.p$W / dat.g.p$K
hist(dat.g.p$p)
