#!/usr/bin/env Rscript
source("R/load.R")
suppressMessages(library(dplyr))

## Build the genus -> family -> order lookup table.

## We do this in two steps.
##   1. Get the mapping of genus -> family from The Plant List
##   2. Get the mapping of family -> order from Zanne et al.
tpl.full <- load.theplantlist()
tpl <- unique(tpl.full[c("genus", "family", "major.clade")])
rownames(tpl) <- NULL

## Genus -> Family -> Order lookup table from Zanne et al., augmented
## by Dave Tank:
lookup1 <- read.csv("data/zae/genus_order_lookup.csv",
                    stringsAsFactors=FALSE)
lookup2 <- read.csv("data/genus_order_lookup_extra.csv",
                    stringsAsFactors=FALSE)
lookup <- rbind(lookup1, lookup2)
rownames(lookup) <- NULL

## There are a handful of problem cases in TPL: these are genera that
## are duplicated but have different families.  None of these are very
## large.
duplicates <- sort(tpl$genus[duplicated(tpl$genus)])

info <- tpl.full                %>%
  filter(genus %in% duplicates) %>%
  group_by(family, genus, major.clade) %>%
  summarise(n=length(species))  %>%
  arrange(genus, desc(n))

## We'll generally keep just the largest group:
keep <- as.data.frame(info)     %>%
  group_by(genus)               %>%
  summarise(family=family[[1]], major.clade=major.clade[[1]]) %>%
  arrange(genus)

## These are the species that we are ignoring (28).  This is tiny
## compared with the other taxonomic errors in the database.
drop <- as.data.frame(info) %>% group_by(genus) %>%
  summarise(family=family[-1], n=n[-1]) %>% arrange(genus)
sum(drop$n)

## But rewrite a couple of these anyway:
## * Malocarpus we have is not a cactus:
keep$family[keep$genus == "Malocarpus"] <- "Zygophyllaceae"
## * Ericaceae record is a tropicos error
keep$family[keep$genus == "Oreocallis"] <- "Proteaceae"
## * I camped under a Washingtonia filifera once (WCK)
keep$family[keep$genus == "Washingtonia"] <- "Arecaceae"

## Join this back in with the plant list, removing the duplicates.
tpl <- rbind(tpl[!(tpl$genus %in% duplicates),], keep)

## Assign order to the Plant List families using this lookup:
tpl$order <- lookup$order[match(tpl$family, lookup$family)]
## There are encoding issues here:
tpl$order[tpl$family == "Isoëtaceae"] <- "Isoetales"

## There are a handful of essentially unplaced families.  For now,
## these get their own pseudo-family
i <- tpl$order == ""
tpl$order[i] <- paste0("UnknownOrder-", tpl$family[i])

tpl <- tpl[c("genus", "family", "order", "major.clade")]

write.csv(tpl, "output/genus_order_lookup.csv", row.names=FALSE)
