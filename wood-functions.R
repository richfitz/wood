load.clean.data <- function(regenerate=FALSE) {
  filename <- "woodiness.rds"
  if ( !regenerate && file.exists(filename) )
    return(readRDS(filename))
  
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

  ## Locate the species names within the Plant List lookup table,
  ## dropping all species that do not exist.
  dat <- dat[dat$species %in% spp$synonym,]
  idx <- match(dat$species, spp$synonym) 

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

  ## Most duplicated species are of a single type (good) but a few
  ## aren't:
  f <- function(x)
    if ( length(unique(x)) == 1 ) x[1] else NA
  dups.i <- which(dat$species %in% dups)
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

  ## All of the genera with data are present in the lookup table:
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
  dat.g$N <- as.integer(spp.known[dat.g$genus])
  dat.g$p <- dat.g$W / dat.g$K

  ## Higher order taxonomy:
  idx <- match(dat.g$genus, lookup$genus)
  dat.g$family <- lookup$family[idx]
  dat.g$order <- lookup$order[idx]

  ## TODO: Hack for now.  what is actually up with these rows
  ## dat.g[which(is.na(dat.g$order)),]
  ##             genus W H K  N family order   p
  ## 1600    Benthamia 0 0 0 35   <NA>  <NA> NaN
  ## 7424 Lepidostemon 0 0 0  6   <NA>  <NA> NaN
  ## 8533     Monniera 0 0 0  2   <NA>  <NA> NaN
  dat.g <- dat.g[!is.na(dat.g$order),]

  rownames(dat.g) <- NULL

  saveRDS(dat.g, filename)
}
