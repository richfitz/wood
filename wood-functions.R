## # Load the woodiness data

## Load the data from the "forest" database.  This cleans up taxonomy,
## duplicated and nonstandard entries, empty records.  It then
## collapses things down to counts by genus, expanded to all known
## genera according to the plant list.  Final columns are
##   "genus", "family", "order", "W", "H", "K", "N", "p"
## where the last 5 columns are known woody, known herbaceous, known
## state, total species, and percentage of known species that are
## woody (i.e., W / K).
load.clean.data <- function(regenerate=FALSE) {
  filename <- "dat.g.rds"
  if ( !regenerate && file.exists(filename) )
    return(readRDS(filename))
  
  read.forest.csv <- function(filename)
    read.csv(file.path(path.forest, filename), stringsAsFactors=FALSE)

  ## Start by getting the woodiness information from the database
  dat <- read.forest.csv("export/speciesTraitData.csv")

  ## Only the columns we care about:
  dat <- data.frame(species=sub(" ", "_", dat$gs),
                    woodiness=dat$woodiness,
                    stringsAsFactors=FALSE)

  ## Filtered by whether or not they have woodiness information
  to.drop.wood.NA <- is.na(dat$woodiness)
  message(sprintf("Dropping %d species with NA woodiness values",
                  sum(to.drop.wood.NA)))
  dat <- dat[!to.drop.wood.NA,]

  ## Score the 633 species with no known information as NA
  to.drop.variable <- !(dat$woodiness %in% c("H", "W"))
  message(sprintf("Dropping %d species with variable woodiness",
                  sum(to.drop.variable)))
  dat <- dat[!to.drop.variable,]

  ## Next, normalise the species names.
  spp <-
    read.forest.csv("taxonomic/spermatophyta_synonyms_PLANTLIST.csv")
  spp$genus.synonym <- sub("_.+", "", spp$synonym)

  ## Locate the species names within the Plant List lookup table,
  ## dropping all species that do not exist.
  to.drop.no.name <- !dat$species %in% spp$synonym
  message(sprintf("Dropping %d species not in Plant List",
                  sum(to.drop.no.name)))
  dat <- dat[!to.drop.no.name,]

  ## Match the species names against the plant list
  idx <- match(dat$species, spp$synonym)

  ## About 90% of names are valid:
  message(sprintf("%2.1f%% of species have a valid name",
                  100*mean(spp$valid[idx])))

  ## Assign the correct/current species name and adding the genus
  dat$species <- spp$species[idx]
  dat$genus   <- spp$genus[idx]

  ## And look to see which species have now got duplicated records due
  ## to synonomy resolution:
  dups <- unique(sort(dat$species[duplicated(dat$species)]))
  message(sprintf("After synonym correction, %d duplicated entries",
                  length(dups)))

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
  to.drop.conflict <- is.na(dups.fixed$woodiness)
  message(sprintf("Dropping %d species due to trait confict after synonym",
                  sum(to.drop.conflict)))
  dups.fixed <- dups.fixed[!to.drop.conflict,]

  ## Drop the duplicated records from the original vector, and add in
  ## the resolved entries here:
  dat <- rbind(dat[-dups.i,], dups.fixed)

  ## So now 'dat' has species names sanitised to the same list that the
  ## plant list uses, and includes genus information.

  ## Next, compare the genus-> order lookup with the list of species
  ## that we have in the plant list.
  lookup <- read.forest.csv("taxonomic/genus_order_lookup.csv")
  lookup <- lookup[c("genus", "family", "order")]
  ## TODO: I've made this change here, because I want to confirm
  ## before changing the lookup table.
  lookup$order[lookup$family == "Adiantaceae"] <- "Polypodiales"

  ## There are a handful of essentially unplaced families.  For now,
  ## these get their own pseudo-family
  i <- lookup$order == ""
  lookup$order[i] <- paste0("UnknownOrder-", lookup$family[i])

  ## All of the genera with data are present in the lookup table:
  ## cases.
  if ( !all(dat$genus %in% spp$genus) ) # should *never* fail
    stop("Genera in data not known from Plant List")
  if ( !all(dat$genus %in% lookup$genus) )
    warning("Genera in data not known from lookup table",
            immediate.=TRUE)

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

  to.drop.no.order <- is.na(dat.g$order)
  if ( any(to.drop.no.order) ) {
    message(sprintf("Dropping %d genera (%d species, %d data) due to taxon fail",
                    sum(to.drop.no.order),
                    sum(dat.g$N[to.drop.no.order]),
                    sum(dat.g$K[to.drop.no.order])))
    dat.g <- dat.g[!to.drop.no.order,]
  }

  to.drop.no.family <- dat.g$family == ""
  if ( any(to.drop.no.family) ) {
    message(sprintf("Dropping %d genera (%d species, %d data) due to taxon fail (family)",
                    sum(to.drop.no.family),
                    sum(dat.g$N[to.drop.no.family]),
                    sum(dat.g$K[to.drop.no.family])))
    dat.g <- dat.g[!to.drop.no.family,]
  }
  rownames(dat.g) <- NULL

  ## This is OK for Pteridales, as it got synonomysed into other
  ## groups.  Where is it coming from though?
  to.drop.no.data.order <-
    which(tapply(dat.g$p, dat.g$order, function(x) all(is.nan(x))))
  if ( any(to.drop.no.data.order) ) {
    message(sprintf("Dropping %d orders because they have no data:\n\t%s",
                    length(to.drop.no.data.order),
                    paste(names(to.drop.no.data.order), collapse=", ")))
    dat.g <- dat.g[!(dat.g$order %in% names(to.drop.no.data.order)),]
  }

  message(sprintf("Final set: %d genera, %d with data, %d species known",
                  nrow(dat.g), sum(dat.g$K > 0), sum(dat.g$K)))

  ## Reorder columns
  cols <- c("genus", "family", "order", "W", "H", "K", "N", "p")
  dat.g <- dat.g[cols]

  saveRDS(dat.g, filename)
  dat.g
}

## # Load the survey data

load.survey <- function() {
  d <- read.csv(file="survey/Plant_survey_final.csv",
                stringsAsFactors=FALSE)
  names(d) <- c("Time", "Estimate", "Familiarity", "Training",
                "Continent", "Country")

  ## Drop the Time and Continent columns:
  d <- d[!(names(d) %in% c("Time", "Continent"))]

  ## Here are the different familiarity and training categories from
  ## "best" to "worst".
  lvl.familiarity <- c("Very Familiar", "Familiar", "Somewhat Familiar",
                       "What's a Plant?")
  d$Familiarity <- factor(d$Familiarity, lvl.familiarity, ordered=TRUE)
  
  lvl.training <-
    c("Postgraduate degree in botany or a related field",
      "Partially complete postgraduate degree in botany or a related field",
      "Undergraduate degree in botany or a related field",
      "Some botany courses at either an undergraduate or postgraduate level",
      "No formal training in botany")             
  d$Training <- factor(d$Training, lvl.training, ordered=TRUE)

  ## Standardise the country names:
  countries <- read.csv("survey/country_coords.csv",
                        stringsAsFactors=FALSE)
  d$Country <- cleanup.country.names(d$Country)
  
  idx <- match(d$Country, countries$Country)
  mssg <- na.omit(d$Country[is.na(idx)])
  if ( length(mssg) > 0 )
    warning("Dropped countries %s", paste(mssg, collapse=", "))
  d <- cbind(d, countries[idx,c("Long", "Lat")])
  d$Tropical <- abs(d$Lat) < 23 + 26/60

  rownames(d) <- NULL
  d
}

## Standardise the given country names into a common list.
cleanup.country.names <- function(x) {
  ## In cases where multiple countries are given, take the first one:
  x <- sub("( and |, | / | & ).+", "", x)
  ## Trim trailing non-alphabetic characters
  x <- sub("[^A-Za-z]+$", "", x)
  ## Translate inconsistent names:
  translate <- list(France="france",
                    "United States"=c("US", "USA"),
                    "United Kingdom"=c("UK", "Scotland"),
                    "Brazil"="Brasil")
  tr <- cbind(to=rep(names(translate), sapply(translate, length)),
              from=unlist(translate))
  i <- match(x, tr[,"from"])
  x[!is.na(i)] <- unname(tr[i[!is.na(i)],"to"])
  x[x == ""] <- NA
  x
}

## This generates the survey results with coordinates added.  Because
## it depends on rgdal, we don't run this very often.  In fact, it's
## not run automatically by wood.R, and the result of running it is
## under version control.
add.coordinates.to.survey <- function() {
  library(rgdal) # will cause error if not installed.
  country <- readOGR('survey/country/country.shp', 'country')
  
  coords <- as.data.frame(coordinates(country)[idx,])
  rownames(coords) <- NULL
  names(coords) <- c("Long", "Lat")

  d <- cbind(d, coords)

  write.csv(d, "survey/survey_results.csv", row.names=FALSE)
  invisible(TRUE)
}

## This downloads a bunch of data from gbif and gets latitude and
## longitude for all countries.  This is saved in the
## survey/country_coords.csv file, so does not need to be run except
## to refresh these coordinates (note that this file is under version
## control, so this is really here only as a reference of what was
## done).
build.country.list <- function() {
  ## Download files if they do not already exist.
  download.maybe <- function(url, dest) {
    if ( length(url) != 1 )
      stop("Scalar URL required")
    dest.file <- file.path(dest, basename(url))
    if ( !file.exists(dest.file) )
      download.file(url, dest.file)
    invisible(TRUE)
  }

  library(rgdal)
  dir.create("survey/country", FALSE)
  ext <- c("dbf", "fix", "ORG.dbf", "prj", "qix", "shp", "shx")
  urls <- paste0("http://ogc.gbif.org/data/data/shapefiles/country.",
                 ext)
  lapply(urls, download.maybe, "survey/country")

  country <- readOGR('survey/country/country.shp', 'country')
  ret <- as.data.frame(coordinates(country))
  names(ret) <- c("Long", "Lat")
  ret <- data.frame(Country=as.character(country@data$CNTRY_NAME),
                    ret)
  write.csv(ret, "survey/country_coords.csv", row.names=FALSE)
}

## # Order level phylogeny

build.order.tree <- function(dat.g, regenerate=FALSE) {
  filename <- "phy.o.rds"
  if ( !regenerate && file.exists(filename) ) {
    phy.o <- readRDS(filename)
  } else {
    mrca.tipset <- diversitree:::mrca.tipset
    drop.tip <- diversitree:::drop.tip.fixed

    phy <- read.tree("large-phylogeny.tre")

    ## Two phylogenetic errors in ferns need fixing:
    phy.order <- dat.g$order[match(sub("_.+$", "", phy$tip.label),
                                   dat.g$genus)]

    ## The Cyatheales and Polypodiales are hopelessly intertwined.
    ## I'm fixing this by dropping all Cyatheales except for
    ## Dicksonia_antarctica (arbitrarily).  According to APWEB, these
    ## groups are reciprocally monophyletic, but the tree does not
    ## resolve this and we didn't enforce it as a constraint.
    to.drop1 <-
      phy$tip.label[which(phy.order == "Cyatheales" &
                          phy$tip.label != "Dicksonia_antarctica")]

    ## Same with Ophioglossales, but using Ophioglossum_lusitanicum as
    ## the placeholder.
    to.drop2 <-
      phy$tip.label[which(phy.order == "Ophioglossales" &
                          phy$tip.label != "Ophioglossum_lusitanicum")]
    
    phy <- drop.tip(phy, c(to.drop1, to.drop2))

    ## Get the genera from this new tree.
    phy.genus <- sub("_.+$", "", phy$tip.label)

    missing.orders <- setdiff(dat.g$order, c(phy$node.label, ""))

    f <- function(x) {
      spp.x <- phy$tip.label[phy.genus %in% dat.g$genus[dat.g$order == x]]
      if ( length(spp.x) > 0 ) {
        node <- mrca.tipset(phy, spp.x)
        desc.x <- descendants.spp(node, phy)
        gen.x <- unique(phy.genus[match(desc.x, phy$tip.label)])
        ret <- unique(dat.g$order[dat.g$genus %in% gen.x])
        attr(ret, "node") <- node
        ret
      } else {
        character(0)
      }
    }

    tmp <- lapply(missing.orders, f)
    n <- sapply(tmp, length)

    ## These are going to be dropped.
    dropped.orders <- missing.orders[n == 0]
    dropped.orders <- dropped.orders[-grep("^Unknown", dropped.orders)]
    if ( length(dropped.orders) > 0 )
      warning(sprintf("Dropping orders: %s",
                      paste(dropped.orders, collapse=", ")))

    tmp.ok <- tmp[n == 1]
    nd <- sapply(tmp.ok, attr, "node")
    names(nd) <- sapply(tmp.ok, "[[", 1)
    n.tip <- length(phy$tip.label)
    i <- nd > n.tip
    phy$node.label[nd[i]] <- names(nd[i])
    phy$tip.label[nd[!i]] <- names(nd[!i])

    if ( any(n > 1) )
      warning("Dropping some orders!")

    nodes <- intersect(unique(dat.g$order), phy$node.label)
    to.keep <- sapply(match(nodes, phy$node.label), function(x)
                      descendants.spp(x, phy)[[1]])
    names(to.keep) <- nodes
    to.keep <- c(to.keep, structure(phy$tip.label[nd[!i]],
                                    names=names(nd[!i])))
    to.keep <- to.keep[-grep("UnknownOrder", names(to.keep))]

    to.drop <- setdiff(phy$tip.label, to.keep)
    phy.o <- drop.tip(phy, to.drop)
    phy.o$tip.label <- names(to.keep)[match(phy.o$tip.label, to.keep)]
    phy.o <- ladderize(phy.o)

    phy.o$n.taxa <- tapply(dat.g$N, dat.g$order, sum)[phy.o$tip.label]
    saveRDS(phy.o, filename)
  }
  
  phy.o
}

## # Utilities

## Evaluate expression 'expr' that produces a figure as a side effect,
## saving the result in a pdf file.
to.pdf <- function(filename, width, height, expr,
                   ..., pointsize=12, verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, width=width, height=height, pointsize=pointsize, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

## Add a label to a plot at a fixed relative location.
label <- function(px, py, lab, ..., adj=c(0, 1)) {
  lab <- LETTERS[lab]
  usr <- par("usr")
  text(usr[1] + px*(usr[2] - usr[1]),
       usr[3] + py*(usr[4] - usr[3]),
       lab, adj=adj, ...)
}

## Identify species descended from a node
descendants.spp <- function(node, phy) {
  i <- diversitree:::descendants(node, phy$edge)
  phy$tip.label[i[i <= length(phy$tip.label)]]
}

## Draw the outline of a histogram
hist.outline <- function(h, col, ..., density=TRUE) {
  dx <- diff(h$mids[1:2])
  xx <- rep(with(h, c(mids - dx/2, mids[length(mids)] + 
                      dx/2)), each = 2)
  yy <- c(0, rep(if (density) h$density else h$counts, each = 2), 0)
  lines(xx, yy, col = col, ...)
}
