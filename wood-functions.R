load.clean.data <- function(regenerate=FALSE) {
  filename <- "woodiness.rds"
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

  saveRDS(dat.g, filename)
  dat.g
}

to.pdf <- function(filename, width, height, expr,
                   ..., pointsize=12, cairo=FALSE, verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  if ( cairo ) require(Cairo)
  dev <- if ( cairo ) CairoPDF else pdf
  dev(filename, width=width, height=height, pointsize=pointsize, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

label <- function(px, py, lab, ..., adj=c(0, 1)) {
  usr <- par("usr")
  text(usr[1] + px*(usr[2] - usr[1]),
       usr[3] + py*(usr[4] - usr[3]),
       lab, adj=adj, ...)
}

## This will roll into trait.plot soon.
trait.plot.cont <- function(tree, dat, cols, lab=names(cols), str=0:1,
                            class=NULL, type="f", w=1/50,
                            legend=length(cols) > 1, cex.lab=.5,
                            font.lab=3, cex.legend=.75, margin=1/4,
                            check=TRUE, quiet=FALSE) {
  if ( type != "f" )
    stop("type != f not yet implemented")
  if ( !is.null(class) && length(class) != length(tree$tip.label) )
    stop("'class' must be a vector along tree$tip.label")
  n <- length(cols)
  if ( n < 1 )
    stop("Need some colours")
  if ( !is.data.frame(dat) ) {
    if ( is.vector(dat) && n == 1 ) {
      nm <- names(dat)
      dat <- matrix(dat)
      rownames(dat) <- nm
    } else {
      stop("dat must be a matrix")
    }
  }
  if ( !all(tree$tip.label %in% rownames(dat)) )
    stop("All taxa must have entries in 'dat' (rownames)")
  if ( n > 1 ) {
    if ( !all(names(cols) %in% names(dat)) )
      stop("Not all colours have data")
    if ( is.null(names(cols)) )
      stop("'cols' must be named")
    dat <- dat[names(cols)]
  }

  dat <- dat[tree$tip.label,,drop=FALSE]

  par(mar=rep(0, 4))
  t <- max(branching.times(tree))
  w <- w * t

  plot2.phylo <- diversitree:::plot2.phylo
  group.label.tip <- diversitree:::group.label.tip
  filled.arcs <- diversitree:::filled.arcs
  if ( is.null(class) ) {
    plt <- plot2.phylo(tree, type="f", show.tip.label=TRUE,
                       label.offset=(n+2)*w, cex=cex.lab)
  } else {
    plt <- plot2.phylo(tree, type="f", show.tip.label=FALSE,
                       label.offset=t*margin)
    group.label.tip(plt, class, "black", "black",
                    offset.bar=w*(n+2), offset.lab=w*(n+3), lwd=1.5,
                    cex=cex.lab, font=font.lab,
                    check=check, quiet=quiet)
  }

  xy <- plt$xy
  theta <- xy$theta[seq_along(tree$tip.label)]
  dt <- diff(sort(theta))[1]/2

  for ( i in seq_along(cols) ) {
    filled.arcs(theta - dt, theta + dt, max(xy$x) + i * w, w,
                cols[[i]](dat[,i]))
  }
  invisible(plt)
}

## Here is a function that converts a value on 0..1 to increasingly
## dark blues.
make.col.function <- function(cols) {
  ramp <- colorRamp(cols)  
  function(x) {
    i <- !is.na(x)
    tmp <- ramp(x[i])
    ret <- rep(NA_character_, length(x))
    ret[i] <- rgb(tmp[,1], tmp[,2], tmp[,3], max=255)
    ret
  }
}

download.maybe <- function(url, dest) {
  if ( length(url) != 1 )
    stop("Scalar URL required")
  dest.file <- file.path(dest, basename(url))
  if ( !file.exists(dest.file) )
    download.file(url, dest.file)
  invisible(TRUE)
}

build.country.list <- function() {
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


## To build family-level tree:
get.class <- function(phy.f) {
  lab <- phy.f$node.label
  keep <- grepl("ales$", lab)
  keep[lab == "Monilophyte"] <- TRUE
  keep[lab %in% c("Gleicheniales", "Schizaeales", "Salviniales")] <- FALSE
  keep[grep("_To_", lab)] <- FALSE
  lab[!keep] <- NA

  descendants <- diversitree:::descendants
  desc.spp <- function(node, phy) {
    n.spp <- length(phy$tip.label)
    ret <- descendants(match(node, phy$node.label) + n.spp, phy$edge)
    phy$tip.label[ret[ret <= n.spp]]
  }
  orders <- lab[!is.na(lab)]
  tmp <- lapply(orders, desc.spp, phy.f)
  class <- data.frame(order=rep(orders, sapply(tmp, length)),
                      family=unlist(tmp), stringsAsFactors=FALSE)
  any(duplicated(unlist(class$family)))

  ## There are a handful of otherws potentially worth keeping, too.
  grp <- class$order[match(phy.f$tip.label, class$family)]
  grp[phy.f$tip.label == "Arecaceae"] <- "Arecacales"
  grp[phy.f$tip.label == "Boraginaceae"] <- "Boraginaceae"
  grp[phy.f$tip.label == "Vitaceae"] <- "Vitales"
  grp[phy.f$tip.label == "Dilleniaceae"] <- "Dilleniaceae"
  grp[phy.f$tip.label == "Escalloniaceae"] <- "Escalloniaceae"
  grp[phy.f$tip.label == "Sabiaceae"] <- "Sabiaceae"
  grp[phy.f$tip.label == "Chloranthaceae"] <- "Chloranthaceae"
  grp[phy.f$tip.label == "Buxaceae"] <- "Buxaceae"
  grp[phy.f$tip.label == "Icacinaceae"] <- "Icacinaceae"
  grp[phy.f$tip.label == "Acoraceae"] <- "Acoraceae"
  grp[phy.f$tip.label == "Paracryphiaceae"] <- "Paracryphiaceae"
  grp
}

build.family.tree <- function(regenerate=FALSE) {
  filename <- "phy.f.rds"
  if ( !regenerate && file.exists(filename) ) {
    phy.f <- readRDS(filename)
  } else {
    browser()
    phy <- read.tree(file.path(path.forest,
                               "taxonomic/trees/spLevelApgBackbone.tre"))
    i <- phy$node.label == ""
    phy$node.label[i] <- sprintf("node.%d", (1:phy$Nnode)[i])
    problems <-
      readLines(file.path(path.forest,
                          "taxonomic/trees/problemSpeciesAPG.txt"))
    problems <- sub(" *$", "", problems) # kill trailing whitespace
    phy <- diversitree:::drop.tip.fixed(phy, problems)
    phy <- compute.brlen(phy, method="Grafen", power=.45)
    families <-
      read.csv(file.path(path.forest, "taxonomic/genus_order_lookup.csv"),
               stringsAsFactors=FALSE)
    genus <- sub("_.+$", "", phy$tip.label)
    family <- families$family[match(genus, families$genus)]
    family[is.na(family)] <- genus[is.na(family)] # Fabaceae x 5
    phy.f <- clades.from.classification(phy, family, check=FALSE)
    ord <- get.class(phy.f)
    names(ord) <- phy.f$tip.label
    phy.f$class <- ord
    saveRDS(phy.f, filename)
  }
  phy.f
}

descendants.spp <- function(node, phy) {
  i <- diversitree:::descendants(node, phy$edge)
  phy$tip.label[i[i <= length(phy$tip.label)]]
}

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
    if (  length(dropped.orders) > 0 )
      warning("Dropping orders: %s", paste(dropped.orders,
                                           collapse=", "))

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

hist.outline <- function(h, col, ..., density=TRUE) {
  dx <- diff(h$mids[1:2])
  xx <- rep(with(h, c(mids - dx/2, mids[length(mids)] + 
                      dx/2)), each = 2)
  yy <- c(0, rep(if (density) h$density else h$counts, each = 2), 0)
  lines(xx, yy, col = col, ...)
}
