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

  ## TODO: Hack for now.  what is actually up with these rows
  ## dat.g[which(is.na(dat.g$order)),]
  ##             genus W H K  N family order   p
  ## 1600    Benthamia 0 0 0 35   <NA>  <NA> NaN
  ## 7424 Lepidostemon 0 0 0  6   <NA>  <NA> NaN
  ## 8533     Monniera 0 0 0  2   <NA>  <NA> NaN
  to.drop.no.order <- is.na(dat.g$order)
  message(sprintf("Dropping %d genera (%d species, %d data) due to taxon fail",
                  sum(to.drop.no.order),
                  sum(dat.g$N[to.drop.no.order]),
                  sum(dat.g$K[to.drop.no.order])))
  dat.g <- dat.g[!to.drop.no.order,]
  rownames(dat.g) <- NULL

  to.drop.no.family <- dat.g$family == ""
  message(sprintf("Dropping %d genera (%d species, %d data) due to taxon fail (family)",
                  sum(to.drop.no.family),
                  sum(dat.g$N[to.drop.no.family]),
                  sum(dat.g$K[to.drop.no.family])))
  dat.g <- dat.g[!to.drop.no.family,]
  rownames(dat.g) <- NULL

  to.drop.no.data.order <-
    which(tapply(dat.g$p, dat.g$order, function(x) all(is.nan(x))))
  warning(sprintf("Dropping %d orders because they have no data:\n\t%s",
                  length(to.drop.no.data.order),
                  paste(names(to.drop.no.data.order), collapse=", ")),
          immediate.=TRUE)
  dat.g <- dat.g[!(dat.g$order %in% names(to.drop.no.data.order)),]

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

load.survey <- function() {
  d <- read.csv(file="survey/Plant_survey_final.csv",
                stringsAsFactors=TRUE)
  ## Remove timestamp column and continent column:
  d <- d[-c(1,5)]

  ## change the colnames
  colnames(d) <- c("Estimate", "Familiarity", "Training", "Country")

  ## Here are the different familiarity and training categories from
  ## "best" to "worst".
  lvl.familiarity <- c("Very Familiar", "Familiar", "Somewhat Familiar",
                       "What's a Plant?")
  lvl.training <-
    c("Postgraduate degree in botany or a related field",
      "Partially complete postgraduate degree in botany or a related field",
      "Undergraduate degree in botany or a related field",
      "Some botany courses at either an undergraduate or postgraduate level",
      "No formal training in botany")             

  d$Familiarity <- factor(d$Familiarity, lvl.familiarity, ordered=TRUE)
  d$Training <- factor(d$Training, lvl.training, ordered=TRUE)

  d
}
