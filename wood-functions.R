load.woodiness.data <- function() {
  readRDS("output/woodiness.rds")
}

load.woodiness.data.genus <- function(regenerate=FALSE,
                                      regenerate.raw=FALSE,
                                      extreme=FALSE) {
  if (identical(extreme, FALSE))
    filename <- "output/dat.g.rds"
  else
    filename <- sprintf("output/dat.g.%s.rds", substr(extreme, 1, 1))
  if (!regenerate && file.exists(filename))
    return(readRDS(filename))
  
  ## So now 'dat' has species names sanitised to the same list that the
  ## plant list uses, and includes genus information.
  dat <- load.woodiness.data(regenerate.raw)

  ## Next, compare the genus-> order lookup with the list of species
  ## that we have in the plant list.
  lookup1 <- read.csv("data/zae/genus_order_lookup.csv",
                      stringsAsFactors=FALSE)
  lookup2 <- read.csv("data/genus_order_lookup_extra.csv",
                      stringsAsFactors=FALSE)
  lookup <- rbind(lookup1, lookup2)
  rownames(lookup) <- NULL

  ## There are a handful of essentially unplaced families.  For now,
  ## these get their own pseudo-family
  i <- lookup$order == ""
  lookup$order[i] <- paste0("UnknownOrder-", lookup$family[i])

  ## All of the genera with data are present in the lookup table:
  ## cases.
  tpl <- read.csv("data/theplantlist/names_accepted.csv",
                  stringsAsFactors=FALSE)
  
  if (!all(dat$genus %in% tpl$genus)) # should *never* fail
    stop("Genera in data not known from Plant List")
  if (!all(dat$genus %in% lookup$genus))
    warning("Genera in data not known from lookup table",
            immediate.=TRUE)

  if (!identical(extreme, FALSE))
    dat$woodiness <-
      summarise.count(parse.count(dat$woodiness.count), extreme)

  ## Collapse these to get counts for all the genera that we know about.
  dat.g <- table(factor(dat$genus, sort(unique(tpl$genus))),
                 factor(dat$woodiness, c("W", "variable", "H")))
  dat.g <- data.frame(genus=rownames(dat.g),
                      W=dat.g[,"W"],
                      V=dat.g[,"variable"],
                      H=dat.g[,"H"],
                      stringsAsFactors=FALSE)
  ## Assuming we drop all 
  dat.g$K <- dat.g$H + dat.g$W
  dat.g$p <- dat.g$W / dat.g$K
  
  ## Include the counts of known species:
  spp.known <- table(tpl$genus)
  dat.g$N <- as.integer(spp.known[dat.g$genus])

  ## Higher order taxonomy:
  idx <- match(dat.g$genus, lookup$genus)
  dat.g$family <- lookup$family[idx]
  dat.g$order <- lookup$order[idx]

  k <- dat.g$W + dat.g$H
  message(sprintf("Final set: %d genera, %d with data, %d species known",
                  nrow(dat.g), sum(k > 0), sum(k)))

  ## Reorder columns
  cols <- c("genus", "family", "order", "W", "V", "H", "N", "K", "p")
  dat.g <- dat.g[cols]
  rownames(dat.g) <- NULL

  dat.g <- dat.g[order(dat.g$order, dat.g$family, dat.g$genus),]

  saveRDS(dat.g, filename)
  dat.g
}

## Check the classification by pulling apart the count.
parse.count <- function(x) {
  res <- t(sapply(strsplit(x, ";", fixed=TRUE), as.integer))
  colnames(res) <- c("H", "V", "W")
  drop(res)
}
summarise.count <- function(x, extreme=FALSE) {
  if (identical(extreme, FALSE)) {
    ans <- ifelse(x[,"W"] > x[,"H"], "W", "H")
    ans[(x[,"W"] == 0 & x[,"H"] == 0 & x[,"V"] > 0) |
        x[,"W"] == x[,"H"]] <- "variable"
  } else if (identical(extreme, "woody")) {
    ans <- ifelse(x[,"W"] > 0 | x[,"V"] > 0, "W", "H")
  } else if (identical(extreme, "herbaceous")) {
    ans <- ifelse(x[,"H"] > 0 | x[,"V"] > 0, "H", "W")
  } else {
    stop("Invalud argument for 'extreme'")
  }
  ans
}


load.survey <- function() {
  d <- read.csv(file="data/survey_results.csv",
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
  countries <- read.csv("data/geo/country_coords.csv",
                        stringsAsFactors=FALSE)
  d$Country <- cleanup.country.names(d$Country)
  
  idx <- match(d$Country, countries$Country)
  mssg <- na.omit(d$Country[is.na(idx)])
  if (length(mssg) > 0)
    warning("Dropped countries %s", paste(mssg, collapse=", "))
  d <- cbind(d, countries[idx,c("Long", "Lat")])
  d$Tropical <- abs(d$Lat) < 23 + 26/60

  rownames(d) <- NULL
  d
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

  write.csv(d, "data/survey_results.csv", row.names=FALSE)
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
    if (length(url) != 1)
      stop("Scalar URL required")
    dest.file <- file.path(dest, basename(url))
    if (!file.exists(dest.file))
      download.file(url, dest.file)
    invisible(TRUE)
  }

  library(rgdal)
  path.raw <- "data/geo/raw"
  dir.create(path.raw, showWarnings=FALSE, recursive=TRUE)
  ext <- c("dbf", "fix", "ORG.dbf", "prj", "qix", "shp", "shx")
  urls <- paste0("http://ogc.gbif.org/data/data/shapefiles/country.",
                 ext)
  lapply(urls, download.maybe, path.raw)

  country <- readOGR(file.path(path.raw, "country.shp"), 'country')
  ret <- as.data.frame(coordinates(country))
  names(ret) <- c("Long", "Lat")
  ret <- data.frame(Country=as.character(country@data$CNTRY_NAME),
                    ret)
  write.csv(ret, "data/geo/country_coords.csv", row.names=FALSE)
}

## # Order level phylogeny

build.order.tree <- function(dat.g, regenerate=FALSE) {
  filename <- "output/phy.o.rds"
  if (!regenerate && file.exists(filename)) {
    phy.o <- readRDS(filename)
  } else {
    mrca.tipset <- diversitree:::mrca.tipset
    drop.tip <- diversitree:::drop.tip.fixed

    phy <- read.tree("data/zae/Vascular_Plants_rooted.dated.tre")

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
      if (length(spp.x) > 0) {
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
    if (length(dropped.orders) > 0)
      warning(sprintf("Dropping orders: %s",
                      paste(dropped.orders, collapse=", ")))

    tmp.ok <- tmp[n == 1]
    nd <- sapply(tmp.ok, attr, "node")
    names(nd) <- sapply(tmp.ok, "[[", 1)
    n.tip <- length(phy$tip.label)
    i <- nd > n.tip
    phy$node.label[nd[i]] <- names(nd[i])
    phy$tip.label[nd[!i]] <- names(nd[!i])

    if (any(n > 1))
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

fig.fraction.on.phylogeny <- function(phy.o, res) {
  ## Higher level taxonomy
  hlt <- read.csv("data/high-level-taxonomy.csv", stringsAsFactors=FALSE)
  phy.group <- hlt$Group[match(phy.o$tip.label, hlt$Order)]
  tmp <- 
    lapply(seq_len(max(phy.o$edge)), function(x)
           if ( x <= length(phy.o$tip.label) ) phy.group[[x]] else
           unique(phy.group[get.descendants(x, phy.o, TRUE)]))
  grp <- sapply(tmp, function(x) if (length(x) == 1) x else "Rest")

  col <- unname(cols.tree[grp])
  col2 <- col[match(phy.o$edge[,2], seq_along(grp))]

  p <- structure(res$order[["p.mean"]], names=res$order$order)
  p <- p[phy.o$tip.label]

  t <- max(branching.times(phy.o))
  offset <- .15

  op <- par(no.readonly=TRUE)
  on.exit(par(op))

  ## Drop orders with < 100 species, execpt for a couple we can fit
  ## in.
  drop <- c("Isoetales",
            "Psilotales",
            "Ophioglossales",
            "Equisetales",
            "Osmundales",
            "Salviniales",
            "Ginkgoales",
            "Welwitschiales",
            "Gnetales",
            "Ephedrales",
            "Nymphaeales",
            "Amborellales",
            "Austrobaileyales",
            "Chloranthales",
            "Ceratophyllales",
            # "Canellales",    # keep
            # "Acorales",      # keep
            # "Petrosaviales", # keep
            "Trochodendrales", # drop } could keep 1/2
            "Gunnerales",      # drop }
            "Crossosomatales", # drop
            "Picramniales",    # drop
            "Huerteales",      # drop
            "Berberidopsidales", # drop (borderline)
            "Paracryphiales", # drop
            "Escalloniales",  # drop
            "Bruniales",      # drop
            #"Garryales"      # keep
            "Schizeales"       # also drop
            )
  tip.color <- ifelse(phy.o$tip.label %in% drop, "#ffffff00", "black")
  plt <- 
    diversitree:::plot2.phylo(phy.o, type="fan", cex=.5, no.margin=TRUE,
                              label.offset=t * .15, font=1,
                              edge.col=col2, tip.color=tip.color,
                              n.taxa=sqrt(phy.o$n.taxa)/10)
  xy <- plt$xy

  r <- max(xy$r)*(1+offset)
  n.tip <- length(phy.o$tip.label)
  xy <- plt$xy[seq_len(n.tip),]
  xy.lab <- data.frame(x=cos(xy$theta)*r,
                       y=sin(xy$theta)*r)
  xrad <- .5 * diff(par("usr")[1:2])/50
  pie <- cbind(p, 1 - p)
  pie.col <- cols.woody

  r <- 3/4
  r0 <- max(xy$r) * (1 + offset * (1-r)/2)
  r2 <- max(xy$r) * (1 + offset * (1 - (1-r)/2))
  r1 <- r0 * p + r2 * (1-p)

  w <- 3

  xx1 <- c(rbind(r0 * cos(xy$theta) + w * cos(xy$theta + pi/2),
                 r0 * cos(xy$theta) - w * cos(xy$theta + pi/2),
                 r1 * cos(xy$theta) - w * cos(xy$theta + pi/2),
                 r1 * cos(xy$theta) + w * cos(xy$theta + pi/2),
                 NA))
  yy1 <- c(rbind(r0 * sin(xy$theta) + w * sin(xy$theta + pi/2),
                 r0 * sin(xy$theta) - w * sin(xy$theta + pi/2),
                 r1 * sin(xy$theta) - w * sin(xy$theta + pi/2),
                 r1 * sin(xy$theta) + w * sin(xy$theta + pi/2),
                 NA))
  
  xx2 <- c(rbind(r2 * cos(xy$theta) + w * cos(xy$theta + pi/2),
                 r2 * cos(xy$theta) - w * cos(xy$theta + pi/2),
                 r1 * cos(xy$theta) - w * cos(xy$theta + pi/2),
                 r1 * cos(xy$theta) + w * cos(xy$theta + pi/2),
                 NA))
  yy2 <- c(rbind(r2 * sin(xy$theta) + w * sin(xy$theta + pi/2),
                 r2 * sin(xy$theta) - w * sin(xy$theta + pi/2),
                 r1 * sin(xy$theta) - w * sin(xy$theta + pi/2),
                 r1 * sin(xy$theta) + w * sin(xy$theta + pi/2),
                 NA))

  polygon(xx1, yy1, border="black", col=pie.col[2], lwd=.3)
  polygon(xx2, yy2, border="black", col=pie.col[1], lwd=.3)

  cex.legend <- 2/3
  str <- str.leg <- setdiff(names(cols.tree), "Rest") # drop backbone
  str.leg[str.leg == "BasalAngiosperms"] <- '"Basal Angiosperms"'
  legend("topleft", str.leg, fill=cols.tree[str],
         cex=cex.legend, bty="n", border=NA)
  legend("topright", names(cols.woody), fill=cols.woody,
         cex=cex.legend, bty="n", border="black")
}

## # Utilities

## Evaluate expression 'expr' that produces a figure as a side effect,
## saving the result in a pdf file.
to.pdf <- function(filename, width, height, expr,
                   ..., pointsize=12, verbose=TRUE) {
  if (verbose)
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, width=width, height=height, pointsize=pointsize, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

## Add a label to a plot at a fixed relative location.
label <- function(px, py, lab, ..., adj=c(0, 1)) {
  lab <- LETTERS[lab]
  usr <- par("usr")
  x <- usr[1] + px*(usr[2] - usr[1])
  y <- usr[3] + py*(usr[4] - usr[3])
  if (par("xlog")) x <- 10^x
  if (par("ylog")) y <- 10^y
  text(x, y, lab, adj=adj, ...)
}

## Identify species descended from a node
descendants.spp <- function(node, phy) {
  i <- diversitree:::descendants(node, phy$edge)
  phy$tip.label[i[i <= length(phy$tip.label)]]
}

## Draw the outline of a histogram
hist.outline <- function(h, ..., density=TRUE) {
  xy <- hist.xy(h, density)
  lines(xy, ...)
}
hist.fill <- function(h, ..., density=TRUE) {
  xy <- hist.xy(h, density)
  polygon(xy, ...)
}

hist.xy <- function(h, density=TRUE) {
  dx <- diff(h$mids[1:2])
  xx <- rep(with(h, c(mids - dx/2, mids[length(mids)] + 
                      dx/2)), each = 2)
  yy <- c(0, rep(if (density) h$density else h$counts, each = 2), 0)
  list(x=xx, y=yy)
}

mix <- function(cols, col2, p) {
  m <- col2rgb(cols)
  m2 <- col2rgb(rep(col2, length=length(cols)))
  m3 <- (m * p + m2 * (1-p))/255
  rgb(m3[1,], m3[2,], m3[3,])
}

log.seq.range <- function(x, n)
  exp(seq(log(min(x)), log(max(x)), length=n))
