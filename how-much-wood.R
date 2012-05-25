## % How much of the world is woody, first analysis.
library(multicore)
library(diversitree)
library(RColorBrewer)
path.forest <- readLines("~/.forest_path")
source("functions.R")

read.forest.csv <- function(filename)
  read.csv(file.path(path.forest, filename), stringsAsFactors=FALSE)

dat <- read.forest.csv("export/speciesTraitData.csv")

## Extract the woodiness data:
dat$woodiness[!(dat$woodiness %in% c("H", "W")) &
              !is.na(dat$woodiness)] <- NA
x <- factor(dat$woodiness)
names(x) <- dat$gs
x <- x[!is.na(x)]

info <-
  read.forest.csv("taxonomic/spermatophyta_richnesses.plantlist.APGIIItax.csv")

## Determine orders:
lookup <- read.forest.csv("taxonomic/genus_order_lookup.csv")

## Currently, some orders are missing for individual genera that
## belong to families where the order relationship is already known.
## Fill these in:
i <- lookup$order_final == "" & lookup$family_final != ""
tmp <- lookup[lookup$order_final != "",]
res <- tmp$order_final[match(lookup$family_final[i],
                             tmp$family_final)]
res[is.na(res)] <- ""
lookup$order_final[i] <- res

## Add order level to info, which also contains per-genus counts.
info$order <- lookup$order_final[match(info$family, lookup$family_final)]
info$order[is.na(info$order) | info$order == ""] <- "UnknownOrder"

## Extract genus for each observation to match against the 'info' table.
genus <- sub(" .+", "", names(x))

## Filter for known-ness (there is a lookup issue here to resolve,
## dropping about 2400 species).  This filters both the genus names
## and the data.
ok <- genus %in% info$genus
x <- x[ok]
genus <- genus[ok]
genus <- factor(genus, levels=sort(unique(info$genus)))

## 'dat.g' is grouped by genus, and contains known numbers of species
## (N), and numbers known to be woody (W) and herbacious (H).
dat.g <- as.data.frame(unclass(table(genus, x)))
idx <- match(rownames(dat.g), info$genus)
dat.g <- cbind(order=info$order[idx],
               family=info$family[idx],
               genus=rownames(dat.g),
               N=info$species[idx],
               dat.g)
dat.g <- dat.g[order(dat.g$order),]
rownames(dat.g) <- NULL

## For some species, we know about more states than we know about
## species; increase the species count in this case.
fix <- with(dat.g, N < H + W)
dat.g$N[fix] <- with(dat.g[fix,], H + W)

## Finally, estimate the fraction of known states are woody:
dat.g$p <- with(dat.g, W / (H + W))

## Here is the simulator.  For each genus, sample unknown species from
## a hypergeometric distribution with parameters sampled from the
## known number of woody, herbacious and known species.  With this
## sampled data, sample proportions of woody species for the genera
## with nothing known from the empirical distribution of woodiness,
## and sample woody/herbacious species from a binomial.
sim <- function(x, nrep) {
  i <- is.na(x$p)
  w <- matrix(NA, nrow(x), nrep)
  n.unk <- sum(i)
  w[!i,] <- t(sapply(which(!i), function(j)
                     rhyper2(nrep, x$H[j], x$W[j], x$N[j])))
  w[i,] <- apply(w[!i,,drop=FALSE] / x$N[!i], 2, function(y)
                 rbinom(n.unk, x$N[i], quantile(y, runif(n.unk))))

  ## This is surprisingly convoluted:
  fam <- as.character(x$family)
  ret <- t(sapply(split(as.data.frame(w), fam), colSums))
  colnames(ret) <- NULL
  ret
}

## Sample the fraction 1000 times.  This takes a while - currently
## about 1.5 minutes.
##+ simulate,cache=TRUE
nrep <- 1000
dat.g.split <- split(dat.g, dat.g$order)
simulated.counts <- mclapply(dat.g.split, sim, nrep)

## Next, build collapsed sets of data by family and by order.
tot.f <- do.call(rbind, simulated.counts)
tot.o <- do.call(rbind, lapply(simulated.counts, colSums))

## Number of species per family and order:
n.f <- tapply(dat.g$N, dat.g$family, sum)[rownames(tot.f)]
n.o <- tapply(dat.g$N, dat.g$order, sum)[rownames(tot.o)]

## From this, compute the per-family and per-order fraction
prop.f <- tot.f / rep(n.f, nrep)
prop.o <- tot.o / rep(n.o, nrep)
prop.all <- colSums(tot.o) / sum(n.o)

## And the mean fraction woody per family and order:
p.f <- rowMeans(prop.f)
p.o <- rowMeans(prop.o)
p.all <- mean(prop.all)

## Here is a function that converts a value on 0..1 to increasingly
## dark blues.
cols <- function(x) {
  tmp <- colorRamp(brewer.pal(9, "Blues"))(x)
  rgb(tmp[,1], tmp[,2], tmp[,3], max=255)
}

## Histogram of woodiness
##+ fig.cap="Estimated fraction of woodiness over vascular plants"
par(mar=c(4.1, 4.1, .5, .5))
hist(prop.all, xlab="Percent woody", ylab="",
     yaxt="n", bty="n", main="", col=cols(p.all))

## Plot the tree with a distribution of fractions of woodiness per
## family around the outside.
phy.f <- get(load("phy.f.Rdata"))

## This tree has order information:
phy.f.ord <- phy.f$class

## Drop tips with no woodiness estimates (due to taxonomic
## differences)
to.drop <- setdiff(phy.f$tip.label, names(p.f))
phy.f <- diversitree:::drop.tip.fixed(phy.f, to.drop)

## And sort the order information by the current taxon labels:
phy.f.ord <- unname(phy.f.ord[phy.f$tip.label])


## Identify the nodes corresponding to orders, and plot the internal
## ones.
##+ fig.cap="Estimated per-family and per-order woodiness proportions"
plt <- trait.plot.cont(phy.f, p.f, list(cols), class=phy.f.ord,
                       w=1/30, font=1)
mrca.tipset <- diversitree:::mrca.tipset
nd <- sapply(sort(unique(phy.f.ord)), function(x)
             mrca.tipset(phy.f, phy.f$tip.label[phy.f.ord==x &
                                                !is.na(phy.f.ord)]))
nd.int <- nd[nd > length(phy.f$tip.label)]
points(plt$xy$xx[nd.int], plt$xy$yy[nd.int], pch=19, cex=2,
       col=cols(p.o[names(nd.int)]))

