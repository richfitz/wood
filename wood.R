## % How much of the world is woody, first analysis.
library(multicore)
library(diversitree)
library(RColorBrewer)
path.forest <- readLines("~/.forest_path")
source("build-family-tree.R")

source("wood-functions.R")
## Data at the level of genus: has taxonomic information and counts of
## number of species that are woody, herbaceous, known, and known to
## exist.
dat.g <- load.clean.data()

## Here is the simulator.  For each genus, sample unknown species from
## a hypergeometric distribution with parameters sampled from the
## known number of woody, herbacious and known species.  With this
## sampled data, sample proportions of woody species for the genera
## with nothing known from the empirical distribution of woodiness,
## and sample woody/herbacious species from a binomial.
sim <- function(x, nrep, mode="replacement") {
  i <- is.na(x$p)
  w <- matrix(NA, nrow(x), nrep)
  n.unk <- sum(i)
  if ( mode == "replacement" ) {
    w[!i,] <- t(sapply(which(!i), function(j)
                       rhyper2(nrep, x$H[j], x$W[j], x$N[j])))
  } else {
    j <- !i
    w[j,] <- x$W[j] + rbinom(sum(j), x$N[j]-x$K[j], x$W[j]/x$K[j])
  }
  w[i,] <- apply(w[!i,,drop=FALSE] / x$N[!i], 2, function(y)
                 rbinom(n.unk, x$N[i], quantile(y, runif(n.unk))))

  ## This is surprisingly convoluted:
  fam <- as.character(x$family)
  ret <- t(sapply(split(as.data.frame(w), fam), colSums))
  colnames(ret) <- NULL
  ret
}

## Sample from the hyper geometric distribution.  For R's rhyper, we
## have to provide the number of white and black balls from the urn,
## and it samples them.  However, we don't know what this is.

## What this does is samples the (estimated) true number of white
## balls 'nn' times and returns that.
rhyper2 <- function(nn, s0, s1, xn, fraction=FALSE) {
  x1 <- seq(s1, xn - s0)
  x0 <- xn - x1
  p1 <- dhyper(s1, x1, x0, s0+s1)
  p1 <- p1 / sum(p1)
  x1[sample(length(p1), nn, TRUE, p1)]
}

## Sample the fraction 1000 times.  This takes a while - currently
## about 1.5 minutes.
##+ simulate,cache=TRUE
nrep <- 1000
dat.g.split <- split(dat.g, dat.g$order)
## simulated.counts <- lapply(dat.g.split, sim, nrep)

## Drop one problem case entirely for now:
nok <- which(sapply(dat.g.split, function(x) all(is.nan(x$p))))
dat.g.split <- dat.g.split[-nok]

## "Weak prior" sampling with replacement
sim.weak <- lapply(dat.g.split, sim, nrep)
## "Strong prior" sampling based on the number of species
sim.strong <- lapply(dat.g.split, sim, nrep, "no-replacement")

## RGF: There is a bunch of repetition here and I'll try and clean it
## up soon.

## Next, build collapsed sets of data by family and by order.
tot.weak.f <- do.call(rbind, sim.weak)
tot.weak.o <- do.call(rbind, lapply(sim.weak, colSums))

tot.strong.f <- do.call(rbind, sim.strong)
tot.strong.o <- do.call(rbind, lapply(sim.strong, colSums))

## Number of species per family and order:
n.f <- tapply(dat.g$N, dat.g$family, sum)[rownames(tot.strong.f)]
n.o <- tapply(dat.g$N, dat.g$order, sum)[rownames(tot.strong.o)]

## From this, compute the per-family and per-order fraction
prop.weak.f <- tot.weak.f / rep(n.f, nrep)
prop.weak.o <- tot.weak.o / rep(n.o, nrep)
prop.weak.all <- colSums(tot.weak.o) / sum(n.o)

prop.strong.f <- tot.strong.f / rep(n.f, nrep)
prop.strong.o <- tot.strong.o / rep(n.o, nrep)
prop.strong.all <- colSums(tot.strong.o) / sum(n.o)

## And the mean fraction woody per family and order:
p.weak.f <- rowMeans(prop.weak.f)
p.weak.o <- rowMeans(prop.weak.o)
p.weak.all <- mean(prop.weak.all)

p.strong.f <- rowMeans(prop.strong.f)
p.strong.o <- rowMeans(prop.strong.o)
p.strong.all <- mean(prop.strong.all)

## Now, look at the distributions of woodiness among families:
label <- function(px, py, lab, ..., adj=c(0, 1)) {
  usr <- par("usr")
  text(usr[1] + px*(usr[2] - usr[1]),
       usr[3] + py*(usr[4] - usr[3]),
       lab, adj=adj, ...)
}

## RGF: I will tidy this tomorrow.
pdf("doc/fraction-by-genus.pdf", height=6, width=6)
par(mfrow=c(2, 1), mar=c(2, 1, .5, .5), oma=c(2, 0, 0, 0))

tmp <- dat.g$p[dat.g$K >= 10] # genera with 10 records
h <- hist(tmp, 50, main="", xaxt="n", yaxt="n", lab="")
mtext("Probability density", 2)
axis(1, tick=TRUE, label=FALSE)
label(.02, .96, "a)")
    
cols <- c("#ff0000ff", "#0000ff66")
hist(p.strong.f, col=cols[1], n=50, freq=FALSE, main="", yaxt="n",
     ylab="", xlab="% woody species in genus")
mtext("Probability density", 2)
hist(p.weak.f, col=cols[2], n=50, freq=FALSE, add=TRUE)
legend("topleft", c("No replacement (strong prior)",
                    "Replacement (weak prior)"),
       fill=cols, bty="n", cex=.75, inset=c(.1, 0))
label(.02, .96, "b)")
dev.off()
