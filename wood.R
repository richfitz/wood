library(diversitree)
source("wood-functions.R")

## Suppress a warning about incompatibility with results from R < 2.2.0
invisible(suppressWarnings(sample(1:250, 1, pr=rep(1, 250), replace=TRUE)))

## We're going to put some partly processed data here.
dir.create("output", FALSE)

## Colours used throughout:
cols.methods <- c(binomial="#a63813",       # red
                  hypergeometric="#4d697f") # blue
cols.tree <- c(Monilophytes="#a63813",     # reddish brown
               Gymnosperms="#21313b",      # dark brown
               BasalAngiosperms="#eeb911", # yellow
               Monocots="#204d14",         # green
               Eudicots="#4d697f",         # light blue
               Rest="gray15")              # dark grey
cols.woody <- c(Woody="#533d0c",           # brown
                Herbaceous="#799321")      # green
cols.woody <- c(Woody="black",
                Herbaceous="white")
cols.shading <- "#eeb911"                  # yellow
cols.tropical <- c(tropical="#ea7518",     # orange
                   temperate="#62a184")    # teal

## The 'load.woodiness.data.genus' function hides a reasonable amount
## of data cleaning required to load the data.  This mostly involves
## matching the woodiness data set from Zanne et al. to our species
## list (derived from The Plant List), cleaning up synonomies, and
## then collapsing down to genus.
##
## The final object has columns
## * genus, family, order -- taxonomic information
## * W, V, H              -- number of species scored as woody,
##                           variable, herbaceous (respectively)
## * N                    -- number of species in the genus
## * K                    -- number of species with known state,
##                           after dropping all "variable" species
## * p                    -- fraction of known species that are woody,
##                           after dropping all "variable" species
dat.g <- load.woodiness.data.genus()

## Here is the simulator.

## For each genus:

## A. If the genus has a valid fraction of species (i.e. K > 0 so p is
## defined), then sample the number of species that are woody from
## either
##    - the binomial distribution (strong prior; assuming known
##      species were sampled with replacement from the pool of
##      species); or
##    - the hypergeometric distribution (weak prior; assuming that
##      known species were sampled without replacement from the pool
##      of species).

## B. If the genus has no valid fraction of species (i.e., K == 0 so p
## is undefined), then sample from the emperical distribution of
## per-genus fractions.  We're going to feed data into this by
## taxonomic order, so this will come from the per-order distribution.
sim <- function(x, nrep, with.replacement=TRUE, p=1/20) {
  ## First, focus on cases where we have a valid estimate of the
  ## fraction of species that are woody (i.e., at least one known
  ## species).
  ok <- !is.na(x$p)

  w <- matrix(NA, nrow(x), nrep)

  ## A: genera with any known species
  if (with.replacement)
    w[ok,] <- x$W[ok] + rbinom(sum(ok), x$N[ok]-x$K[ok], x$W[ok]/x$K[ok])
  else
    w[ok,] <- t(sapply(which(ok), function(i)
                       rhyper2(nrep, x$H[i], x$W[i], x$N[i])))

  ## B: genera with no known species
  n.unk <- sum(!ok)
  w[!ok,] <- apply(w[ok,,drop=FALSE] / x$N[ok], 2, function(y)
                   rbinom(n.unk, x$N[!ok], quantile(y, runif(n.unk))))
  
  rownames(w) <- x$genus

  summarise.sim(w, x[c("order", "family", "genus", "N")])
}

## This collects up the results at different taxonomic levels.
summarise <- function(x, p=1/20)
  structure(c(mean(x), quantile(x, c(p/2, 1-p/2))),
            names=c("mean", "lower", "upper"))
summarise.sim <- function(w, info) {
  order <- info$order[[1]]

  ## Genus is easy;
  w.g <- cbind(info, t(apply(w, 1, summarise)))

  ## Family is a pain:
  w.f <- do.call(rbind,
                 lapply(split(as.data.frame(w), info$family), colSums))
  w.f <- t(apply(w.f, 1, summarise))
  w.f <- data.frame(order=order, family=rownames(w.f),
                    N=as.integer(tapply(info$N, info$family, sum)),
                    w.f, stringsAsFactors=TRUE)
  rownames(w.f) <- NULL
  
  ## Order is easy; we are guaranteed to have just one order here, so:
  w.o <- data.frame(order=order, N=sum(info$N),
                    t(summarise(colSums(w))), stringsAsFactors=FALSE)

  ret <- list(genus=w.g, family=w.f, order=w.o)
  attr(ret, "total") <- colSums(w)
  ret
}

rhyper2 <- function(nn, s0, s1, xn, fraction=FALSE) {
  x1 <- seq(s1, xn - s0)
  x0 <- xn - x1
  p1 <- dhyper(s1, x1, x0, s0+s1)
  p1 <- p1 / sum(p1)
  x1[sample(length(p1), nn, TRUE, p1)]
}

do.simulation <- function(dat.g, nrep, with.replacement) {
  f <- function(level) {
    ret <- do.call(rbind, lapply(res, "[[", level))
    rownames(ret) <- NULL
    ret[c("p.mean", "p.lower", "p.upper")] <-
      ret[c("mean", "lower", "upper")] / ret[["N"]]
    ret
  }

  res <- lapply(split(dat.g, dat.g$order),
                sim, nrep, with.replacement)
  total <- rowSums(sapply(res, attr, "total"))
  overall <- summarise(total)
  overall.p <- overall / sum(dat.g$N)

  list(genus=f("genus"), family=f("family"), order=f("order"),
       overall=overall, overall.p=overall.p, total=total)
}

res.b <- do.simulation(dat.g, 1000, TRUE)  # binomial - with replacement
res.h <- do.simulation(dat.g, 1000, FALSE) # hypergeometric - without

fig.distribution.raw <- function(res.b, res.h) {
  n.spp <- sum(res.b$order$N)
  p.b <- res.b$total / n.spp * 100
  p.h <- res.h$total / n.spp * 100
  
  r <- range(p.b, p.h)
  br <- seq(r[1], r[2], length.out=30)

  h.b <- hist(p.b, br, plot=FALSE)
  h.h <- hist(p.h, br, plot=FALSE)

  xlim <- c(42, 50)
  ylim <- range(h.b$density, h.h$density)

  cols <- cols.methods
  
  op <- par(mar=c(4.1, 4.1, .5, .5))
  on.exit(par(op))
  plot(h.b, col=cols[1], xlim=xlim, ylim=ylim, freq=FALSE, yaxt="n",
       ylab="",
       xlab="Percentage of woody species among all vascular plants",
       main="")
  box(bty="l")
  lines(h.h, col=cols[2], freq=FALSE)
  mtext("Probability density", 2, line=.5)
}

##+ distribution_raw,fig.cap="Distribution of simulated woodiness percentage"
fig.distribution.raw(res.b, res.h)

## Now, look at the distributions of woodiness among families:
fig.fraction.by.group <- function(res.b, res.h, dat.g, level="genus") {
  op <- par(mfrow=c(2, 1), mar=c(2, 2, .5, .5), oma=c(2, 0, 0, 0))
  on.exit(par(op))
  lwd <- 1.5

  n.br <- c(genus=50, family=40, order=30)[[level]]
  tmp <- aggregate(dat.g[c("W", "K")], dat.g[level], sum)
  tmp <- tmp[tmp$K >= 10,] # at least 10 records per group
  h <- hist(100 * tmp$W / tmp$K, n.br, plot=FALSE)

  plot(NA, xlim=c(0, 100), ylim=range(0, h$density),
       xaxt="n", yaxt="n", bty="l", xlab="", ylab="")
  mtext("Probability density", 2, line=.5)
  axis(1, tick=TRUE, label=FALSE)
  label(.02, .96, 1)
  hist.outline(h, "black", lwd=lwd)  
  
  cols <- cols.methods

  x.b <- res.b[[level]]
  x.h <- res.h[[level]]
  h.b <- hist(100*x.b$p.mean[x.b$N >= 10], n=n.br, plot=FALSE)
  h.h <- hist(100*x.h$p.mean[x.h$N >= 10], n=n.br, plot=FALSE)
  ylim <- range(h.b$density, h.h$density)
  plot(NA, xlim=c(0, 100), ylim=ylim,
       xlab="", ylab="", yaxt="n", bty="n", bty="l")
  mtext("Probability density", 2, line=.5)
  mtext(paste("Percentage of woody species in", level), 1, outer=TRUE,
        line=.5)

  hist.outline(h.b, cols[1], lwd=lwd)
  hist.outline(h.h, cols[2], lwd=lwd)
  
  legend("topleft", c("Replacement (strong prior)",
                      "No replacement (weak prior)"),
         col=cols, lty=1, bty="n", cex=.85, inset=c(.1, 0), lwd=lwd)
  label(.02, .96, 2)
}

##+ fraction_by_genus,fig.cap="Fraction of woodiness by genus"
fig.fraction.by.group(res.b, res.h, dat.g, "genus")
##+ fraction_by_family,fig.cap="Fraction of woodiness by family"
fig.fraction.by.group(res.b, res.h, dat.g, "family")
##+ fraction_by_order,fig.cap="Fraction of woodiness by order"
fig.fraction.by.group(res.b, res.h, dat.g, "order")

## Phylogeny at the level of order
phy.o <- build.order.tree(dat.g)

##+ fraction_phy_binomial,fig.cap="Woodiness percentage by order"
fig.fraction.on.phylogeny(phy.o, res.b)

##+ fraction_phy_hypergeometric,fig.cap="Woodiness percentage by order"
fig.fraction.on.phylogeny(phy.o, res.h)

## # Survey:

## If you want to redo the lat/long calculations, this is the function
## to run:
##   `build.country.list()`
## before running `load.survey()` below.
## This downloads some files from GBIF and extracts country-level
## lat/long values from it to rebuild `geo/country_coords.csv`.  This
## file is already available in the repository, so this will rarely be
## needed.  It is how this file was generated, however.
d.survey <- load.survey()

## Convert estimates to normal using logit transformation:
d.survey$Estimate.logit <- boot::logit(d.survey$Estimate / 100)

## Model with training and familiarity as factors:
res <- lm(Estimate.logit ~ Training + Familiarity, d.survey)
summary(res)
anova(res)

## Regression against |latitude|:
res.lat <- lm(Estimate.logit ~ abs(Lat), d.survey)
anova(res.lat)
summary(res.lat)

## Here is the fitted result:
##+ survey_by_latidude,fig.cap="Fitted latitude survey regression"
plot(Estimate.logit ~ abs(Lat), d.survey)
abline(res.lat)

## As a categorical tropical/non-tropical variable:
res.tro <- lm(Estimate.logit ~ Tropical, d.survey)
anova(res.tro)
summary(res.tro)

fig.survey.results <- function(d.survey, res.b, res.h) {
  ci <- 100*cbind(res.b$overall.p, res.h$overall.p)
  cols <- cols.methods
  cols.tr <- diversitree:::add.alpha(cols, .5)

  op <- par(no.readonly=TRUE)
  on.exit(par(op))

  layout(rbind(1:2), widths=c(4, 5))
  par(mar=c(6.5, 2, .5, .5), oma=c(0, 2, 0, 0))
  plot(Estimate ~ Familiarity, d.survey, col=cols.shading, axes=FALSE,
       xlab="", ylab="", bty="l",
       ylim=c(0, 100))
  axis(2, las=1)
  text(1:4, -5, levels(d.survey$Familiarity),
       srt=-55, xpd=NA, adj=c(0, NA), cex=.85)
  mtext("Estimate of percentage woodiness", 2, line=2.75)
  label(.02, .96, 1)

  usr <- par("usr")
  rect(usr[1], ci["lower",], usr[2], ci["upper",], col=cols.tr,
       border=NA)
  abline(h=ci["mean",], col=cols)

  plot(Estimate ~ Training, d.survey, col=cols.shading, axes=FALSE,
       xlab="", ylab="", bty="l", ylim=c(0, 100))
  axis(2, las=1)
  xl <- c("Postgrad","Part postgrad","Undergrad","Part undergrad", "None")
  text(1:5, -5, xl,
       srt=-55, xpd=TRUE, adj=c(0, NA), cex=.85) 
  label(.02, .96, 2) 

  usr <- par("usr")
  rect(usr[1], ci["lower",], usr[2], ci["upper",], col=cols.tr,
       border=NA)
  abline(h=ci["mean",], col=cols)
}

##+ survey_results,fig.cap="Survey results by predictor"
fig.survey.results(d.survey, res.b, res.h)

fig.survey.distribution <- function(d.survey, res.b, res.h) {
  op <- par(mfrow=c(2, 1), mar=c(2, 4, .5, .5), oma=c(2, 0, 0, 0))
  on.exit(par(op))
  lwd <- 1.5

  ci <- 100*cbind(res.b$overall.p, res.h$overall.p)
  hist(d.survey$Estimate, xlim=c(0, 100), las=1, col=cols.shading,
       xaxt="n", xlab="", ylab="Number of responses", main="")
  box(bty="l")
  axis(1, label=FALSE)
  label(.02, .96, 1)

  usr <- par("usr")
  rect(ci["lower",], usr[3], ci["upper",], usr[4],
       col=diversitree:::add.alpha(cols.methods, .5), border=NA)
  abline(v=ci["mean",], col=cols.methods)

  h.tropical <- hist(d.survey$Estimate[d.survey$Tropical], plot=FALSE)
  h.temperate <- hist(d.survey$Estimate[!d.survey$Tropical], plot=FALSE)

  ylim <- range(h.tropical$counts, h.temperate$counts)
  plot(NA, xlim=c(0, 100), ylim=ylim, las=1, xlab="",
       ylab="Number of responses", bty="n", bty="l")
  mtext("Estimate of percentage woodiness", 1, outer=TRUE, line=.5)

  hist.outline(h.tropical,  cols.tropical[1], lwd=lwd, density=FALSE)
  hist.outline(h.temperate, cols.tropical[2], lwd=lwd, density=FALSE)

  usr <- par("usr")
  rect(ci["lower",], usr[3], ci["upper",], usr[4],
       col=diversitree:::add.alpha(cols.methods, .5), border=NA)
  abline(v=ci["mean",], col=cols.methods)

  label(.02, .96, 2)

  legend("topright", c("Tropical", "Temperate"), lwd=lwd,
         col=cols.tropical, bty="n", cex=.75)
}

##+ survey_distribution,fig.cap="Distribution of survey results"
fig.survey.distribution(d.survey, res.b, res.h)

## Produce versions for publication:
if (!interactive()) {
  to.pdf("doc/figs/fraction-by-genus.pdf", 6, 6,
         fig.fraction.by.group(res.b, res.h, dat.g, "genus"))
  to.pdf("doc/figs/fraction-by-family.pdf", 6, 6,
         fig.fraction.by.group(res.b, res.h, dat.g, "family"))
  to.pdf("doc/figs/fraction-by-order.pdf", 6, 6,
         fig.fraction.by.group(res.b, res.h, dat.g, "order"))

  to.pdf("doc/figs/fraction-on-phylogeny.pdf", 6, 6,
         fig.fraction.on.phylogeny(phy.o, res.b))

  to.pdf("doc/figs/fraction-on-phylogeny-supp.pdf", 6, 6,
         fig.fraction.on.phylogeny(phy.o, res.h))

  to.pdf("doc/figs/distribution-raw.pdf", 6, 4,
         fig.distribution.raw(res.b, res.h))

  to.pdf("doc/figs/survey-results.pdf", 6, 4,
         fig.survey.results(d.survey, res.b, res.h))

  to.pdf("doc/figs/survey-distribution.pdf", 6, 5,
         fig.survey.distribution(d.survey, res.b, res.h))
}
