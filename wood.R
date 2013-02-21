## % How much of the world is woody, first analysis.
library(multicore)
library(diversitree)
library(RColorBrewer)
path.forest <- readLines("~/.forest_path")
source("wood-functions.R")

## Colours used throughout:
cols.methods <- c(strong="#a63813", weak="#4d697f") # red, blue
cols.tree <- c(Monilophytes="#a63813",       # reddish brown
               Gymnosperms="#21313b", # dark brown
               BasalAngiosperms="#eeb911", # yellow (basal)
               Monocots="#204d14",   # green
               Eudicots="#4d697f",   # light blue
               Rest="gray15")
cols.woody <- c(Woody="#533d0c", Herbaceous="#799321")
cols.shading <- "#eeb911"
cols.tropical <- c(tropical="#ea7518", temperate="#62a184")


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
sim <- function(x, nrep, mode="weak") {
  i <- is.na(x$p)
  w <- matrix(NA, nrow(x), nrep)
  n.unk <- sum(i)
  if ( mode == "weak" ) {
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

## Sample from the hypergeometric distribution.  For R's rhyper, we
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

## Caching simulation run:
do.simulation <- function(type, nrep, regenerate=FALSE) {
  filename <- sprintf("sim.%s.rds", type)
  if ( !regenerate && file.exists(filename) )
    return(readRDS(filename))
  dat.g.split <- split(dat.g, dat.g$order)
  set.seed(1) # repeatable seed
  sim <- lapply(dat.g.split, sim, nrep, type)
  saveRDS(sim, filename)
  sim
}

## "Weak prior" sampling with replacement
sim.weak <- do.simulation("weak", 1000)
## "Strong prior" sampling based on observed fractions of species.
sim.strong <- do.simulation("strong", 1000)

## Need to process the simulated data to produce some summary
## statistics.  We care about:
##
## * overall woodiness fraction
## * per family woodiness fraction
## * per order woodiness fraction
process <- function(samples, dat.g) {
  samples.f <- do.call(rbind, samples)
  ## Temporary for now:
  samples.f <- samples.f[!duplicated(rownames(samples.f)),]
  samples.o <- do.call(rbind, lapply(samples, colSums))

  n.f <- tapply(dat.g$N, dat.g$family, sum)
  n.o <- tapply(dat.g$N, dat.g$order, sum)

  ## From this, compute the per-family and per-order fraction
  nrep <- ncol(samples[[1]])
  prop.f <- samples.f / rep(n.f[rownames(samples.f)], nrep)
  prop.o <- samples.o / rep(n.o[rownames(samples.o)], nrep)
  prop.all <- colSums(samples.o) / sum(n.o)

  ## And the mean fraction woody per family and order:
  f <- function(x) {
    ret <- c(mean(x), quantile(x, c(1, 39)/40))
    names(ret) <- c("mean", "lower", "upper")
    ret
  }

  list(overall=f(prop.all),
       family=as.data.frame(t(apply(prop.f, 1, f))),
       order=as.data.frame(t(apply(prop.o, 1, f))),
       distribution=prop.all)
}

res.strong <- process(sim.strong, dat.g)
res.weak   <- process(sim.weak,   dat.g)

## Now, look at the distributions of woodiness among families:
fig.fraction.by.genus <- function(res.strong, res.weak) {
  par(mfrow=c(2, 1), mar=c(2, 2, .5, .5), oma=c(2, 0, 0, 0))
  lwd <- 1.5

  tmp <- dat.g$p[dat.g$K >= 10] # genera with 10 records
  ## TODO: Do we want the distribution of % woody-per-family here too?
  ## Is that something that can actually be easily computed (probably
  ## not).
  h <- hist(100*tmp, 50, plot=FALSE)
  plot(NA, xlim=c(0, 100), ylim=range(0, h$density),
       xaxt="n", yaxt="n", bty="l", xlab="", ylab="")
  mtext("Probability density", 2, line=.5)
  axis(1, tick=TRUE, label=FALSE)
  label(.02, .96, "a)")
  hist.outline(h, "black", lwd=lwd)  
  
  cols <- cols.methods

  h.strong <- hist(100*res.strong$family$mean, n=50, plot=FALSE)
  h.weak   <- hist(100*res.weak$family$mean, n=50, plot=FALSE)
  ylim <- range(h.strong$density, h.weak$density)
  plot(NA, xlim=c(0, 100), ylim=ylim,
       xlab="", ylab="", yaxt="n", bty="n", bty="l")
  mtext("Probability density", 2, line=.5)
  mtext("Percentage of woody species in genus", 1, outer=TRUE,
        line=.5)

  hist.outline(h.strong, cols[1], lwd=lwd)
  hist.outline(h.weak,   cols[2], lwd=lwd)
  
  legend("topleft", c("No replacement (strong prior)",
                      "Replacement (weak prior)"),
         col=cols, lty=1, bty="n", cex=.85, inset=c(.1, 0), lwd=lwd)
  label(.02, .96, "b)")
}

fig.fraction.on.phylogeny <- function(res) {
  phy.o <- build.order.tree(dat.g)
  ## Higher level taxonomy
  hlt <- read.csv("high-level-taxonomy.csv", stringsAsFactors=FALSE)
  phy.group <- hlt$Group[match(phy.o$tip.label, hlt$Order)]
  tmp <- 
    lapply(seq_len(max(phy.o$edge)), function(x)
           if ( x <= length(phy.o$tip.label) ) phy.group[[x]] else
           unique(phy.group[get.descendants(x, phy.o, TRUE)]))
  grp <- sapply(tmp, function(x) if (length(x) == 1) x else "Rest")

  col <- unname(cols.tree[grp])
  col2 <- col[match(phy.o$edge[,2], seq_along(grp))]

  p <- structure(res$order[["mean"]], names=rownames(res$order))
  p <- p[phy.o$tip.label]

  t <- max(branching.times(phy.o))
  offset <- .15

  plt <- 
    diversitree:::plot2.phylo(phy.o, type="fan", cex=.5, no.margin=TRUE,
                              label.offset=t * .15, font=1,
                              edge.col=col2,
                              n.taxa=log1p(phy.o$n.taxa))
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

  polygon(xx1, yy1, border=NA, col=pie.col[2])
  polygon(xx2, yy2, border=NA, col=pie.col[1])

  cex.legend <- 2/3
  str <- str.leg <- setdiff(names(cols.tree), "Rest") # drop backbone
  str.leg[str.leg == "BasalAngiosperms"] <- '"Palaeodiocots"'
  legend("topleft", str.leg, fill=cols.tree[str],
         cex=cex.legend, bty="n", border=NA)
  legend("topright", names(cols.woody), fill=cols.woody,
         cex=cex.legend, bty="n", border=NA)
}

## If you want to redo the lat/long calculations, this is the function
## to run:
##   build.country.list()
## before running load.survey() below.
## This downloads some files from GBIF and extracts country-level
## lat/long values from it to rebuild survey/country_coords.csv.  This
## file is already available in the repository, so this will rarely be
## needed.  It is how this file was generated, however.
d.survey <- load.survey()

## Convert estimates to normal using logit transformation
d.survey$Estimate.logit <- boot::logit(d.survey$Estimate / 100)
res <- lm(Estimate.logit ~ Training + Familiarity, d.survey)
summary(res)
anova(res)

res.lat <- lm(Estimate.logit ~ abs(Lat), d.survey)
anova(res.lat)
summary(res.lat)

res.tro <- lm(Estimate.logit ~ Tropical, d.survey)

## Here is the fitted result:
## plot(Estimate.logit ~ abs(Lat), d.survey)
## abline(res.lat)

fig.survey.results <- function(d.survey, res.strong, res.weak) {
  ci <- 100*cbind(res.strong$overall, res.weak$overall)
  cols <- cols.methods
  cols.tr <- diversitree:::add.alpha(cols, .5)
  
  layout(rbind(1:2), widths=c(4, 5))
  par(mar=c(6.5, 2, .5, .5), oma=c(0, 2, 0, 0))
  plot(Estimate ~ Familiarity, d.survey, col=cols.shading, axes=FALSE,
       xlab="", ylab="", bty="l",
       ylim=c(0, 100))
  axis(2, las=1)
  text(1:4, -5, levels(d.survey$Familiarity),
       srt=-55, xpd=NA, adj=c(0, NA), cex=.85)
  mtext("Estimate of percentage woodiness", 2, line=2.75)
  label(.02, .96, "a)")

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
  label(.02, .96, "b)") 

  usr <- par("usr")
  rect(usr[1], ci["lower",], usr[2], ci["upper",], col=cols.tr,
       border=NA)
  abline(h=ci["mean",], col=cols)
}

fig.distribution.raw <- function(res.strong, res.weak) {
  r <- range(res.strong$distribution, res.weak$distribution)
  br <- seq(r[1], r[2], length.out=30)*100

  h.strong <- hist(100*res.strong$distribution, br, plot=FALSE)
  h.weak   <- hist(100*res.weak$distribution,   br, plot=FALSE)

  xlim <- c(42, 50)
  ylim <- range(h.strong$density, h.weak$density)

  cols <- cols.methods
  
  par(mar=c(4.1, 4.1, .5, .5))
  plot(h.strong, col=cols[1], xlim=xlim, ylim=ylim, freq=FALSE, yaxt="n",
       ylab="",
       xlab="Percentage of woody species among all vascular plants",
       main="")
  box(bty="l")
  lines(h.weak, col=cols[2], freq=FALSE)
  mtext("Probability density", 2, line=.5)
}

fig.survey.distribution <- function(d.survey, res.strong, res.weak) {
  par(mfrow=c(2, 1), mar=c(2, 4, .5, .5), oma=c(2, 0, 0, 0))
  lwd <- 1.5

  ci <- 100*cbind(res.strong$overall, res.weak$overall)
  hist(d.survey$Estimate, xlim=c(0, 100), las=1, col=cols.shading,
       xaxt="n", xlab="", ylab="Number of responses", main="")
  box(bty="l")
  axis(1, label=FALSE)
  label(.02, .96, "a)")

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

  label(.02, .96, "b)")

  legend("topright", c("Tropical", "Temperate"), lwd=lwd,
         col=cols.tropical, bty="n", cex=.75)
}


to.pdf("doc/figs/fraction-by-genus.pdf", 6, 6,
       fig.fraction.by.genus(res.strong, res.weak))

to.pdf("doc/figs/fraction-on-phylogeny.pdf", 6, 6,
       fig.fraction.on.phylogeny(res.strong))

to.pdf("doc/figs/fraction-on-phylogeny-supp.pdf", 6, 6,
       fig.fraction.on.phylogeny(res.weak))

to.pdf("doc/figs/distribution-raw.pdf", 6, 4,
       fig.distribution.raw(res.strong, res.weak))

to.pdf("doc/figs/survey-results.pdf", 6, 4,
       fig.survey.results(d.survey, res.strong, res.weak))

to.pdf("doc/figs/survey-distribution.pdf", 6, 5,
       fig.survey.distribution(d.survey, res.strong, res.weak))



