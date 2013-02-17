## % How much of the world is woody, first analysis.
library(multicore)
library(diversitree)
library(RColorBrewer)
path.forest <- readLines("~/.forest_path")
source("wood-functions.R")

## Data at the level of genus: has taxonomic information and counts of
## number of species that are woody, herbaceous, known, and known to
## exist.
dat.g <- load.clean.data(TRUE)

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

## "Weak prior" sampling with replacement
sim.weak <- lapply(dat.g.split, sim, nrep)
## "Strong prior" sampling based on the number of species
sim.strong <- lapply(dat.g.split, sim, nrep, "no-replacement")

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
       order=as.data.frame(t(apply(prop.o, 1, f))))
}

res.strong <- process(sim.strong, dat.g)
res.weak   <- process(sim.weak,   dat.g)

## Now, look at the distributions of woodiness among families:
fig.fraction.by.genus <- function(res.strong, res.weak) {
  par(mfrow=c(2, 1), mar=c(2, 1, .5, .5), oma=c(2, 0, 0, 0))

  tmp <- dat.g$p[dat.g$K >= 10] # genera with 10 records
  h <- hist(tmp, 50, main="", xaxt="n", yaxt="n", lab="")
  mtext("Probability density", 2)
  axis(1, tick=TRUE, label=FALSE)
  label(.02, .96, "a)")
  
  cols <- c("#ff0000ff", "#0000ff66")
  hist(res.strong$family$mean, col=cols[1], n=50, freq=FALSE,
       main="", yaxt="n", xlab="", ylab="")
  mtext("Probability density", 2)
  mtext("% woody species in genus", 1, outer=TRUE, line=.5)
  hist(res.weak$family$mean, col=cols[2], n=50, freq=FALSE, add=TRUE)
  legend("topleft", c("No replacement (strong prior)",
                      "Replacement (weak prior)"),
         fill=cols, bty="n", cex=.75, inset=c(.1, 0))
  label(.02, .96, "b)")
}

fig.fraction.on.phylogeny <- function(res) {
  phy.f <- build.family.tree()

  ## This tree has order information:
  phy.f.ord <- phy.f$class

  ## Drop tips with no woodiness estimates (due to taxonomic
  ## differences)
  to.drop <- setdiff(phy.f$tip.label, rownames(res$family))
  phy.f <- diversitree:::drop.tip.fixed(phy.f, to.drop)

  ## And sort the order information by the current taxon labels:
  phy.f.ord <- unname(phy.f.ord[phy.f$tip.label])

  cols <- make.col.function(brewer.pal(9, "Blues")[-1])
  plt <- trait.plot.cont(phy.f, res$family["mean"], list(cols),
                         class=phy.f.ord, w=1/30, font=1, margin=1/3)
  mrca.tipset <- diversitree:::mrca.tipset
  nd <- sapply(sort(unique(phy.f.ord)), function(x)
               mrca.tipset(phy.f, phy.f$tip.label[phy.f.ord==x &
                                                  !is.na(phy.f.ord)]))
  nd.int <- nd[nd > length(phy.f$tip.label)]
  points(plt$xy$xx[nd.int], plt$xy$yy[nd.int], pch=19, cex=2,
         col=cols(res$order[names(nd.int),]$mean))
}

d.survey <- load.survey()

## Convert estimates to normal using logit transformation
d.survey$Estimate.logit <- boot::logit(d.survey$Estimate / 100)
res <- lm(Estimate.logit ~ Training + Familiarity, d.survey)
summary(res)
anova(res)

fig.survey.results <- function(d.survey, res.strong, res.weak) {
  ci <- 100*cbind(res.strong$overall, res.weak$overall)

  layout(rbind(1:2), widths=c(4, 5))
  par(mar=c(6.5, 2, .5, .5), oma=c(0, 2, 0, 0))
  plot(Estimate ~ Familiarity, d.survey, col="lightgrey", axes=FALSE,
       xlab="", ylab="", bty="l",
       ylim=c(0, 100))
  axis(2, las=1)
  text(1:4, -5, levels(d.survey$Familiarity),
       srt=-55, xpd=NA, adj=c(0, NA), cex=.85)
  mtext("Estimate of proportion woodiness", 2, line=2.75)
  label(.02, .96, "a)")

  usr <- par("usr")
  rect(usr[1], ci["lower",], usr[2], ci["upper",], col="#77777777",
       border=NA)
  abline(h=ci["mean",], lty=c(1, 2))

  plot(Estimate ~ Training, d.survey, col="lightgrey", axes=FALSE,
       xlab="", ylab="", bty="l", ylim=c(0, 100))
  axis(2, las=1)
  xl <- c("Postgrad","Part postgrad","Undergrad","Part undergrad", "None")
  text(1:5, -5, xl,
       srt=-55, xpd=TRUE, adj=c(0, NA), cex=.85) 
  label(.02, .96, "b)") 

  usr <- par("usr")
  rect(usr[1], ci["lower",], usr[2], ci["upper",], col="#77777777",
       border=NA)
  abline(h=ci["mean",], lty=c(1, 2)) 
}

to.pdf("doc/figs/fraction-by-genus.pdf", 6, 6,
       fig.fraction.by.genus(res.strong, res.weak))

to.pdf("doc/figs/fraction-on-phylogeny.pdf", 6, 6,
       fig.fraction.on.phylogeny(res.strong))

to.pdf("doc/figs/fraction-on-phylogeny-supp.pdf", 6, 6,
       fig.fraction.on.phylogeny(res.weak))

to.pdf("doc/figs/survey-results.pdf", 6, 4,
       fig.survey.results(d.survey, res.strong, res.weak))

