## Sample from the hyper geometric distribution.  For R's rhyper, we
## have to provide the number of white and black balls from the urn,
## and it samples them.  However, we don't know what this is.

## What this does is samples the true number of white balls 'nn' times
## and returns that.
rhyper2 <- function(nn, s0, s1, xn, fraction=FALSE) {
  x1 <- seq(s1, xn - s0)
  x0 <- xn - x1
  p1 <- dhyper(s1, x1, x0, s0+s1)
  p1 <- p1 / sum(p1)
  x1[sample(length(p1), nn, TRUE, p1)]
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
