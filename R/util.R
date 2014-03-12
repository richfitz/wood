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
