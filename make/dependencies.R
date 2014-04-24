#!/usr/bin/env Rscript
pkgs <- as.data.frame(read.dcf(".packrat/packrat.lock"),
                      stringsAsFactors=FALSE)
pkgs <- pkgs[!is.na(pkgs$Package),]
pkgs.installed <- installed.packages()

# First, find out which packages are missing:
pkgs.missing <- setdiff(pkgs$Package, rownames(pkgs.installed))

if (length(pkgs.missing) > 0) {
  # Missing and on cran -- just install these:
  pkgs.missing.info <- pkgs[match(pkgs.missing, pkgs$Package),]
  pkgs.missing.cran <-
    pkgs.missing.info$Package[pkgs.missing.info$Source == "CRAN"]
  if (length(pkgs.missing.cran) > 0) {
    install.packages(pkgs.missing.cran)
  }

  # Then get the github ones:
  pkgs.missing.info.github <-
    pkgs.missing.info[pkgs.missing.info$Source == "github",]
  pkgs.missing.github <-
    paste(pkgs.missing.info.github$GithubUsername,
          pkgs.missing.info.github$GithubRepo, sep="/")
  for (p in pkgs.missing.github) {
    devtools::install_github(p)
  }
}

# Now check which packages are out of date.  By this point, everything
# is installed so we're cool.
pkgs.installed <- installed.packages()

pkgs.missing <- setdiff(pkgs$Package, rownames(pkgs.installed))
if (length(pkgs.missing) > 0) {
  stop("Still have missing packages somehow: ",
       paste(pkgs.missing, collapse=", "))
}

target.version <- numeric_version(pkgs$Version)

i <- match(pkgs$Package, rownames(pkgs.installed))
installed.version <-
  numeric_version(numeric_version(pkgs.installed[i, "Version"]))

old <- installed.version < target.version
if (any(old)) {
  str <- sprintf("\t%s: installed %s, expected %s", pkgs$Package[old],
                 target.version[old], installed.version[old])
  message("Some packages are old and should probably be upgraded:")
  message(paste(str, collapse="\n"))
  message("Commands to upgrade:")

  pkgs.old <- pkgs[old,,drop=FALSE]
  old.cran <- pkgs.old$Package[pkgs.old$Source == "CRAN"]

  old.github <-
    paste(pkgs.old$GithubUsername[pkgs.old$Source == "github"],
          pkgs.old$GithubRepo[pkgs.old$Source == "github"],
          sep="/")

  in.quotes <- function(x) sprintf('"%s"', x)

  str <- character(0)
  if (length(old.cran) > 0) {
    str <- c(str, sprintf("install.packages(c(%s))",
                          paste(in.quotes(old.cran), collapse=", ")))
  }
  if (length(old.github) > 0) {
    str <- c(str, sprintf("devtools::install_github(%s)",
                          in.quotes(old.github)))
  }

  message(paste0("\t", str, collapse="\n"))

  message("NOTE: things may well work fine with the installed versions")
}
