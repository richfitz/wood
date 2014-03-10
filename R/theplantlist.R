library(RCurl, quietly=TRUE)

plant.list.groups <- function()
  c("angiosperm", "gymnosperm", "pteridophyte")

plant.list.url <- function(family, group) {
  group <- match.arg(tolower(group), plant.list.groups())
  sprintf("http://www.theplantlist.org/1.1/browse/%s/%s/%s.csv",
          toupper(substr(group, 1, 1)), family, family)
}

plant.list <- function(family, group) {
  getURLContent(plant.list.url(family, group))  
}

plant.list.csv <- function(family, path="data/theplantlist") {
  file.path(path, "acceptedNames1.1", paste0(family, ".csv"))
}

plant.list.get <- function(family, group,
                           regenerate=FALSE, verbose=TRUE) {
  file.out <- plant.list.csv(family)
  if (!regenerate && file.exists(file.out)) {
    if (verbose)
      message(sprintf("Skipping %s (%s) -- already exists",
                      family, group))
    invisible(FALSE)
  } else {
    message(sprintf("Fetching %s (%s)", family, group))
    dat <- plant.list(family, group)
    writeLines(dat, file.out)
  }
}

plant.list.get.group <- function(group, ..., path="data/theplantlist") {
  group <- match.arg(tolower(group), plant.list.groups())
  families <- readLines(file.path(path, sprintf("families/%s.txt", group)))
  for (f in families) {
    plant.list.get(f, group, ...)
  }
}

load.names <- function(filename) {
  keep <- c(genus="Genus", species="Species", family="Family",
            major.clade="Major.group", status="Taxonomic.status.in.TPL")
  dat <- read.csv(filename, stringsAsFactors=FALSE)
  dat <- dat[dat$Genus != "", keep]
  names(dat) <- names(keep)
  
  dat$species <- gsub('"', "", dat$species, fixed=TRUE)
  dat
}
