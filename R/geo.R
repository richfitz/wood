## This downloads a bunch of data from gbif and gets latitude and
## longitude for all countries.  This is saved in the
## data/geo/country_coords.csv file, so does not need to be run except
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
