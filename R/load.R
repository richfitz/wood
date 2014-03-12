load.genus.order.lookup <- function() {
  read.csv("output/genus_order_lookup.csv", stringsAsFactors=FALSE)
}

load.woodiness <- function() {
  readRDS("output/woodiness.rds")
}

load.woodiness.genus <- function(extreme=FALSE) {
  readRDS(filename.woodiness.genus(extreme))
}

load.theplantlist <- function() {
  tpl <- read.csv("data/theplantlist/names_accepted.csv",
                  stringsAsFactors=FALSE)
  
  ## For ease later, we prefer Asteraceae/Fabaceae to
  ## Compositae/Leguminosae and fix a spelling mistake for
  ## Dryopteridaceae:
  tr <- c("Compositae"     = "Asteraceae",
          "Leguminosae"    = "Fabaceae",
          "Dryopteridacae" = "Dryopteridaceae")
  i <- match(tpl$family, names(tr))
  j <- !is.na(i)
  tpl$family[j] <- unname(tr[i[j]])
  tpl[names(tpl) != "status"]
}

filename.woodiness.genus <- function(extreme=FALSE) {
  if (identical(extreme, FALSE))
    "output/dat.g.rds"
  else if (extreme %in% c("woody", "herbaceous"))
    sprintf("output/dat.g.%s.rds", substr(extreme, 1, 1))
  else
    stop("Invalid argument for 'extreme'")
}

load.survey <- function() {
  d <- read.csv(file="data/survey_results.csv",
                stringsAsFactors=FALSE)
  names(d) <- c("Time", "Estimate", "Familiarity", "Training",
                "Continent", "Country")

  ## Drop the Time and Continent columns:
  d <- d[!(names(d) %in% c("Time", "Continent"))]

  ## Here are the different familiarity and training categories from
  ## "best" to "worst".
  lvl.familiarity <- c("Very Familiar", "Familiar", "Somewhat Familiar",
                       "What's a Plant?")
  d$Familiarity <- factor(d$Familiarity, lvl.familiarity, ordered=TRUE)
  
  lvl.training <-
    c("Postgraduate degree in botany or a related field",
      "Partially complete postgraduate degree in botany or a related field",
      "Undergraduate degree in botany or a related field",
      "Some botany courses at either an undergraduate or postgraduate level",
      "No formal training in botany")             
  d$Training <- factor(d$Training, lvl.training, ordered=TRUE)

  ## Standardise the country names:
  countries <- read.csv("data/geo/country_coords.csv",
                        stringsAsFactors=FALSE)
  d$Country <- cleanup.country.names(d$Country)
  
  idx <- match(d$Country, countries$Country)
  mssg <- na.omit(d$Country[is.na(idx)])
  if (length(mssg) > 0)
    warning("Dropped countries %s", paste(mssg, collapse=", "))
  d <- cbind(d, countries[idx,c("Long", "Lat")])
  d$Tropical <- abs(d$Lat) < 23 + 26/60

  rownames(d) <- NULL
  d
}

cleanup.country.names <- function(x) {
  ## In cases where multiple countries are given, take the first one:
  x <- sub("( and |, | / | & ).+", "", x)
  ## Trim trailing non-alphabetic characters
  x <- sub("[^A-Za-z]+$", "", x)
  ## Translate inconsistent names:
  translate <- list(France="france",
                    "United States"=c("US", "USA"),
                    "United Kingdom"=c("UK", "Scotland"),
                    "Brazil"="Brasil")
  tr <- cbind(to=rep(names(translate), sapply(translate, length)),
              from=unlist(translate))
  i <- match(x, tr[,"from"])
  x[!is.na(i)] <- unname(tr[i[!is.na(i)],"to"])
  x[x == ""] <- NA
  x
}
