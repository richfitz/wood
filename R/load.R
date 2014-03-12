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
