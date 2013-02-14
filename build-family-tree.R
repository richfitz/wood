#!/usr/bin/env Rscript

## Temporary script to build the family tree.

library(diversitree)
path.forest <- readLines("~/.forest_path")

get.class <- function(phy.f) {
  lab <- phy.f$node.label
  keep <- grepl("ales$", lab)
  keep[lab == "Monilophyte"] <- TRUE
  keep[lab %in% c("Gleicheniales", "Schizaeales", "Salviniales")] <- FALSE
  keep[grep("_To_", lab)] <- FALSE
  lab[!keep] <- NA

  descendants <- diversitree:::descendants
  desc.spp <- function(node, phy) {
    n.spp <- length(phy$tip.label)
    ret <- descendants(match(node, phy$node.label) + n.spp, phy$edge)
    phy$tip.label[ret[ret <= n.spp]]
  }
  orders <- lab[!is.na(lab)]
  tmp <- lapply(orders, desc.spp, phy.f)
  class <- data.frame(order=rep(orders, sapply(tmp, length)),
                      family=unlist(tmp), stringsAsFactors=FALSE)
  any(duplicated(unlist(class$family)))

  ## There are a handful of otherws potentially worth keeping, too.
  grp <- class$order[match(phy.f$tip.label, class$family)]
  grp[phy.f$tip.label == "Arecaceae"] <- "Arecacales"
  grp[phy.f$tip.label == "Boraginaceae"] <- "Boraginaceae"
  grp[phy.f$tip.label == "Vitaceae"] <- "Vitales"
  grp[phy.f$tip.label == "Dilleniaceae"] <- "Dilleniaceae"
  grp[phy.f$tip.label == "Escalloniaceae"] <- "Escalloniaceae"
  grp[phy.f$tip.label == "Sabiaceae"] <- "Sabiaceae"
  grp[phy.f$tip.label == "Chloranthaceae"] <- "Chloranthaceae"
  grp[phy.f$tip.label == "Buxaceae"] <- "Buxaceae"
  grp[phy.f$tip.label == "Icacinaceae"] <- "Icacinaceae"
  grp[phy.f$tip.label == "Acoraceae"] <- "Acoraceae"
  grp[phy.f$tip.label == "Paracryphiaceae"] <- "Paracryphiaceae"
  grp
}

build.family.tree <- function(regenerate=FALSE) {
  filename <- "phy.f.rds"
  if ( !regenerate && file.exists(filename) ) {
    phy.f <- readRDS(filename)
  } else {
    phy <- read.tree(file.path(path.forest,
                               "taxonomic/trees/spLevelApgBackbone.tre"))
    i <- phy$node.label == ""
    phy$node.label[i] <- sprintf("node.%d", (1:phy$Nnode)[i])
    problems <-
      readLines(file.path(path.forest,
                          "taxonomic/trees/problemSpeciesAPG.txt"))
    problems <- sub(" *$", "", problems) # kill trailing whitespace
    phy <- diversitree:::drop.tip.fixed(phy, problems)
    phy <- compute.brlen(phy, method="Grafen", power=.45)
    families <-
      read.csv(file.path(path.forest, "taxonomic/genus_order_lookup.csv"),
               stringsAsFactors=FALSE)
    genus <- sub("_.+$", "", phy$tip.label)
    family <- families$family[match(genus, families$genus)]
    family[is.na(family)] <- genus[is.na(family)] # Fabaceae x 5
    phy.f <- clades.from.classification(phy, family, check=FALSE)
    ord <- get.class(phy.f)
    names(ord) <- phy.f$tip.label
    phy.f$class <- ord
    saveRDS(phy.f, filename)
  }
  phy.f
}

if ( !interactive() ) {
  build.family.tree()
}
