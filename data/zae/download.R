#!/usr/bin/env Rscript

# These are the files to download
prefix <- 'http://datadryad.org/bitstream/handle/10255/'
suffix <- c('dryad.55304/Spermatophyta_Genera.csv?sequence=2',
            'dryad.55447/GlobalWoodinessDatabase.csv?sequence=1',
            'dryad.55548/PhylogeneticResources.zip?sequence=1')
urls <- paste0(prefix, suffix)

download.maybe <- function(url, refetch=FALSE, path=".") {
  dest <- file.path(path, sub("\\?.+$", "", basename(url)))
  if (refetch || !file.exists(dest))
    download.file(url, dest)
  dest
}

paths <- lapply(urls, download.maybe)

# Then unzip the tree that we use:
if (!file.exists("Vascular_Plants_rooted.dated.tre"))
  unzip('PhylogeneticResources.zip',
        'PhylogeneticResources/Vascular_Plants_rooted.dated.tre',
        junkpaths=TRUE)

# Load the lookup table -- we need to expand this a bit:
if (!file.exists("genus_order_lookup.csv")) {
  lookup <- read.csv("Spermatophyta_Genera.csv", stringsAsFactors=FALSE)
  names(lookup)[1] <- "genus"
  lookup <- lookup[c("genus", "family", "order")]

  # Two genera need families set:
  lookup$family[lookup$genus == "Peltanthera"] <- "Gesneriaceae"
  lookup$family[lookup$genus == "Brachynema" ] <- "Olacaceae"

  write.csv(lookup, "genus_order_lookup.csv", row.names=FALSE)
}
