#!/usr/bin/env Rscript
path <- "data/zae"

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

paths <- lapply(urls, download.maybe, path=path)

# Then unzip the tree that we use:
tmp <- tempdir()
unzip(file.path(path, 'PhylogeneticResources.zip'),
      'PhylogeneticResources/Vascular_Plants_rooted.dated.tre',
      junkpaths=TRUE, exdir=tmp)
invisible(file.rename(file.path(tmp, "Vascular_Plants_rooted.dated.tre"),
                      file.path(path, "Vascular_Plants_rooted.dated.tre")))

# Load the lookup table -- we need to expand this a bit:
lookup <- read.csv(file.path(path, "Spermatophyta_Genera.csv"),
                   stringsAsFactors=FALSE)
names(lookup)[1] <- "genus"
lookup <- lookup[c("genus", "family", "order")]

# Two genera need families set:
lookup$family[lookup$genus == "Peltanthera"] <- "Gesneriaceae"
lookup$family[lookup$genus == "Brachynema" ] <- "Olacaceae"

write.csv(lookup, file.path(path, "genus_order_lookup.csv"),
          row.names=FALSE)
