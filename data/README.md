This directory contains a few files that are special to this project:

* `survey_results.csv`: survey responses containing a guess as to what
  fraction of species are woody, along with information about botany
  training, familiarity and location.
* `high-level-taxonomy.csv`: Lookup table, manually constructed, that
  links plant orders to high level clades.  This is used only for
  colouring the tree figures.
* `genus_order_lookup_extra.csv`: Lookup table with extra
  genus/family/order relationships, for non-seed plants.  This is
  merged with a similar lookup file from Zanne et al. (see below)
* `synonyms.csv`: Synonyms for species in the woodiness database.
* `geo/country_coords.csv`: nominal latitude and longitude for 259
  countries, for use with the survey results.  This can be regenerated
  by deleting the file and running `make data/geo/country_coords.csv`

Running `make data-raw` (optionally running `make theplantlist-cache-unpack` first to avoid hammering the plant list) will download and create create some additional files:

# directory `zae`

Files from Zanne et al. (2015) [doi:10.1038/nature12872](http://doi.org/10.1038/nature12872).

* `zae/GlobalWoodinessDatabase.csv`: list of species with their woodiness
  status.  See `woodinessReadme.txt` for more information.
* `zae/Spermatophyta_Genera.csv`: Taxonomic information relating plant
  genera, families and orders.  This is used for computing the numbers
  of woody species at the family and order level.
* `zae/genus_order_lookup.csv`: Generated from `Spermatophyta_Genera.csv`,
  with just the `genus`, `family` and `order` columns included.
* `zae/Vascular_Plants_rooted.dated.tre`: very large phylogeny
  of plants.  This is used only for display purposes.

# directory `theplantlist`

* `theplantlist/names_accepted.csv`: List of accepted names (and
  families) from [The Plant List](http://theplantlist.org)
* `theplantlist/families`: Lists of families contained within The
  Plant List (used to create `theplantlist/names_accepted.csv`)
* `theplantlist/acceptedNames1.1`: Raw data downloaded from the Plant List.
