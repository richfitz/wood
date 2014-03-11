# How Much Wood?

Avoid hammering TPL's servers:

```
make downloaded-data-bulk-fetch
make downloaded-data-unpack
```

Then simply

```
make
```

## Requirements

To generate the report, we depend on the non-CRAN package [sowsear](https://github.com/richfitz/sowsear).  The easiest way to install that is with [devtools](https://github.com/hadley/devtools)

```
library(devtools)
install_github("richfitz/sowsear")
```

(install `devtools` with `install.packages("devtools")` if you don't already have it).

The other dependencies are for the packages `ape` and `diversitree`, both of which are on CRAN.  To recreate the geographic data (in `data/geo/country_coords.csv`) the `rgdal` package is required.
