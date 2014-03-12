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

We require a few packages, namely `dplyr`, `diversitree`, `RCurl` and `knitr`.  These can all be installed off CRAN.

To generate the report, we depend on the non-CRAN package [sowsear](https://github.com/richfitz/sowsear).  The easiest way to install that is with [devtools](https://github.com/hadley/devtools)

```
library(devtools)
install_github("richfitz/sowsear")
```

(install `devtools` with `install.packages("devtools")` if you don't already have it).

To recreate the geographic data (in `data/geo/country_coords.csv`) the
`rgdal` package is required.
