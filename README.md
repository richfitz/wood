# How Much Wood?

[![Build Status](https://travis-ci.org/richfitz/wood.png?branch=master)](https://travis-ci.org/richfitz/wood)

Avoid hammering TPL's servers:

```
make downloaded-data-bulk-fetch
make downloaded-data-unpack
```

Then simply

```
make
```

This will build processed versions of the data in the `output` directory.  It then converts the file `wood.R` to a [knitr](http://yihui.name/knitr/) script (`wood.Rnw`) and runs `knitr` on this to generate `wood.md` (in [markdown](http://daringfireball.net/projects/markdown/)) and the figures for the paper (in `doc/figs`).

The `wood.md` file is turned into a little html report of the analysis (`wood.html`).

The actual manuscript is in `doc/wood-ms.pdf`.  Compiling this requires LaTeX to be installed.

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
