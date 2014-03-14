# How much of the world is woody?

[Richard G. FitzJohn](http://www.zoology.ubc.ca/~fitzjohn),
[Matthew W. Pennell](mwpennell.wordpress.com),
[Amy E. Zanne](mwpennell.wordpress.com),
[Peter F. Stevens](http://www.missouribotanicalgarden.org/plant-scence/research-staff-article/487/stevens-p-f.aspx),
[David C. Tank](http://www.phylodiversity.net/dtank/), 
[William K. Cornwell](http://www.phylodiversity.net/wcornwell/)

This repository contains all the code and data used in the manuscript.

[![Build Status](https://travis-ci.org/richfitz/wood.png?branch=master)](https://travis-ci.org/richfitz/wood)

## Fetching the data

There are two ways of fetching the required data (see `data/README.md` for information on the data that we depend on).

**Directly**

This downloads data from [Dryad](http://datadryad.org) and from [The Plant List](http://www.theplantlist.org)

```
make data-raw
```

**Avoid hammering TPL**

This fetches a set of data that I've archived.

```
make downloaded-data-unpack
```

This route allows you to delete all the data (`make clean downloaded-data-delete`) and easily rerun the analysis (`make downloaded-data-unpack`).

## Running the analysis

To run the analysis, run the command

```
make
```

This will build processed versions of the data in the `output` directory.  It then converts the file `wood.R` to a [knitr](http://yihui.name/knitr/) script (`wood.Rnw`) and runs `knitr` on this to generate `wood.md` (in [markdown](http://daringfireball.net/projects/markdown/)) and the figures for the paper (in `doc/figs`).

The `wood.md` file is turned into a little html report of the analysis (`wood.html`).

The actual manuscript is in `doc/wood-ms.pdf`.  Compiling this requires LaTeX to be installed.

## Manually running everything

If you don't have `make` installed, then you can compile everything by running

```
source("make/manual.R")
```

(this needs to be run from within R, with the working directory set to the same as this file.  If you use [Rstudio](http://rstudio.com), then opening the file `wood.Rproj` sets the working directory for you.

This will not compile the manuscript `doc/wood-ms.tex` to pdf; if you have LaTeX installed you will need to do that in whatever way you normally would on your system.  However, all figures in the manuscript will be created in `doc/figs`.

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
