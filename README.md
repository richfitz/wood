# How much of the world is woody?

[Richard G. FitzJohn](http://www.zoology.ubc.ca/~fitzjohn),
[Matthew W. Pennell](http://mwpennell.github.io),
[Amy E. Zanne](http://phylodiversity.net/azanne/),
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
make theplantlist-cache-unpack
```

This route allows you to delete all the data (`make purge`) and easily rerun the analysis (`make theplantlist-cache-unpack all`) without redownloading the data.

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

At present, we depend on the github version of diversitree; install that with 

```
library(devtools)
install_github("richfitz/diversitree")
```

To recreate the geographic data (in `data/geo/country_coords.csv`) the
`rgdal` package is required.

## Using a known set of working packages with `packrat`

Version rot means that while the analysis works now, it may not work in a few years when packages have been updated and changed their APIs.  To guard against this, we have archived a set of known working packages using [packrat](https://github.com/rstudio/packrat).

We didn't want to use packrat *all the time* (our package use is hopefully straightforward enough that a plain installation should work) and we didn't want to bog down the repository with about 20MB of package sources (especially as there are stable canonical sources for almost all packages because CRAN retains sources indefinitely).  As such there is a fairly unfortunate, and likey fragile, bootstrapping procedure for enabling packrat that we have bodged together.

Run

```
make packrat-enable
```

which will download the known set of working sources from our [releases page](https://github.com/richfitz/wood/releases) and copy files over from the `.packrat` directory.  This puts packrat into the state that packrat assumes the project is always in.  Packrat then goes through and compiles all the packages and installs them locally into a directory `library`.  This process can take a while!

To disable packrat (putting the project back to using system-installed packages) run

```
make packrat-disable
```

To update the set of known working packages you can use the normal packrat tools and then run

```
make packrat-update
```

which copies local changes into the `.packrat` directory.  These files can then be committed, though the remote `tar.gz` file would then need updating to share these changes.

To record a set of system-installed packages as working, run

```
make packrat-refresh
```

(note that this also sets the project up to use packrat, so running `make packrat-disable` afterwards is probably wise).

See `make/packrat.mk` for more information on our approach here, which does not fit neatly within packrat's scope.  It's possible that by the time using the archived packages is necessary, better systems for doing this will exist.
