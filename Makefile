RSCRIPT_PKGS := $(shell Rscript -e 'library(methods);writeLines(Sys.getenv("R_DEFAULT_PACKAGES"))')
RSCRIPT = Rscript --default-packages="${RSCRIPT_PKGS},methods"

DATA_RAW = data/zae/genus_order_lookup.csv data/theplantlist/names_accepted.csv
DATA_PROCESSED = output/woodiness.rds output/dat.g.rds \
	output/dat.g.w.rds output/dat.g.h.rds output/phy.o.rds
REPORT = wood.Rmd wood.md wood.html

all: wood.html doc/wood-ms.pdf

data-raw: ${DATA_RAW}
data-processed: ${DATA_PROCESSED}

wood.Rmd: wood.R
	${RSCRIPT} -e "library(sowsear); sowsear('$<', 'Rmd')"
wood.md: wood.Rmd ${DATA_PROCESSED}
	${RSCRIPT} -e "library(knitr); knit('$<')"
wood.html: wood.md
	${RSCRIPT} -e "library(markdown);\
	 opts <- setdiff(markdownHTMLOptions(TRUE), 'base64_images');\
	 markdownToHTML('$<', '$@', options=opts, stylesheet='stylesheet.css', template='.template.html')"

doc/wood-ms.pdf: wood.md
	make -C doc

# This target will never run because it depends on nothing.  But if
# data/geo/country_coords.csv is deleted then this will regenerate
# that file.  It's here for reference only, really.
data/geo/country_coords.csv:
	${RSCRIPT} make/data-geo-country_coords.csv.R

data/zae/genus_order_lookup.csv: make/data-zae.R
	${RSCRIPT} $<

data/theplantlist/names_accepted.csv: make/data-theplantlist.R
	${RSCRIPT} $<

output/genus_order_lookup.csv: make/output-genus_order_lookup.csv.R ${DATA_RAW}
	${RSCRIPT} $<

output/woodiness.rds: make/output-woodiness.rds.R R/load.R R/build.R ${DATA_RAW}
	${RSCRIPT} $<

DATA_GENUS_DEPS = output/genus_order_lookup.csv output/woodiness.rds \
	R/load.R R/build.R

output/dat.g.rds: make/output-dat.g.rds.R ${DATA_GENUS_DEPS}
	${RSCRIPT} $<
output/dat.g.w.rds: make/output-dat.g.w.rds.R ${DATA_GENUS_DEPS}
	${RSCRIPT} $<
output/dat.g.h.rds: make/output-dat.g.h.rds.R ${DATA_GENUS_DEPS}
	${RSCRIPT} $<

output/phy.o.rds: make/output-phy.o.rds.R output/dat.g.rds
	${RSCRIPT} $<

# Cache downloads
RELEASE=v1.0

THEPLANTLIST_CONTENTS = data/theplantlist/acceptedNames1.1
THEPLANTLIST_CACHE = .theplantlist-cache.tar.gz
THEPLANTLIST_URL = https://github.com/richfitz/wood/releases/download/${RELEASE}/theplantlist-cache.tar.gz

DRYAD_CONTENTS = data/zae/Spermatophyta_Genera.csv \
  data/zae/PhylogeneticResources.zip \
  data/zae/GlobalWoodinessDatabase.csv
DRYAD_CACHE = .dryad-cache.tar.gz
DRYAD_URL = https://github.com/richfitz/wood/releases/download/${RELEASE}/dryad-cache.tar.gz

${THEPLANTLIST_CACHE}:
	curl -L -o $@ ${THEPLANTLIST_URL}
theplantlist-cache-fetch: ${THEPLANTLIST_CACHE}
theplantlist-delete:
	rm -rf ${THEPLANTLIST_CONTENTS}
theplantlist-cache-unpack: ${THEPLANTLIST_CACHE}
	tar -zxf $<
theplantlist-cache-archive:
	tar -zcf ${THEPLANTLIST_CACHE} ${THEPLANTLIST_CONTENTS}

${DRYAD_CACHE}:
	curl -L -o $@ ${DRYAD_URL}
dryad-cache-fetch: ${DRYAD_CACHE}
dryad-delete:
	rm -rf ${DRYAD_CONTENTS}
dryad-cache-unpack: ${DRYAD_CACHE}
	tar -zxf $<
dryad-cache-archive:
	tar -zcf ${DRYAD_CACHE} ${DRYAD_CONTENTS}

cache-fetch: theplantlist-cache-fetch dryad-cache-fetch
cache-unpack: theplantlist-cache-unpack dryad-cache-unpack

ARCHIVES = wood-supporting.tar.gz wood-analysis.tar.gz

# Packrat support.
PACKRAT_SOURCES_URL = https://github.com/richfitz/wood/releases/download/${RELEASE}/packrat.sources.tar.gz
include make/packrat.mk

release-files: ${ARCHIVES} ${DRYAD_CACHE} ${THEPLANTLIST_CACHE}
	rm -rf release
	mkdir -p release
	cp ${THEPLANTLIST_CACHE} release/theplantlist-cache.tar.gz
	cp ${DRYAD_CACHE} release/dryad-cache.tar.gz
	mv ${ARCHIVES} release
	cp ${PACKRAT_SOURCES_ARCHIVE} release

wood-supporting.tar.gz: ./make/wood-supporting.tar.gz.sh doc/wood-ms.pdf
	$<
wood-analysis.tar.gz: ./make/wood-analysis.tar.gz.sh doc/wood-ms.pdf
	$<

clean:
	rm -f data/theplantlist/names_accepted.csv
	rm -f data/zae/Vascular_Plants_rooted.dated.tre \
		data/zae/genus_order_lookup.csv
	rm -rf data/geo/raw
	rm -f output/*.rds output/*.csv
	rm -rf output/results
	rm -rf ${REPORT} figure cache
	make -C doc clean
	rm -f doc/figs/[a-z]*.pdf
	rm -f ${ARCHIVES}
	rm -rf release

purge: clean
	rm -rf ${THEPLANTLIST_CONTENTS}
	rm -f ${DRYAD_CONTENTS}

deps:
	${RSCRIPT} make/dependencies.R

# Lots of phony targets not listed...
.PHONY: all clean purge data-raw data-processed deps
