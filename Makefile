DATA_RAW = data/zae/genus_order_lookup.csv data/theplantlist/names_accepted.csv
DATA_PROCESSED = output/woodiness.rds output/dat.g.rds \
	output/dat.g.w.rds output/dat.g.h.rds output/phy.o.rds
REPORT = wood.Rmd wood.md wood.html

all: wood.html doc/wood-ms.pdf

data-raw: ${DATA_RAW}
data-processed: ${DATA_PROCESSED}

wood.Rmd: wood.R
	Rscript -e "library(sowsear); sowsear('$<', 'Rmd')"
wood.md: wood.Rmd ${DATA_PROCESSED}
	Rscript -e "library(knitr); knit('$<')"
wood.html: wood.md
	Rscript -e "library(markdown);\
	 opts <- setdiff(markdownHTMLOptions(TRUE), 'base64_images');\
	 markdownToHTML('$<', '$@', options=opts)"

doc/wood-ms.pdf: wood.md
	make -C doc

# This target will never run because it depends on nothing.  But if
# data/geo/country_coords.csv is deleted then this will regenerate
# that file.  It's here for reference only, really.
data/geo/country_coords.csv:
	Rscript make/data-geo-country_coords.csv.R

data/zae/genus_order_lookup.csv: make/data-zae.R
	Rscript $<

data/theplantlist/names_accepted.csv: make/data-theplantlist.R
	Rscript $<

output/genus_order_lookup.csv: make/output-genus_order_lookup.csv.R ${DATA_RAW}
	Rscript $<

output/woodiness.rds: make/output-woodiness.rds.R R/load.R R/build.R ${DATA_RAW}
	Rscript $<

DATA_GENUS_DEPS = output/genus_order_lookup.csv output/woodiness.rds \
	R/load.R R/build.R

output/dat.g.rds: make/output-dat.g.rds.R ${DATA_GENUS_DEPS}
	Rscript $<
output/dat.g.w.rds: make/output-dat.g.w.rds.R ${DATA_GENUS_DEPS}
	Rscript $<
output/dat.g.h.rds: make/output-dat.g.h.rds.R ${DATA_GENUS_DEPS}
	Rscript $<

output/phy.o.rds: make/output-phy.o.rds.R output/dat.g.rds
	Rscript $<

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
	rm -f wood-supporting.zip

# Save on some farting about with data:
DOWNLOADED_DATA =               	   \
	data/theplantlist/acceptedNames1.1 \
        data/zae/PhylogeneticResources.zip \
	data/zae/Spermatophyta_Genera.csv  \
	data/zae/GlobalWoodinessDatabase.csv
DOWNLOADED_DATA_SAVE = .downloaded_data.tar.gz

downloaded-data-delete:
	rm -fr ${DOWNLOADED_DATA}
downloaded-data-save:
	tar -zcf ${DOWNLOADED_DATA_SAVE} ${DOWNLOADED_DATA}
downloaded-data-unpack:
	tar -zxf ${DOWNLOADED_DATA_SAVE}
downloaded-data-bulk-fetch:
	curl -o ${DOWNLOADED_DATA_SAVE} http://www.zoology.ubc.ca/~fitzjohn/files/wood_data.tar.gz

wood-supporting.zip: ./make/wood-supporting.zip.sh
	$<

.PHONY: all clean data-raw data-processed
