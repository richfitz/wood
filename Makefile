DATA_RAW = data/zae/genus_order_lookup.csv data/theplantlist/names_accepted.csv
DATA_PROCESSED = output/woodiness.rds output/dat.g.rds \
	output/dat.g.w.rds output/dat.g.h.rds

all: wood.html wood.pdf doc/wood-ms.pdf

data-raw: ${DATA_RAW}
data-processed: ${DATA_PROCESSED}

wood.Rmd: wood.R
	Rscript -e "library(sowsear); sowsear('wood.R', 'Rmd')"
wood.md: wood.Rmd ${DATA_PROCESSED}
	Rscript -e "library(knitr); knit('wood.Rmd')"
wood.html: wood.md
	pandoc wood.md -o wood.html --standalone --highlight-style=tango
wood.pdf: wood.md
	pandoc wood.md -o wood.pdf

doc/wood-ms.pdf: wood.pdf
	make -C doc wood-ms.pdf

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
output/dat.g.h.rds: make/output-dat.g.w.rds.R ${DATA_GENUS_DEPS}
	Rscript $<

clean:
	rm -f output/*.rds output/*.csv
	make -C doc clean

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

.PHONY: all clean data-raw data-processed
