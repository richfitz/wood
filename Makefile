DATA_RAW = data/zae/genus_order_lookup.csv data/theplantlist/names_accepted.csv
DATA_PROCESSED = output/woodiness.rds output/dat.g.rds \
	output/dat.g.w.rds output/dat.g.h.rds

data-raw: ${DATA_RAW}
data-processed: ${DATA_PROCESSED}

data/zae/genus_order_lookup.csv: data/zae/download.R
	Rscript $<

data/theplantlist/names_accepted.csv: data/theplantlist/build.R
	Rscript $<

output/woodiness.rds: ${DATA} wood-functions.R
	Rscript R/make-output-woodiness.rds.R

output/dat.g.rds: output/woodiness.rds wood-functions.R
	Rscript R/make-output-dat.g.rds.R
output/dat.g.w.rds: output/woodiness.rds wood-functions.R
	Rscript R/make-output-dat.g.w.rds.R
output/dat.g.h.rds: output/woodiness.rds wood-functions.R
	Rscript R/make-output-dat.g.h.rds.R

clean:
	rm -f output/*.rds

.PHONY: data-raw data-processed
