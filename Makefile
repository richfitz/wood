DATA = data/zae/genus_order_lookup.csv data/theplantlist/names_accepted.csv

data: ${DATA}

data/zae/genus_order_lookup.csv: data/zae/download.R
	Rscript $<

data/theplantlist/names_accepted.csv: data/theplantlist/build.R
	Rscript $<
