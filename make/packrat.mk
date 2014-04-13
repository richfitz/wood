# Because packrat isn't already complicated enough I wanted to be able
# to turn it on and off.  So that typing something like
#
#   make packrat-enable
#   make packrat-disable
#
# will toggle whether packrat is enabled or not.
#
# This is a bit peculiar because it doesn't really fit into the use
# case for packrat (in which projects are completely isolated from one
# another always).  Instead, I want to document a set of working
# packages but also just presume that version rot won't terribly
# affect us.  An experiment, really.
#
# So the basic idea is that someone cloning the project gets just our
# files without any package management.  If they run into trouble then
# they could
#   export USE_PACKRAT=1
# and then that would do all its buisness installing things.
#
# This is probably really fragile!  And packrat feels very new at the
# moment, so it's not totally clear that this really helps anything.
# However, the bootstrap script is nice, and at the least we'll end up
# with a big mess of sources that people can work through.
#
# What would be nice is if we could just point at the lockfile and say
# "download these things".  Will be interesting to see how these
# change over the next year or so.
#
# Note that This will clobber the .Rprofile and .Renviron files in the
# project directory!
#
# I think that rbundler would have also been a decent solution, but
# suffers similar issues, plus the bootstraping is less graceful.  But
# the package is simpler!

PACKRAT_INIT = "library(methods);if(exists('initPackrat')) initPackrat()"
PACKRAT_SOURCES_ARCHIVE = .packrat/packrat.sources.tar.gz

packrat-enable: ${PACKRAT_SOURCES_ARCHIVE}
	cp .packrat/.Renviron .packrat/.Rprofile .packrat/packrat.lock .
	tar -zxf ${PACKRAT_SOURCES_ARCHIVE}
	Rscript -e ${PACKRAT_INIT}
packrat-disable:
	rm -f .Renviron .Rprofile packrat.lock
	rm -rf packrat.sources
packrat-purge: packrat-disable
	rm -rf library
packrat-update:
	cp .Renviron .Rprofile packrat.lock .packrat
	tar -zcf ${PACKRAT_SOURCES_ARCHIVE} packrat.sources
packrat-refresh: packrat-purge
	Rscript -e "library(methods); packrat::bootstrap()"
	Rscript -e ${PACKRAT_INIT}
	make packrat-update

${PACKRAT_SOURCES_ARCHIVE}:
	curl -L -o $@ ${PACKRAT_SOURCES_URL}
packrat-fetch: ${PACKRAT_SOURCES_ARCHIVE}
packrat-fetch-force:
	rm -f ${PACKRAT_SOURCES_ARCHIVE}
	make packrat-fetch

packrat-perhaps:
ifeq ($(USE_PACKRAT),1)
	make packrat-enable
endif
