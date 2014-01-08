#!/bin/sh

BASE=wood
FILE=wood.R
STR="library(sowsear);sowsear(\"$FILE\", \"Rmd\");knit(\"${FILE}md\")"
Rscript -e "$STR"

pandoc $BASE.md -o $BASE.html --standalone --highlight-style=tango
pandoc $BASE.md -o $BASE.pdf
