#!/bin/sh

# Build the supporting information
DEST=wood-analysis
rm -rf $DEST
mkdir $DEST
cp -R figure wood.Rmd wood.md wood.html $DEST

tar -zcf $DEST.tar.gz $DEST
rm -rf $DEST
