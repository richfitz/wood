#!/bin/sh

# Build the supporting information
DEST=wood-supporting
rm -rf $DEST
mkdir $DEST
mkdir -p $DEST/output
cp -R output/results $DEST/output

# Copy our data over, leaving external data behind.  This should be
# easier!
mkdir $DEST/data
tar -cf - `git ls-files --directory data` | tar -xpf - -C $DEST

tar -zcf $DEST.tar.gz $DEST
rm -rf $DEST
