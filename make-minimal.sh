#!/bin/sh
rm -rf wood-minimal

mkdir wood-minimal
cp wood.* wood-functions.R wood-minimal

# Required data only:
mkdir wood-minimal/data
mkdir wood-minimal/data/geo
mkdir wood-minimal/data/zae

cp $(git ls-files --directory data)     wood-minimal/data
cp $(git ls-files --directory data/geo) wood-minimal/data/geo
cp $(git ls-files --directory data/zae) wood-minimal/data/zae
