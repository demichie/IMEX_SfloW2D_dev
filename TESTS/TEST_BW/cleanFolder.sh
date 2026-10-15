#!/bin/bash

echo "Cleaning folder form output of previous runs of IMEX_SfloW2D"

rm -rf OUTPUT/TEST_BW_400

file="TEST_BW*"

rm -f $file

file="*bak"

rm -f $file

file="dem_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

file="dem_interfaces_esri.asc"

if [ -f $file ] ; then
    rm $file
fi
