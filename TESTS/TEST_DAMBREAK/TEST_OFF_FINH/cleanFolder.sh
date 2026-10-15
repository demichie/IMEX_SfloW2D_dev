#!/bin/bash

echo "Cleaning folder form output of previous runs of IMEX_SfloW2D"

rm -rf OUTPUT/dambreak2D

file="dambreak2D*"

rm -f $file
