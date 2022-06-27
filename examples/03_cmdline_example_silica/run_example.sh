#!/bin/bash

# if sportran was pip installed use:
SPORTRAN="sportran-analysis"
# otherwise use:
#SPORTRAN="../../sportran/analysis.py"

${SPORTRAN} ../data/Silica/Silica.dat --input-format table -k flux1 -C heat -u metal -t 1.0 --VOLUME 3130.431110818 --param-from-input-file-column Temp TEMPERATURE -w 0.1 --FSTAR 28.0 -r
