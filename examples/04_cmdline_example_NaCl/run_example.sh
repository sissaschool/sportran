#!/bin/bash

# if sportran was pip installed use:
SPORTRAN="sportran-analysis"
# otherwise use:
#SPORTRAN="../../sportran/analysis.py"

${SPORTRAN} ../data/NaCl/NaCl.dat --input-format table -k flux -j 'vcm[1]' -C heat -u metal -t 5.0 --VOLUME 65013.301261 --param-from-input-file-column Temp TEMPERATURE -w 0.1 --FSTAR 14.0 -r
