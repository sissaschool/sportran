#!/bin/bash

# if sportran was pip installed use:
sportran-analysis ../data/NaCl.dat -k flux -j 'vcm[1]' -t 5.0 --VOLUME 65013.301261 --param-from-input-file-column Temp TEMPERATURE -w 0.1 --FSTAR 14.0 -r

# if sportran is not pip installed use:
#../../sportran/analysis.py ../data/NaCl.dat -k flux -j 'vcm[1]' -t 5.0 --VOLUME 65013.301261 --param-from-input-file-column Temp TEMPERATURE -w 0.1 --FSTAR 14.0 -r
