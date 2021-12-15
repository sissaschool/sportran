#!/bin/bash

# if sportran was pip installed use:
sportran-analysis ../data/NaCl.dat -k flux -j 'vcm[1]' -t 5.0 -V 65013.301261 -w 0.1 --FSTAR 14.0 -r

# if sportran is not pip installed use:
#../../sportran/analysis.py ../data/NaCl.dat -k flux -j 'vcm[1]' -t 5.0 -V 65013.301261 -w 0.1 --FSTAR 14.0 -r
