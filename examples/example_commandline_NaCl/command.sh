#!/bin/bash

# if thermocepstrum was pip installed use:
thermocepstrum-analysis ../data/NaCl.dat -k flux -j 'vcm[1]' -t 5.0 -V 65013.301261 -w 0.1 --FSTAR 14.0 -r

# if thermocepstrum is not pip installed use:
#../../thermocepstrum/analysis.py ../data/NaCl.dat -k flux -j 'vcm[1]' -t 5.0 -V 65013.301261 -w 0.1 --FSTAR 14.0 -r
