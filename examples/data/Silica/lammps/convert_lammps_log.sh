#!/usr/bin/env bash

# Convert the silica.out LAMMPS log file into:
#  - Silica.dat :  a TableFile (a table-style text file with data sorted in columns)
#  - Silica.npy :  a NumPy binary pickle file

LAMMPS_FILE='silica.out'
SPORTRAN_DIR='../../../../sportran'

# convert to table format
${SPORTRAN_DIR}/i_o/extract_lammps_thermo.sh -f ${LAMMPS_FILE} -d 'DUMP_RUN' > Silica.dat

# convert to numpy binary
python ${SPORTRAN_DIR}/i_o/read_lammps_log.py ${LAMMPS_FILE} Silica.npy -s silica_216_1000K.init -d 'NVE RUN' 
