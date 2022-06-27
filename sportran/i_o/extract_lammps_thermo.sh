#!/usr/bin/env bash

function usage
{
   echo -e "This script extracts the thermo log from a lammps log file.
It assumes that a line such as
   dump DUMP_RUN_KEY all custom file.lammpstrj id type x y z
or simply
   #dump DUMP_RUN_KEY
is present in the file before the 'run' command.
DUMP_RUN can be set with the -d flag.
Note: DUMP_RUN cannot contain spaces

Usage:  extract_current.sh -f LOG_FILE [-h]
         -f LOG_FILE LAMMPS log file name
         -d DUMP_RUN_KEY label - used to identify the production run
         --start-time    initial time to add to the Time column
         --start-step    initial step to add to the Step column
         --last-step     last step to read
         --skip-first    skip header and first step (useful to concatenate files)
         -h  --help      print this help
Example:  extract_current.sh -f water.log -d DUMP_RUN > thermo.dat
"
}

LOG_FILE=
DUMP_RUN_KEY='DUMP_RUN'
VERBOSITY=1
START_STEP=0
START_TIME=0
SKIP_FIRST=0
LAST_STEP=-1

while [ $# -gt 0 ]; do
    case $1 in
        -f )             shift; LOG_FILE=$1
                         ;;
        -d )             shift; DUMP_RUN_KEY=$1
                         ;;
        -v )             VERBOSITY=0
                         ;;
        --start-time )   shift; START_TIME=$1
                         ;;
        --start-step )   shift; START_STEP=$1
                         ;;
        --last-step )    shift; LAST_STEP=$1
                         ;;
        --skip-first )   SKIP_FIRST=1
                         ;;
        -h | --help )    usage
                         exit
                         ;;
        * )              usage
                         exit 1
    esac
    shift
done

if [[ ( (-z "$LOG_FILE") ) ]]; then
   usage
   exit 1
fi

echo " Reading file:  "${LOG_FILE} > "/dev/stderr"
echo " DUMP RUN key:  "${DUMP_RUN_KEY} > "/dev/stderr"

awk -vDUMPKEY=${DUMP_RUN_KEY} -vVERB=${VERBOSITY} -vSTART_TIME=${START_TIME} -vSTART_STEP=${START_STEP} -vSKIP_FIRST=${SKIP_FIRST} -vLAST_STEP=${LAST_STEP} '
BEGIN {
  flag = 0       # flag = 1 read column headers, flag = 2 read data
  nfluxcol = 0
  stepcol = 0
  timecol = 0
  count = 0
  if (LAST_STEP >= 0 )
    laststepflag = 1
  else
    laststepflag = 0
}
(($1 == "dump" || $1 == "#dump") && $2 == DUMPKEY) || (($2 == "dump" && $1 == "#") && $3 == DUMPKEY){
  ## production run found
  printf(" dump %s  found at line %d.\n", DUMPKEY, NR) > "/dev/stderr"
  flag = 1
}
flag == 1 && $1 == "Step" {
  ## here you can read column headers
  flag = 2
  for (i=1; i<=NF; i++) {
    printf("Found: %16s at line %d, column %d\n", $i, NR, i) > "/dev/stderr"
    if ($i == "Step") {
      stepcol = i
    }
    else if ($i == "Time") {
      timecol = i
    }
  }
  if (!SKIP_FIRST)
    print $0
  next
}
flag == 2 {
  if (count==0 && SKIP_FIRST) {
    count++
    next
  }
  else {
    ## read data until Loop
    if ($1 == "Loop") {
      flag = 0
      exit 0
    }
    else if ($1 == "WARNING:") {
      flag = 0
      printf("WARNING found at line %d.\n", NR) > "/dev/stderr"
    }
    else if (laststepflag && ($stepcol > LAST_STEP)) {  # Last step reached
      exit 0
    }
    else {      ## print all
      for (i=1; i<=NF; i++) {
        if (i == stepcol) {
          printf("%d ", $i + START_STEP)
        }
        else if (i == timecol) {
          printf("%g ", $i + START_TIME)
        }
        else {
          printf("%s ", $i)
        }
      }
      printf("\n")
      count++
    }
  }
}
END {
  if ( count == 0 )
    printf("ERROR: no steps found.\n") > "/dev/stderr"
  else
    printf("DONE:  %d steps of data read.\n", count) > "/dev/stderr"
}
' $LOG_FILE
