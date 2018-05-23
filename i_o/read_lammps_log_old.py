################################################################################
###
###   ReadLAMMPSlogfile
###
################################################################################

import numpy as np
from thermocepstrum.i_o import __path__
from subprocess import check_output


def is_vector_variable( string ):
   bracket = string.rfind('[')
   if (bracket == -1):
      bracket = 0
   return bracket


def Read_LAMMPSlogfile(filename, DUMP_RUN='DUMP_RUN', group_vectors=True, even_NSTEPS=True, GUI=False):
   """example:
   data = read_LAMMPSlogfile('log.lammps', '_dump_run')
   """
   if GUI:
      from ipywidgets import FloatProgress
      from IPython.display import display
      global FloatProgress, display

   def extract_thermo():
      script = __path__[0] + '/extract_lammps_thermo.sh'
      out = check_output([script, '-f', filename, '-d', DUMP_RUN]).splitlines()
      all_ckeys = read_ckeys(out[0])  # read column header
      return out, all_ckeys


   def read_ckeys(strarray):
      """Read the column keys. If group_vectors=True the vector ckeys are grouped togheter"""
      all_ckeys = {}
      values = np.array(strarray.split())
      for i in range(len(values)):
         if group_vectors:
            bracket = is_vector_variable( values[i] )  # position of left square bracket
         else:
            bracket = 0
         if (bracket == 0):     # the variable is a scalar
            key = values[i]
            if (key[:2] == 'c_'):    # remove 'c_' if present
               key = key[2:]
            all_ckeys[key] = [i]
         else:                  # the variable is a vector
            key = values[i][:bracket]       # name of vector
            if (key[:2] == 'c_'):    # remove 'c_' if present
               key = key[2:]
            vecidx = int(values[i][bracket+1:-1])  # current index
            if key in all_ckeys:     # if this vector is already defined, add this component
               if (vecidx > all_ckeys[key].size):
                  all_ckeys[key] = np.resize(all_ckeys[key], vecidx)
               all_ckeys[key][vecidx-1] = i
            else:               # if it is not, define a vector
               all_ckeys[key] = np.array( [0]*vecidx )
               all_ckeys[key][-1] = i
      print ' #####################################'
      print '  all_ckeys = ', all_ckeys
      print ' #####################################'
      return all_ckeys


   def initialize_dic(NSTEPS, ckey):
      """Initialize the data dictionary once the ckeys have been set."""
      data = {}
      for key, idx in ckey.iteritems():
         data[key] = np.zeros( (NSTEPS, len(idx)) )
      return data


   def read_datalines(filedata, ckey, even_NSTEPS=True, GUI=False):
      if GUI:
         progbar = FloatProgress(min=0, max=100)
         display(progbar)
      NSTEPS = len(filedata) - 1
      data = initialize_dic(NSTEPS, ckey)

      progbar_step = max(100000, int(0.005*NSTEPS))
      for step, line in enumerate(filedata[1:]):
         if len(line) == 0:  # EOF
            print "Warning:  reached EOF."
            break
         values = np.array(line.split())
         for key, idx in ckey.iteritems():  # save the selected columns
            data[key][step,:] = np.array(map(float, values[idx]))
         if ( (step+1)%progbar_step == 0 ):
            if GUI:
               progbar.value = float(step+1)/NSTEPS*100.;
               progbar.description = "{:6.2f}%".format(progbar.value)
            else:
               print "    step = {:9d} - {:6.2f}% completed".format(step+1, float(step+1)/NSTEPS*100.)
      if GUI:
         progbar.close()
      # check number of steps read, keep an even number of steps
      if (step + 1 < NSTEPS):
         if (step == 0):
            print "WARNING:  no step read."
            return
         else:
            print "Warning:  less steps read."
            NSTEPS = step + 1
      if even_NSTEPS:
         if (NSTEPS%2 == 1):
            NSTEPS = NSTEPS - 1
      for key, idx in ckey.iteritems():  # free memory not used
         data[key] = data[key][:NSTEPS,:]
      print "  ( %d ) steps read." % (NSTEPS)
      return data


   outdata, ckey = extract_thermo()
   data = read_datalines(outdata, ckey, even_NSTEPS, GUI)
   return data

