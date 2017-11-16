################################################################################
###
###   ReadLAMMPSDump - v0.1.5 - August 28th, 2017
###
################################################################################
###
###  a package to read LAMMPS Dump files
###  (it assumes that the data column names and the number of atoms do not change)
### 
################################################################################

## example:
##   import read_lammps_dump as rd
##   data = rd.LAMMPS_Dump(filename)
##    

import numpy as np
from time import time

def is_string( string ):
  try:
    float( string )
  except ValueError:
    return True
  return False

def is_vector_variable( string ):
  bracket = string.rfind('[')
  if (bracket == -1):
    bracket = 0
  return bracket

def file_length( filename ):
  i = -1
  with open(filename) as f:
    for i, l in enumerate(f,1):
      pass
  return i

def get_volume( filename ):
  f = open(filename, 'r')
  line = f.readline()
  while( line ):
    if "BOX BOUNDS" in line:
      xlo, xhi = map(float, f.readline().split())
      ylo, yhi = map(float, f.readline().split())
      zlo, zhi = map(float, f.readline().split())
      break
    line = f.readline()
  f.close()
  volume = (xhi-xlo)*(yhi-ylo)*(zhi-zlo)
  return volume

def get_natoms( filename ):
  f = open(filename, 'r')
  line = f.readline()
  while( line ):
    if "NUMBER OF ATOMS" in line:
      natoms = int(f.readline())
      break
    line = f.readline()
  f.close()
  return natoms


class LAMMPS_Dump(object):
   """A LAMMPS_Dump file that can be read in blocks."""

   def __init__(self, *args, **kwargs):
      #*******
      if (len(args) > 0):
         self.filename = args[0]
         if (len(args) == 2):
            self.select_ckeys = args[1]
         else:
            self.select_ckeys = None
      else:
         raise ValueError('No file given.')
      group_vectors = kwargs.get('group_vectors', True)
      preload_timesteps = kwargs.get('preload', True)
      self.GUI = kwargs.get('GUI', False)
      if self.GUI:
         from ipywidgets import FloatProgress
         from IPython.display import display
         global FloatProgress, display
      
      self.open_file()
      self.read_ckeys(group_vectors, preload_timesteps)
      self.ckey = None
      #self.MAX_NSTEPS = data_length(self.filename)
      #print "Data length = ", self.MAX_NSTEPS
      return

   def __repr__(self):
      msg = 'LAMMPS_Dump:\n' + \
            '  filename:      {}\n'.format(self.filename) + \
            '  all_ckeys:     {}\n'.format(self.all_ckeys) + \
            '  select_ckeys:  {}\n'.format(self.select_ckeys) + \
            '  used ckey:     {}\n'.format(self.ckey) + \
            '  all_timesteps:    {}\n'.format(self.all_timesteps) + \
            '  select_timesteps: {}\n'.format(self.select_timesteps) + \
            '  used timesteps:   {}\n'.format(self.timestep) + \
            '  start pos:     {}\n'.format(self.start_byte) + \
            '  current pos:   {}\n'.format(self.file.tell()) + \
            '  FIRST TIMESTEP: {}\n'.format(self.FIRST_TIMESTEP) + \
            '  LAST TIMESTEP:  {}\n'.format(self.LAST_TIMESTEP) + \
            '  DELTA TIMESTEP: {}\n'.format(self.DELTA_TIMESTEP) + \
            '  current step:   {}\n'.format(self.current_timestep)
      return msg
      
   def open_file(self):
      """Open the file."""
      try:
         self.file = open(self.filename, 'r')
      except:
         raise ValueError('File does not exist.')
      return


   def read_ckeys(self, group_vectors=True, preload_timesteps=True):
      """Read the column keys. If group_vectors=True the vector ckeys are grouped togheter"""
      self.start_byte = self.file.tell()
      self.all_ckeys = {}
      self.all_timesteps = []
      self.preload_timesteps = preload_timesteps
      while True:
         line = self.file.readline()
         if len(line) == 0:  # EOF
            raise RuntimeError('Reached EOF, no ckeys found.')
         values = np.array(line.split())
         if (values[0] == 'ITEM:'):
            if (values[1] == 'TIMESTEP'):
               self.current_timestep = int(self.file.readline())
               self.FIRST_TIMESTEP = self.current_timestep
               self.all_timesteps.append(self.current_timestep)
            # facoltativo:
            elif ((values[1] == 'NUMBER') and values[2] == 'OF' and values[3] == 'ATOMS'):
               self.NATOMS = int(self.file.readline())
            elif ((values[1] == 'BOX') and values[2] == 'BOUNDS'):
               self.BOX_BOUNDS_TYPE = values[3:6]
               xbox = self.file.readline().split()
               ybox = self.file.readline().split()
               zbox = self.file.readline().split()
               self.BOX_BOUNDS = np.array([xbox, ybox, zbox], dtype='float')
            elif (values[1] == 'ATOMS'):
               for i in xrange(2, len(values)):
                  if group_vectors:
                     bracket = is_vector_variable( values[i] )  # get position of left square bracket
                  else:
                     bracket = 0
                  if (bracket == 0):     # the variable is a scalar
                     key = values[i]
                     if (key[:2] == 'c_'):    # remove 'c_' if present
                        key = key[2:]
                     self.all_ckeys[key] = [i - 2]   # -2 offset
                  else:                  # the variable is a vector
                     key = values[i][:bracket]       # name of vector
                     if (key[:2] == 'c_'):    # remove 'c_' if present
                        key = key[2:]
                     vecidx = int(values[i][bracket+1:-1])  # current index
                     if key in self.all_ckeys:     # if this vector is already defined, add this component
                        if (vecidx > self.all_ckeys[key].size):
                           self.ckeys[key] = np.resize(self.all_ckeys[key], vecidx)
                        self.all_ckeys[key][vecidx-1] = i - 2   # -2 offset!
                     else:               # if it is not, define a vector
                        self.all_ckeys[key] = np.array( [0]*vecidx )
                        self.all_ckeys[key][-1] = i - 2   # -2 offset!
               #self.start_byte = self.file.tell()
               break
         #else:
         #   self.header += line

      if self.preload_timesteps:
         # get the list of time steps
         while True:
            line = self.file.readline()
            if len(line) == 0:  # EOF
               break
            if (line == 'ITEM: TIMESTEP\n'):
               self.current_timestep = int(self.file.readline())
               self.all_timesteps.append(self.current_timestep)

         self.LAST_TIMESTEP = self.all_timesteps[-1]
         self.DELTA_TIMESTEP = self.all_timesteps[1] - self.FIRST_TIMESTEP
         self.TOT_TIMESTEPS = self.all_timesteps.size
         self.all_timesteps = np.array(self.all_timesteps)
      else:
         print " ** No timesteps pre-loaded. Be careful in the selection. **"
         # get the first 2 timesteps
         while (len(self.all_timesteps) < 2):
            line = self.file.readline()
            if len(line) == 0:  # EOF
               break
            if (line == 'ITEM: TIMESTEP\n'):
               self.current_timestep = int(self.file.readline())
               self.all_timesteps.append(self.current_timestep)

         self.LAST_TIMESTEP = None
         self.DELTA_TIMESTEP = self.all_timesteps[1] - self.FIRST_TIMESTEP
         self.TOT_TIMESTEPS = None
         self.all_timesteps = None


      # go back to the first timestep
      self.gototimestep(0)
      self.start_byte = 0
      print '  all_ckeys      = ', self.all_ckeys
      print '  TOT_TIMESTEPS  = ', self.TOT_TIMESTEPS
      print '  FIRST_TIMESTEP = ', self.FIRST_TIMESTEP
      print '  DELTA_TIMESTEP = ', self.DELTA_TIMESTEP
      print '  LAST_TIMESTEP  = ', self.LAST_TIMESTEP
      print '  all_timesteps  = ', self.all_timesteps
      return


   def set_ckey(self, select_ckeys=None):
      """Set the ckeys to read from the selected, checking the available ones."""
      if select_ckeys is not None:
         self.select_ckeys = select_ckeys
      self.ckey = {}
      if self.select_ckeys is None:   # take all ckeys
         self.ckey = self.all_ckeys
      else:
         for key in self.select_ckeys:  # take only the selected ckeys
            value = self.all_ckeys.get(key, None)
            if value is not None:
               self.ckey[key] = value[:]  # copy all indexes (up to max dimension for vectors)
            else:
               print "Warning: ", key, "key not found."
      if (len(self.ckey) == 0):
         raise KeyError("No ckey set. Check selected keys.")
      else:
         print "  ckey = ", self.ckey
      return


   def set_timesteps(self, selection, start_step=0):
      """Set the timesteps to read from the selected, checking the available ones.
      INPUT:   scalar, start  -->  Number of steps to read from start index"""      ### MODIFICA CON INPUT: vettore [first,last,step] oppure scalare [NSTEPS]
      if (start_step == 0):
         start_step = self.FIRST_TIMESTEP
      if np.isscalar(selection) or (len(selection) == 1):   # select N steps from start
         first = start_step
         last  = self.DELTA_TIMESTEP*selection + self.FIRST_TIMESTEP
         step  = None
      elif (len(selection) == 2):
         first = selection[0]
         last  = selection[1]
         step  = None
      elif (len(selection) == 3):
         first = selection[0]
         last  = selection[1]
         step  = selection[2]
      if step is None:
         step = self.DELTA_TIMESTEP
      self.timestep = []
      self.select_timesteps = np.arange(first, last, step) # selected timesteps
      if self.preload_timesteps:
         for step in self.select_timesteps:
            if step in self.all_timesteps:
               self.timestep.append(step)    # make list of available selected-timesteps
            else:
               print "Warning: timestep # {:d} not found.".format(step)
      else:
         self.timestep = self.select_timesteps   # select all (be careful)
      self.nsteps = len(self.timestep)    # number of available steps
      if (self.nsteps == 0):
         raise ValueError("No timestep set. Check selected timesteps.")
      else:
         print "  nsteps   = ", self.nsteps
         print "  timestep = ", self.timestep
      return


   def initialize_dic(self):
      """Initialize the data dictionary once the ckeys and timesteps have been set."""
      if self.ckey is None:
         raise ValueError('ckey not set.')
      if self.timestep is None:
         raise ValueError('timestep not set.')
      self.data = [dict() for i in xrange(self.nsteps)]
      for istep in xrange(self.nsteps):
         for key, idx in self.ckey.iteritems():
            self.data[istep][key] = np.zeros( (self.NATOMS, len(idx)) )
      return


   def gototimestep(self, start_step):
      """Go to the start_step-th line in the time series (assumes step=1).
         start_step = -1  -->  ignore, continue from current step
                       0  -->  go to FIRST timestep
                       N  -->  go to N-th timestep"""
      if (start_step >= 0):
         if (start_step <= self.current_timestep):  # if start_step is before/equal the current step
            self.file.seek(self.start_byte)         #  --> start over
         if (start_step == 0):
            start_step = self.FIRST_TIMESTEP

         # search until start_step is found     ***** MAY BE IMPROVED KNOWING THE N OF LINES TO SKIP ******
         while True:
            line = self.file.readline()
            if len(line) == 0:  # EOF
               raise Warning("Warning (gototimestep):  reached EOF. Timestep NOT FOUND.")
            if (line == 'ITEM: TIMESTEP\n'):
               self.current_timestep = int(self.file.readline())
               if (self.current_timestep == start_step):
                  while (self.file.readline().find('ITEM: ATOMS') < 0):  # jump to the data part
                     pass
                  break
      return

   
   def read_timesteps(self, selection, select_ckeys=None):
      """Read selected keys of file, within the provided range.
      e.g.  read_steps(4, 10, 2)  reads steps n. 4, 6, 8"""
      if self.GUI:
         progbar = FloatProgress(min=0, max=100)
         display(progbar)
      start_time = time()
      self.set_ckey(select_ckeys)                  # set the ckeys to read      --> ckey
      self.set_timesteps(selection)                # set the timesteps to read  --> timestep
      self.initialize_dic()                        # allocate dictionary        --> data

      # extract the steps from the file
      progbar_step = max(1000, int(0.005*self.nsteps))
      atomid_col = self.all_ckeys['id'][0]
      for istep,step in enumerate(self.timestep):
         self.gototimestep(step)                   # jump to the desired step
         self.data[istep]['TIMESTEP'] = step
         for nat in xrange(self.NATOMS):           # read data (may be unsorted)
            line = self.file.readline()
            if len(line) == 0:   # EOF
               raise Warning("Warning:  reached EOF.")
            values = np.array(line.split())
            for key, idx in self.ckey.iteritems():   # save the selected columns
               atomid = int(values[atomid_col]) - 1  # current atom index (in LAMMPS it starts from 1)
               self.data[istep][key][atomid,:] = np.array(map(float, values[idx]))
         if ( (istep+1)%progbar_step == 0 ):
            if self.GUI:
               progbar.value = float(istep+1)/self.nsteps*100.;
               progbar.description = "%g %%" % progbar.value
            else:
               print "    step = {:9d} - {:6.2f}% completed".format(istep+1, float(istep+1)/self.nsteps*100.)
      if self.GUI:
         progbar.close()
      # check number of steps read, keep an even number of steps
      if (istep + 1 < self.nsteps):      # (should never happen)
         if (istep == 0):
            print "WARNING:  no step read."
            return
         else:
            print "Warning:  less steps read."
            self.nsteps = istep + 1
      print "  ( %d ) steps read." % (self.nsteps)
      print "DONE.  Elapsed time: ", time()-start_time, "seconds"
      return self.data

