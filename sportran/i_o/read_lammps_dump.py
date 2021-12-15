# -*- coding: utf-8 -*-

################################################################################
###
###   ReadLAMMPSDump - v0.1.8 - May 03, 2018
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
from sportran.utils import log


def is_string(string):
    try:
        float(string)
    except ValueError:
        return True
    return False


def is_vector_variable(string):
    bracket = string.rfind('[')
    if (bracket == -1):
        bracket = 0
    return bracket


def file_length(filename):
    i = -1
    with open(filename) as f:
        for i, l in enumerate(f, 1):
            pass
    return i


def get_volume(filename):
    f = open(filename, 'r')
    line = f.readline()
    while (line):
        if 'BOX BOUNDS' in line:
            xlo, xhi = list(map(float, f.readline().split()))
            ylo, yhi = list(map(float, f.readline().split()))
            zlo, zhi = list(map(float, f.readline().split()))
            break
        line = f.readline()
    f.close()
    volume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo)
    return volume


def get_natoms(filename):
    f = open(filename, 'r')
    line = f.readline()
    while (line):
        if 'NUMBER OF ATOMS' in line:
            natoms = int(f.readline())
            break
        line = f.readline()
    f.close()
    return natoms


class LAMMPS_Dump(object):
    """
    A LAMMPS_Dump file that can be read in blocks.
    example:
      traj = LAMMPS_Dump(filename, preload=False)  -->> do not preload list of steps (suggested if the file is big)
      traj.read_timesteps(10, start_step=0, select_ckeys=['id,xu,yu,vu']) -->>   Read first 10 timesteps, only the specified columns
      traj.read_timesteps(10, select_ckeys=['id,xu,yu,vu']) -->>   Read the next 10 timesteps, only the specified columns (DELTA_TIMESTEP is assumed)
      traj.read_timesteps((10,30))      -->>  Read from TIMESTEP 10 to 30
      traj.read_timesteps((10,30,2))    -->>  Read every 2 steps from TIMESTEP 10 to 30
      print(traj.data)
    """

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
        self._quiet = kwargs.get('quiet', False)
        self._GUI = kwargs.get('GUI', False)
        if self._GUI:
            from ipywidgets import FloatProgress
            from IPython.display import display
            global FloatProgress, display

        self._open_file()
        self._read_ckeys(group_vectors, preload_timesteps)
        self.ckey = None
        #self.MAX_NSTEPS = data_length(self.filename)
        #log.write_log("Data length = ", self.MAX_NSTEPS)
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
              '  start pos:     {}\n'.format(self._start_byte) + \
              '  current pos:   {}\n'.format(self.file.tell()) + \
              '  FIRST TIMESTEP: {}\n'.format(self.FIRST_TIMESTEP) + \
              '  LAST TIMESTEP:  {}\n'.format(self.LAST_TIMESTEP) + \
              '  DELTA TIMESTEP: {}\n'.format(self.DELTA_TIMESTEP) + \
              '  current step:   {}\n'.format(self.current_timestep)
        return msg

    def _open_file(self):
        """Open the file."""
        try:
            self.file = open(self.filename, 'r')
        except:
            raise ValueError('File does not exist.')
        return

    def _read_ckeys(self, group_vectors=True, preload_timesteps=True):
        """Read the column keys. If group_vectors=True the vector ckeys are grouped togheter"""
        self._start_byte = self.file.tell()
        self.all_ckeys = {}
        self.all_timesteps = []
        self.preload_timesteps = preload_timesteps
        while True:
            line = self.file.readline()
            if len(line) == 0:   # EOF
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
                    for i in range(2, len(values)):
                        if group_vectors:
                            bracket = is_vector_variable(values[i])   # get position of left square bracket
                        else:
                            bracket = 0
                        if (bracket == 0):   # the variable is a scalar
                            key = values[i]
                            if (key[:2] == 'c_'):   # remove 'c_' if present
                                key = key[2:]
                            self.all_ckeys[key] = [i - 2]   # -2 offset
                        else:   # the variable is a vector
                            key = values[i][:bracket]   # name of vector
                            if (key[:2] == 'c_'):   # remove 'c_' if present
                                key = key[2:]
                            vecidx = int(values[i][bracket + 1:-1])   # current index
                            if key in self.all_ckeys:   # if this vector is already defined, add this component
                                if (vecidx > self.all_ckeys[key].size):
                                    self.ckeys[key] = np.resize(self.all_ckeys[key], vecidx)
                                self.all_ckeys[key][vecidx - 1] = i - 2   # -2 offset!
                            else:   # if it is not, define a vector
                                self.all_ckeys[key] = np.array([0] * vecidx)
                                self.all_ckeys[key][-1] = i - 2   # -2 offset!
                    #self._start_byte = self.file.tell()
                    break
            #else:
            #   self.header += line

        if self.preload_timesteps:
            # get the list of time steps
            while True:
                line = self.file.readline()
                if len(line) == 0:   # EOF
                    break
                if (line == 'ITEM: TIMESTEP\n'):
                    self.current_timestep = int(self.file.readline())
                    self.all_timesteps.append(self.current_timestep)

            self.LAST_TIMESTEP = self.all_timesteps[-1]
            self.DELTA_TIMESTEP = self.all_timesteps[1] - self.FIRST_TIMESTEP
            self.TOT_TIMESTEPS = len(self.all_timesteps)
            self.all_timesteps = np.array(self.all_timesteps)
        else:
            log.write_log(' ** No timesteps pre-loaded. Be careful in the selection. **')
            # get the first 2 timesteps
            while (len(self.all_timesteps) < 2):
                line = self.file.readline()
                if len(line) == 0:   # EOF
                    break
                if (line == 'ITEM: TIMESTEP\n'):
                    self.current_timestep = int(self.file.readline())
                    self.all_timesteps.append(self.current_timestep)

            self.LAST_TIMESTEP = None
            self.DELTA_TIMESTEP = self.all_timesteps[1] - self.FIRST_TIMESTEP
            self.TOT_TIMESTEPS = None
            self.all_timesteps = None

        # go back to the first timestep
        self.gototimestep(0)   # compute_first = True
        self._start_byte = 0
        log.write_log('  all_ckeys      = ', self.all_ckeys)
        log.write_log('  TOT_TIMESTEPS  = ', self.TOT_TIMESTEPS)
        log.write_log('  FIRST_TIMESTEP = ', self.FIRST_TIMESTEP)
        log.write_log('  DELTA_TIMESTEP = ', self.DELTA_TIMESTEP)
        log.write_log('  LAST_TIMESTEP  = ', self.LAST_TIMESTEP)
        log.write_log('  all_timesteps  = ', self.all_timesteps)
        return

    def _set_ckey(self, select_ckeys=None):
        """
        Set the ckeys to read from the selected, checking the available ones.
        If select_ckeys is not passed, then use the already selected ones, or all the available ones if no selection
        was previously made.
        """
        if select_ckeys is not None:
            self.select_ckeys = select_ckeys
        self.ckey = {}
        if self.select_ckeys is None:   # take all ckeys
            self.ckey = self.all_ckeys
        else:
            for key in self.select_ckeys:   # take only the selected ckeys
                value = self.all_ckeys.get(key, None)
                if value is not None:
                    self.ckey[key] = value[:]   # copy all indexes (up to max dimension for vectors)
                else:
                    log.write_log('Warning: ', key, 'key not found.')
        if (len(self.ckey) == 0):
            raise KeyError('No ckey set. Check selected keys.')
        else:
            if not self._quiet:
                log.write_log('  ckey = ', self.ckey)
        return

    def _set_timesteps(self, selection, start_step=-1):
        """Set the timesteps to read from the selected, checking the available ones.
      INPUT:  N              -->  Read the next N steps (DELTA_TIMESTEP is assumed)
              N, start_step=30  -->  Read N steps from the TIMESTEP 30
                                  if compute_first=True, read the current step as well
              (10,30)        -->  Read from TIMESTEP 10 to 30
              (10,30,2)      -->  Read every 2 steps from TIMESTEP 10 to 30"""
        if (start_step == -1):
            if self._compute_current_step:
                start_step = self.current_timestep
            else:
                start_step = self.current_timestep + self.DELTA_TIMESTEP
        elif (start_step == 0):
            start_step = self.FIRST_TIMESTEP
        if np.isscalar(selection) or (len(selection) == 1):   # select N steps from start one
            first = start_step
            last = self.DELTA_TIMESTEP * selection + start_step
            step = None
        elif (len(selection) == 2):
            first = selection[0]
            last = selection[1]
            step = None
        elif (len(selection) == 3):
            first = selection[0]
            last = selection[1]
            step = selection[2]
        if step is None:
            step = self.DELTA_TIMESTEP
        elif (step % self.DELTA_TIMESTEP != 0):
            log.write_log('Warning: step is not a multiple of the detected DELTA_TIMESTEP. You may get errors.')
        if (first % step != 0):
            first += step - first % step   # round first step to the next in the list

        self.timestep = []
        self.select_timesteps = np.arange(first, last, step)   # selected timesteps
        if self.preload_timesteps:
            for step in self.select_timesteps:
                if step in self.all_timesteps:
                    self.timestep.append(step)   # make list of available selected-timesteps
                else:
                    log.write_log('Warning: timestep # {:d} not found.'.format(step))
        else:
            self.timestep = self.select_timesteps   # use all the selected (be careful)
        self.nsteps = len(self.timestep)   # number of available steps
        if (self.nsteps == 0):
            raise ValueError('No timestep set. Check selected timesteps.')
        else:
            if not self._quiet:
                log.write_log('  nsteps   = ', self.nsteps)
                log.write_log('  timestep = ', self.timestep)
        return

    def _initialize_dic(self):
        """Initialize the data dictionary once the ckeys and timesteps have been set."""
        if self.ckey is None:
            raise ValueError('ckey not set.')
        if self.timestep is None:
            raise ValueError('timestep not set.')
        self.data = [dict() for i in range(self.nsteps)]
        for istep in range(self.nsteps):
            for key, idx in self.ckey.items():
                if (key == 'element'):   # this should be improved
                    self.data[istep][key] = np.zeros((self.NATOMS, len(idx)), dtype='S8')
                else:
                    self.data[istep][key] = np.zeros((self.NATOMS, len(idx)), dtype='float64')
        return

    def _gototimestep(self, start_step, fast_check=True):
        """
        Go to the start_step-th line in the time series (assumes step=1).
          start_step = -1  -->  ignore, continue from current step
                        0  -->  go to FIRST timestep
                        N  -->  go to N-th timestep
          fast_check = True --> assumes the TIMESTEP are a monotonously increasing.
                                If the the start_step is passed and not found then stop.
        """
        if (start_step >= 0):
            if (start_step <= self.current_timestep):
                # or (self.current_timestep == -1):  # if start_step is before/equal the current step
                self.file.seek(self._start_byte)   #  --> start over
            if (start_step == 0):   # or (self.current_timestep == -1):
                goto_step = self.FIRST_TIMESTEP
            else:
                goto_step = start_step

            # search until start_step is found     ***** MAY BE IMPROVED KNOWING THE N OF LINES TO SKIP ******
            while True:
                line = self.file.readline()
                if len(line) == 0:   # EOF
                    raise EOFError('Warning (gototimestep):  reached EOF. Timestep {} NOT FOUND.'.format(goto_step))
                if (line == 'ITEM: TIMESTEP\n'):
                    self.current_timestep = int(self.file.readline())
                    if (self.current_timestep == goto_step):
                        while (self.file.readline().find('ITEM: ATOMS') < 0):   # jump to the data part
                            pass
                        break
                    if (fast_check) and (self.current_timestep > goto_step):
                        raise Warning(
                            'Warning (gototimestep):  Timestep {} NOT FOUND up to current_step = {}. (To force check the whole trajectory set fast_check=False)'
                            .format(goto_step, self.current_timestep))
        else:
            pass
        return

    def gototimestep(self, start_step, fast_check=True):
        """
        Go to the start_step-th line in the time series (assumes step=1).
          start_step = -1  -->  ignore, continue from current step
                        0  -->  go to FIRST timestep
                        N  -->  go to N-th timestep
          fast_check = True --> assumes the TIMESTEP are a monotonously increasing.
                                If the the start_step is passed and not found then stop.
        """
        ## user-called function
        self._compute_current_step = True
        self._gototimestep(start_step, fast_check)
        return

    def read_timesteps(self, selection, start_step=-1, select_ckeys=None, fast_check=True):
        """
        Read selected keys of file, within the provided range.
        Examples:
            read_timesteps(10, start_step=0, select_ckeys=['id,xu,yu,vu']) -->>   Read first 10 timesteps, only the specified columns
            read_timesteps(10, select_ckeys=['id,xu,yu,vu']) -->>   Read the next 10 timesteps, only the specified columns (DELTA_TIMESTEP is assumed)
            read_timesteps((10,30))      -->>  Read from TIMESTEP 10 to 30
            read_timesteps((10,30,2))    -->>  Read every 2 steps from TIMESTEP 10 to 30
        """
        if self._GUI:
            progbar = FloatProgress(min=0, max=100)
            display(progbar)
        start_time = time()
        self._set_ckey(select_ckeys)   # set the ckeys to read      --> ckey
        self._set_timesteps(selection, start_step)   # set the timesteps to read  --> timestep
        self._initialize_dic()   # allocate dictionary        --> data

        # extract the steps from the file
        progbar_step = max(1000, int(0.005 * self.nsteps))
        atomid_col = self.all_ckeys['id'][0]
        for istep, step in enumerate(self.timestep):
            self._gototimestep(step, fast_check)   # jump to the desired step,
            self.data[istep]['TIMESTEP'] = step
            for nat in range(self.NATOMS):   # read data (may be unsorted)
                line = self.file.readline()
                if len(line) == 0:   # EOF
                    raise EOFError('Warning:  reached EOF.')
                values = np.array(line.split())
                for key, idx in self.ckey.items():   # save the selected columns
                    atomid = int(values[atomid_col]) - 1   # current atom index (in LAMMPS it starts from 1)
                    if (key == 'element'):   # this should be improved
                        self.data[istep][key][atomid, :] = np.array(list(map(str, values[idx])))
                    else:
                        self.data[istep][key][atomid, :] = np.array(list(map(float, values[idx])))
            if ((istep + 1) % progbar_step == 0):
                if self._GUI:
                    progbar.value = float(istep + 1) / self.nsteps * 100.
                    progbar.description = '%g %%' % progbar.value
                else:
                    log.write_log('    step = {:9d} - {:6.2f}% completed'.format(istep + 1,
                                                                                 float(istep + 1) / self.nsteps * 100.))
        if self._GUI:
            progbar.close()
        # check number of steps read, keep an even number of steps
        if (istep + 1 < self.nsteps):   # (should never happen)
            if (istep == 0):
                log.write_log('WARNING:  no step read.')
                return
            else:
                log.write_log('Warning:  less steps read.')
                self.nsteps = istep + 1
        if not self._quiet:
            log.write_log('  ( %d ) steps read.' % (self.nsteps))
            log.write_log('DONE.  Elapsed time: ', time() - start_time, 'seconds')
        self._compute_current_step = False   # next time do not compute the current_step
        return self.data
