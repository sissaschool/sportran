# -*- coding: utf-8 -*-

################################################################################
###
###   ReadTABLEfile
###
################################################################################
###
###  a package that reads a table-style file and organizes it into a dictionary according to the column headers.
###  LAMMPS-style vector variables header are grouped together (only if group_vector = True).
###  If the name starts with "c_" or "v_", this is stripped away.
###    e.g.    c_flux[0] c_flux[1] c_flux[2]  -->  placed in 'flux' key
###  Comments lines are ignored.
###
###  Lines are read SEQUENTIALLY with the method read_datalines.
###  If a start_step is not specified the file is read from the current position.
###  This allows one to read the file in blocks.
###
################################################################################
###  The input file should look like this:
###
###  # COMMENT LINE
###  # COMMENT LINE
###  # COMMENT LINE
###  Step Temp TotEng Press c_flux[1] c_flux[2] c_flux[3] c_stress[1] c_stress[2]
###  0 257.6477 -1085.7346 -1944.803 -129.20254 124.70804 -200.42864 -64.236389 -134.0399
###  1 247.37505 -1085.734 -1909.333 -133.77141 124.25897 -103.27461 -61.022597 -83.17237
###  2 238.37359 -1087.9214 -1874.56 -138.58616 115.84038 -5.7728078 -58.471318 -74.51758
###  etc.
###
################################################################################
###   example:
###      current = TableFile(filename)
###      current.read_datalines(NSTEPS=100, start_step=0, select_ckeys=['Step', 'Temp', 'flux'])
###      print(current.data)
################################################################################

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


def data_length(filename):
    i = 0
    with open(filename) as f:
        while is_string(f.readline().split()[0]):   # skip text lines
            pass
        for i, l in enumerate(f, 2):
            pass
    return i


class TableFile(object):
    """
    A table-style file that can be read in blocks.

    Example:
      current = TableFile(filename)
      current.read_datalines(NSTEPS=100, select_ckeys=['Step', 'Temp', 'flux'])
      print(current.data)

    Variables (columns) are organized into a dictionary according to the column headers.
    LAMMPS-style vector variables header are grouped together (only if group_vector = True).
    If the name starts with "c_" or "v_", this is stripped away.
      e.g.    c_flux[0] c_flux[1] c_flux[2]  -->  placed in 'flux0' key
    Comments lines are ignored.

    #############################################################################
    The input file should look like this:

    # COMMENT LINE
    # COMMENT LINE
    # COMMENT LINE
    Step Temp TotEng Press c_flux[1] c_flux[2] c_flux[3] c_stress[1] c_stress[2]
    0 257.6477 -1085.7346 -1944.803 -129.20254 124.70804 -200.42864 -64.236389 -134.0399
    1 247.37505 -1085.734 -1909.333 -133.77141 124.25897 -103.27461 -61.022597 -83.17237
    2 238.37359 -1087.9214 -1874.56 -138.58616 115.84038 -5.7728078 -58.471318 -74.51758
    etc.
    """

    def __init__(self, *args, **kwargs):
        """ LAMMPS_Current(filename, select_ckeys) """
        if (len(args) > 0):
            self.filename = args[0]
            if (len(args) == 2):
                self.select_ckeys = args[1]
            else:
                self.select_ckeys = None
        else:
            raise ValueError('No file given.')
        group_vectors = kwargs.get('group_vectors', True)
        self._GUI = kwargs.get('GUI', False)
        self._print_elapsed = kwargs.get('print_elapsed', True)
        if self._GUI:
            from ipywidgets import FloatProgress
            from IPython.display import display
            global FloatProgress, display

        self._open_file()
        self._read_ckeys(group_vectors)
        self.ckey = None
        self.MAX_NSTEPS = data_length(self.filename)
        log.write_log('Data length = ', self.MAX_NSTEPS)
        return

    def __repr__(self):
        msg = 'TableFile:\n' + \
              '  filename:      {}\n'.format(self.filename) + \
              '  all_ckeys:     {}\n'.format(self.all_ckeys) + \
              '  select_ckeys:  {}\n'.format(self.select_ckeys) + \
              '  used ckey:     {}\n'.format(self.ckey) + \
              '  start pos:     {}\n'.format(self._start_byte) + \
              '  current pos:   {}\n'.format(self.file.tell())
        return msg

    def _open_file(self):
        """Open the file."""
        try:
            self.file = open(self.filename, 'r')
        except:
            raise ValueError('File does not exist.')
        return

    def _read_ckeys(self, group_vectors=True):
        """Read the column keys. If group_vectors=True the vector ckeys are grouped togheter"""
        self.all_ckeys = {}
        self.header = ''
        while True:
            line = self.file.readline()
            if len(line) == 0:   # EOF
                raise RuntimeError('Reached EOF, no ckeys found.')
            values = np.array(line.split())
            # text line: read variables names and save indexes in ckey
            if (is_string(values[0]) and (values[0].find('#') < 0)):
                self.header += line[:-1]
                for i in range(len(values)):
                    if group_vectors:
                        bracket = is_vector_variable(values[i])   # position of left square bracket
                    else:
                        bracket = 0
                    if (bracket == 0):   # the variable is a scalar
                        key = values[i]
                        if (key[:2] == 'c_'):   # remove 'c_' if present
                            key = key[2:]
                        self.all_ckeys[key] = [i]
                    else:   # the variable is a vector
                        key = values[i][:bracket]   # name of vector
                        if (key[:2] == 'c_'):   # remove 'c_' if present
                            key = key[2:]
                        vecidx = int(values[i][bracket + 1:-1])   # current index
                        if key in self.all_ckeys:   # if this vector is already defined, add this component
                            if (vecidx > self.all_ckeys[key].size):
                                self.all_ckeys[key] = np.resize(self.all_ckeys[key], vecidx)
                            self.all_ckeys[key][vecidx - 1] = i
                        else:   # if it is not, define a vector
                            self.all_ckeys[key] = np.array([0] * vecidx)
                            self.all_ckeys[key][-1] = i
                self._start_byte = self.file.tell()
                break
            else:
                self.header += line
        log.write_log(self.header)
        log.write_log(' #####################################')
        log.write_log('  all_ckeys = ', sorted(self.all_ckeys.items(), key=lambda kv: kv[0]))
        log.write_log(' #####################################')
        return

    def _set_ckey(self, select_ckeys=None, max_vector_dim=None):
        """Set the ckeys that have been selected, checking the available ones."""
        if select_ckeys is not None:
            self.select_ckeys = select_ckeys
        self.ckey = {}
        if self.select_ckeys is None:   # take all ckeys
            self.ckey = self.all_ckeys
        else:
            for key in self.select_ckeys:   # take only the selected ckeys
                value = self.all_ckeys.get(key, None)
                if value is not None:
                    self.ckey[key] = value[:max_vector_dim]   # copy all indexes (up to max dimension for vectors)
                else:
                    log.write_log('Warning: ', key, 'key not found.')
        if (len(self.ckey) == 0):
            raise KeyError('No ckey set. Check selected keys.')
        else:
            log.write_log('  ckey = ', sorted(self.ckey.items(), key=lambda kv: kv[0]))
        return

    def _initialize_dic(self, NSTEPS=None):
        """Initialize the data dictionary once the ckeys have been set."""
        if self.ckey is None:
            raise ValueError('ckey not set.')
        if NSTEPS is None:
            NSTEPS = 1
            self.dic_allocated = False
        else:
            self.dic_allocated = True
        self.NSTEPS = 0
        self.data = {}
        for key, idx in self.ckey.items():
            self.data[key] = np.zeros((NSTEPS, len(idx)))
        return

    def gotostep(self, start_step):
        """Go to the start_step-th line in the time series (assumes step=1).
         start_step = -1  -->  ignore, continue from current step
                       0  -->  go to start step
                       N  -->  go to N-th step"""
        if (start_step >= 0):
            self.file.seek(self._start_byte)
            for i in range(start_step):   # advance of start_step-1 lines
                self.file.readline()
        return

    def read_datalines(self, NSTEPS=0, start_step=-1, select_ckeys=None, max_vector_dim=None, even_NSTEPS=True):
        """Read NSTEPS steps of file, starting from start_step, and store only
      the selected ckeys.
      INPUT:
        NSTEPS         -> number of steps to read (default: 0 -> reads all the file)
        start_step  = -1 -> continue from current step (default)
                       0 -> go to start step
                       N -> go to N-th step
        select_ckeys   -> an array with the column keys you want to read (see all_ckeys for a list)
        max_vector_dim -> when reading vectors read only this number of components (None = read all components)
        even_NSTEPS    -> round the number of steps to an even number (default: True)
      OUTPUT:
        data    ->  a dictionary with the selected-column steps
      """
        if self._GUI:
            progbar = FloatProgress(min=0, max=100)
            display(progbar)
        start_time = time()
        if (NSTEPS == 0):
            NSTEPS = self.MAX_NSTEPS
        self._set_ckey(select_ckeys, max_vector_dim)   # set the ckeys to read
        self._initialize_dic(NSTEPS)   # allocate dictionary
        self.gotostep(start_step)   # jump to the starting step

        # read NSTEPS of the file
        progbar_step = max(100000, int(0.005 * NSTEPS))
        for step in range(NSTEPS):
            line = self.file.readline()
            if len(line) == 0:   # EOF
                log.write_log('Warning:  reached EOF.')
                break
            values = np.array(line.split())
            for key, idx in self.ckey.items():   # save the selected columns
                self.data[key][step, :] = np.array(list(map(float, values[idx])))
            if ((step + 1) % progbar_step == 0):
                if self._GUI:
                    progbar.value = float(step + 1) / NSTEPS * 100.
                    progbar.description = '{:6.2f}%'.format(progbar.value)
                else:
                    log.write_log('    step = {:9d} - {:6.2f}% completed'.format(step + 1,
                                                                                 float(step + 1) / NSTEPS * 100.))

        if self._GUI:
            progbar.close()
        # check number of steps read, keep an even number of steps
        if (step + 1 < self.NSTEPS):
            if (step == 0):
                log.write_log('WARNING:  no step read.')
                return
            else:
                log.write_log('Warning:  less steps read.')
                self.NSTEPS = step + 1
        if even_NSTEPS:
            if (NSTEPS % 2 == 1):
                NSTEPS = NSTEPS - 1
        for key, idx in self.ckey.items():   # free memory not used
            self.data[key] = self.data[key][:NSTEPS, :]
        log.write_log('  ( %d ) steps read.' % (NSTEPS))
        self.NSTEPS = NSTEPS
        if self._print_elapsed:
            log.write_log('DONE.  Elapsed time: ', time() - start_time, 'seconds')
        return self.data
