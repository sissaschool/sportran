#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
###
###   ReadLAMMPSLogFile
###
################################################################################
###
###  A package that reads a LAMMPS Log file and organizes it into a dictionary according to the column headers.
###  LAMMPS-style vector variables header are grouped together (only if group_vector = True).
###  If the name starts with "c_" or "v_", this is stripped away.
###    e.g.    c_flux[0] c_flux[1] c_flux[2]  -->  placed in 'flux0' key
###
###  The output of a specific run command is identified by a string 'run_keyword'.
###  When this keyword is found the next output block is read (it is supposed to
###  begin with a 'Step' in the first column).
###  The reading stops when and 'end_keyword' is found (default: 'Loop time').
###
###  Lines are read SEQUENTIALLY with the method read_datalines.
###  If a start_step is not specified the file is read from the current position.
###  This allows one to read the file in blocks.
###
################################################################################
###  Example of LAMMPS Log file:
###
###  fix NVE all nve
###  # PRODUCTION RUN
###  run 1000
###  Per MPI rank memory allocation (min/avg/max) = 4.45 | 4.456 | 4.458 Mbytes
###  ...
###  ...
###  Step Temp TotEng Press c_flux[1] c_flux[2] c_flux[3] c_stress[1] c_stress[2]
###  0 257.6477 -1085.7346 -1944.803 -129.20254 124.70804 -200.42864 -64.236389 -134.0399
###  1 247.37505 -1085.734 -1909.333 -133.77141 124.25897 -103.27461 -61.022597 -83.17237
###  2 238.37359 -1087.9214 -1874.56 -138.58616 115.84038 -5.7728078 -58.471318 -74.51758
###  ...
###  Loop time of 110.158 on 20 procs for 400000 steps with 1728 atoms
###
################################################################################
###  Example script:
###     data = LAMMPSLogFile('lammps.log', run_keyword='PRODUCTION RUN')
###     data.read_datalines(NSTEPS=100, start_step=0, select_ckeys=['Step', 'Temp', 'flux'])
###     print(data.data)
###
###     # to save data into a Numpy file:
###     save_hc_npz(data, ['flux'], 'lammps.data', 'flux.npz')
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


class LAMMPSLogFile(object):
    """
  A package that reads a LAMMPS Log file and organizes it into a dictionary according to the column headers.
  LAMMPS-style vector variables header are grouped together (only if group_vector = True).
  If the name starts with "c_" or "v_", this is stripped away.
    e.g.    c_flux[0] c_flux[1] c_flux[2]  -->  placed in 'flux0' key

  The output of a specific run command is identified by a string 'run_keyword'.
  When this keyword is found the next output block is read (it is supposed to
  begin with a 'Step' in the first column).
  The reading stops when and 'end_keyword' is found (default: 'Loop time').

  Lines are read SEQUENTIALLY with the method read_datalines.
  If a start_step is not specified the file is read from the current position.
  This allows one to read the file in blocks.

#############################################################################
  Example of LAMMPS Log file:

  fix NVE all nve
  # PRODUCTION RUN
  run 1000
  Per MPI rank memory allocation (min/avg/max) = 4.45 | 4.456 | 4.458 Mbytes
  ...
  ...
  Step Temp TotEng Press c_flux[1] c_flux[2] c_flux[3] c_stress[1] c_stress[2]
  0 257.6477 -1085.7346 -1944.803 -129.20254 124.70804 -200.42864 -64.236389 -134.0399
  1 247.37505 -1085.734 -1909.333 -133.77141 124.25897 -103.27461 -61.022597 -83.17237
  2 238.37359 -1087.9214 -1874.56 -138.58616 115.84038 -5.7728078 -58.471318 -74.51758
  ...
  Loop time of 110.158 on 20 procs for 400000 steps with 1728 atoms

#############################################################################
  Example script:
     data = LAMMPSLogFile(filename, run_keyword='PRODUCTION RUN')
     data.read_datalines(NSTEPS=100, start_step=0, select_ckeys=['Step', 'Temp', 'flux'])
     print(data.data)

     # to save data into a Numpyz file:
     save_hc_npz(data, ['flux'], 'lammps.data', 'flux.npz')
#############################################################################
    """

    def __init__(self, *args, **kwargs):
        """LAMMPSLogFile(filename, run_keyword, select_ckeys)"""
        if (len(args) > 0):
            self.filename = args[0]
            if (len(args) >= 2):
                self.run_keyword = args[1]
                if (len(args) == 3):
                    self.select_ckeys = args[2]
                else:
                    self.select_ckeys = None
            else:
                self.run_keyword = None
        else:
            raise ValueError('No file given.')
        self.run_keyword = kwargs.get('run_keyword', None)
        self.endrun_keyword = kwargs.get('endrun_keyword', 'Loop time')
        group_vectors = kwargs.get('group_vectors', True)
        self._GUI = kwargs.get('GUI', False)
        if self.run_keyword is None:
            raise ValueError('Please specify run_keyword.')
        if self._GUI:
            from ipywidgets import FloatProgress
            from IPython.display import display
            global FloatProgress, display

        self._open_file()
        self.MAX_NSTEPS = file_length(self.filename)
        self._read_ckeys(self.run_keyword, group_vectors)
        self.ckey = None
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

    def _read_ckeys(self, run_keyword, group_vectors=True):
        """Seek the line containing 'run_keyword'. Read the column keys. If group_vectors=True the vector ckeys are grouped togheter."""
        self.all_ckeys = {}
        nlines = 0
        while True:
            line = self.file.readline()
            nlines += 1
            if len(line) == 0:   # EOF
                raise RuntimeError('Reached EOF, run_keyword was not found!')
            # check if run_keyword string is found
            if run_keyword in line:
                log.write_log('  run_keyword found at line {:d}.'.format(nlines))
                break
        while True:
            line = self.file.readline()
            nlines += 1
            if len(line) == 0:   # EOF
                raise RuntimeError('Reached EOF, no ckeys found.')
            values = np.array(line.split())
            # find the column headers line
            if (len(values) and (is_string(values[0])) and (values[0] == 'Step')):
                log.write_log('  column headers found at line {:d}. Reading data...'.format(nlines))
                for i in range(len(values)):
                    if group_vectors:
                        bracket = is_vector_variable(values[i])   # position of left square bracket
                    else:
                        bracket = 0
                    if (bracket == 0):   # the variable is a scalar
                        key = values[i]
                        if (key[:2] == 'c_'):   # remove 'c_' if present
                            key = key[2:]
                        self.all_ckeys[key] = np.array([i])
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
                self.MAX_NSTEPS -= nlines
                break
        self.NALLCKEYS = np.concatenate(list(self.all_ckeys.values())).size
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
        """
        Go to the start_step-th line in the time series (assumes step=1).
          start_step = -1  -->  ignore, continue from current step
                        0  -->  go to start step
                        N  -->  go to N-th step
        """
        if (start_step >= 0):
            self.file.seek(self._start_byte)
            for i in range(start_step):   # advance of start_step-1 lines
                self.file.readline()
        return

    def read_datalines(self, NSTEPS=0, start_step=-1, select_ckeys=None, max_vector_dim=None, even_NSTEPS=True):
        """
        Read NSTEPS steps of file, starting from start_step, and store only the selected ckeys.

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
            if (len(line) == 0):   # EOF
                log.write_log('Warning:  reached EOF.')
                break
            if self.endrun_keyword in line:   # end-of-run keyword
                log.write_log('  endrun_keyword found.')
                step -= 1
                break
            values = np.array(line.split())
            if (values.size != self.NALLCKEYS):
                log.write_log('Warning:  line with wrong number of columns found. Stopping here...')
                log.write_log(line)
                break
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
        if (step + 1 < NSTEPS):
            if (step == 0):
                log.write_log('WARNING:  no step read.')
                return
            else:
                if (NSTEPS != self.MAX_NSTEPS):   # if NSTEPS was specified
                    log.write_log('Warning:  less steps read.')
                NSTEPS = step + 1   # the correct number of read steps
        # even the number of steps
        if even_NSTEPS:
            if (NSTEPS % 2 == 1):
                NSTEPS = NSTEPS - 1
                log.write_log('  Retaining an even number of steps (even_NSTEPS=True).')
        for key, idx in self.ckey.items():   # free the memory not used
            self.data[key] = self.data[key][:NSTEPS, :]
        log.write_log('  ( %d ) steps read.' % (NSTEPS))
        self.NSTEPS = NSTEPS
        log.write_log('DONE.  Elapsed time: ', time() - start_time, 'seconds')
        return self.data


def save_hc_npz(lammpslogfile, select_ckeys, lammps_structurefilename, outfilename):
    """
     Takes a LAMMPSLogFile object, a LAMMPS structure data file (optional), takes
     the desired columns and save data into a Numpyz file.

     # example to save data into a Numpyz file:
     Save_HC_npz(lammpslogfile_object, ['flux'], 'lammps.data', 'flux.npz')
   """

    def get_box(filename):
        """Return the box edges and volume from a LAMMPS Data file."""
        box = np.zeros((2, 3))
        with open(filename, 'r') as f:
            while True:
                line = f.readline().split()
                if 'xlo' in line:
                    box[:, 0] = [float(line[0]), float(line[1])]
                if 'ylo' in line:
                    box[:, 1] = [float(line[0]), float(line[1])]
                if 'zlo' in line:
                    box[:, 2] = [float(line[0]), float(line[1])]
                    break
        volume = np.prod(box[1, :] - box[0, :])
        return box, volume

    if not isinstance(lammpslogfile, LAMMPSLogFile):
        raise ValueError('lammpslogfile is not a LAMMPSLogFile object.')

    dic = {}
    if 'Temp' not in lammpslogfile.ckey:
        raise RuntimeError('Temp not found.')
    dic['Temp_ave'] = np.mean(lammpslogfile.data['Temp'])
    dic['Temp_std'] = np.std(lammpslogfile.data['Temp'])

    for key in select_ckeys:
        if key in lammpslogfile.ckey:
            dic[key] = lammpslogfile.data[key]
        else:
            raise ValueError('ckey not found.')

    box, volume = get_box(lammps_structurefilename)
    dic['box'] = box
    dic['Volume'] = volume

    dic['DT'] = lammpslogfile.data['Step'][1, 0] - lammpslogfile.data['Step'][0, 0]
    if 'Time' in lammpslogfile.ckey:
        dic['DT_TIMEUNITS'] = lammpslogfile.data['Time'][1, 0] - lammpslogfile.data['Time'][0, 0]

    log.write_log('These keys will be saved in file \"{:}\" :'.format(outfilename))
    log.write_log(' ', list(dic.keys()))
    np.savez(outfilename, **dic)
    return


################################################################################
def main():
    """
    This script extracts the desired columns from a LAMMPS Log file and saves them into a Numpyz file for later use.

    Example:
    start reading when "PRODUCTION RUN" is found, read "flux1" and "Press" columns from log.lammps
    log file and structure.data data file (containing structure):
       python read_lammps_log.py  log.lammps structure.data out.npz -k flux1 Press -d "PRODUCTION RUN"
    """

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('lammps_logfile', help='lammps log file to read')
    parser.add_argument('lammps_structurefile', help='lammps structure data file to read')
    parser.add_argument('output_file', help='numpy file that will be saved')
    parser.add_argument('-k', '--ckeys', nargs='+', dest='ckeys', help='list of column keys to read')
    parser.add_argument('-d', '--runkey', help='RUN keyword')
    parser.add_argument('-e', '--endrunkey', help='ENDRUN keyword (default: "Loop time")')
    args = parser.parse_args()

    logfile = LAMMPSLogFile(args.lammps_logfile, run_keyword=args.runkey)
    select_ckeys = args.ckeys[:]
    if 'Step' not in args.ckeys:
        select_ckeys += ['Step']
    if 'Time' not in args.ckeys:
        if 'Time' in logfile.all_ckeys:
            select_ckeys += ['Time']
    if 'Temp' not in args.ckeys:
        select_ckeys += ['Temp']
    logfile.read_datalines(select_ckeys=select_ckeys)
    save_hc_npz(logfile, args.ckeys, args.lammps_structurefile, args.output_file)
    return 0


if __name__ == '__main__':
    main()
