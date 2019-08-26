# -*- coding: future_fstrings -*-

"""
--------------------------------------------
    Thermocepstrum graphic user interface

    Control unit file
--------------------------------------------

This file contains the functions that the GUI will directly call to start new
multithread operations and to make requests to the thermocepstrum calculus unit.
"""

import os
from . import settings

import thermocepstrum as tc
import numpy as np
from thermocepstrum_gui.utils import Graph

try:
    from thermocepstrum.utils.utils import PrintMethod
except ImportError:
    from thermocepstrum_gui.utils.utils import PrintMethod

# Init print method
log = PrintMethod()


# -------- DATA SECTION --------

class Data:
    """
    This class is used as a structure to
    contains all the variables that the software use.
    """

    def __init__(self):
        self.CURRENT_FILE = ''
        self.loaded = False
        self.changes = False

        self.inputformat = None
        self.jfile = None
        self.jdata = None
        self.j = None
        self.xf = None
        self.axis = None

        self.keys = None
        self.description = None

        self._temperature = 0.0
        self.temperature_std = 0.0
        self._volume = 0
        self._DT_FS = 0.0
        self.currents = None

        self._fstar = 0.0
        self.old_fstar = 0.0
        self._psd_filter_width = 0.1
        self._units = 'metal'

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        self._verify_changes(value, self.temperature)
        self._temperature = value
        print(self.changes, ' temperature')

    @property
    def volume(self):
        return self._volume

    @volume.setter
    def volume(self, value):
        self._verify_changes(value, self.volume)
        self._volume = value
        print(self.changes, ' volume')

    @property
    def DT_FS(self):
        return self._DT_FS

    @DT_FS.setter
    def DT_FS(self, value):
        self._verify_changes(value, self.DT_FS)
        self._DT_FS = value
        print(self.changes, ' DT_FS')

    @property
    def psd_filter_width(self):
        return self._psd_filter_width

    @psd_filter_width.setter
    def psd_filter_width(self, value):
        self._verify_changes(value, self.psd_filter_width)
        self._psd_filter_width = value
        print(self.changes, ' filter')

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, value):
        self._verify_changes(value, self.units)
        self._units = value
        print(self.changes, ' units')

    @property
    def fstar(self):
        return self._fstar

    @fstar.setter
    def fstar(self, value):
        self._verify_changes(value, self.fstar)
        self._fstar = value
        print(self.changes, ' units')

    def _verify_changes(self, val1, val2):
        if not self.changes:
            if val1 != val2:
                self.changes = True

# todo: add header property
# todo: check if the behaviour is correct

try:
    data
except:
    data = Data()

# -------- GRAPH SECTION --------


"""
This section contains the functions that deal with
the graph manager.
"""


# Setup graph manager
gm = Graph.GraphManager()


def set_graph(axis_, func, **kwargs):
    """
    This function request to the graph manager to generate
    a new graph to be displayed in the axis_ screen.

    :param axis_: the location to draw the new graph.
    :param func: the function to generate the graph.
    :param kwargs: other arguments to pass to the file manager.
    :return axis: the new graph.
    """

    axis = func(axis=axis_, external_object=data, **kwargs)
    return axis


# -------- FILES AND LOAD SECTION --------


"""
In this section there are the functions that are used to load the data
that the software will use and some functions to work with files.

The starter settings are loaded from the .ini file that contains the
values that need to be stored.
"""


def get_file_size(path):
    """
    This function calculate the file size.

    :param path: the location of the file.
    :return: the file size.
    """
    file_size = os.path.getsize(path)
    if file_size >= 1000000:
        file_size //= 1000000
        return f"{file_size} KB"
    else:
        file_size //= 1000
        return f"{file_size} MB"


def load_path():
    """
    This function loads from the .ini file the paths
    where the file will be stored.
    """

    # Verify that the file exists otherwise generate it
    if not os.path.exists('thcp.ini'):
        set_defaults()

    # Read the file
    with open('thcp.ini', 'r') as file:
        settings_file = file.read().split('\n')

        # Load the paths
        if settings_file[0]:
            for setting in settings_file:
                if setting:
                    var, val = setting.split(':')
                    if var == 'DP':
                        settings.DATA_PATH = val
                        if not os.path.exists(settings.DATA_PATH):
                            os.mkdir(settings.DATA_PATH)
                    elif var == 'OP':
                        settings.OUTPUT_PATH = val
                        if not os.path.exists(settings.OUTPUT_PATH):
                            os.mkdir(settings.OUTPUT_PATH)
                    elif var == 'LP':
                        settings.LOG_PATH = val
                        if not os.path.exists(settings.LOG_PATH):
                            os.mkdir(settings.LOG_PATH)
        else:
            set_defaults()


def set_defaults():
    """
    This function is used to set the default values of the .ini file
    """
    # Create the file
    with open('thcp.ini', 'w') as file:
        defaults = (('DP', 'data'),
                    ('OP', 'outputs'),
                    ('LP', 'logs'))

        # Write the defaults paths
        for setting in defaults:
            try:
                os.mkdir(os.path.join(settings.BASE_PATH, setting[1]))
            except:
                pass
            file.write(f"{setting[0]}:{os.path.join(settings.BASE_PATH, setting[1])}\n")

    load_path()


def load_keys(inputfile):
    """
    This function is used to load the header keys of the selected file.

    :param inputfile: the path of the selected file.
    :return:
    """

    jfile = tc.i_o.TableFile(inputfile, group_vectors=True)
    return jfile.all_ckeys


def load_data(inputfile, input_format, _selected_keys, temperature=None, NSTEPS=0, START_STEP=0,
              run_keyword='', units=None, DT_FS=None, volume=None, psd_filter_w=None, axis_=None,
              structurefile=None, _descriptions=[]):

    selected_keys = _selected_keys.copy()
    descriptions = _descriptions.copy()
    data.temperature = temperature
    data.volume = volume
    data.DT_FS = DT_FS
    data.inputformat = input_format
    data.psd_filter_width = psd_filter_w

    if input_format == 'table':
        jfile = tc.i_o.TableFile(inputfile, group_vectors=True)
        data.jfile = jfile
        jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
        data.jdata = jfile.data
    elif input_format == 'dict':
        data.jdata = np.load(inputfile)
    elif input_format == 'lammps':
        jfile = tc.i_o.LAMMPSLogFile(inputfile, run_keyword=run_keyword)
        jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
        data.jdata = jfile.data
    else:
        raise NotImplemented('input format not implemented.')

    ## Define currents
    # print(selected_keys, jindex)

    if descriptions.count('Temperature') == 1:
        temperature = get_temp(data.jdata, get_cor_index(selected_keys, descriptions, 'Temperature')[0])
        data.temperature = temperature
        i = get_cor_index(selected_keys, descriptions, 'Temperature')[1]
        del descriptions[i]
        del selected_keys[i]
    if volume is None:
        volume = get_volume(data.jdata, structurefile)

    if NSTEPS == 0:
        NSTEPS = data.jdata[list(data.jdata.keys())[0]].shape[0]
    if True:    # jindex is None:
        heat_current, i = get_cor_index(selected_keys, descriptions, 'Energy current')
        del descriptions[i]
        del selected_keys[i]
        currents_headers = [heat_current]
        for other in selected_keys:
            currents_headers.append(other)

        currents = np.array([data.jdata[key][START_STEP:(START_STEP + NSTEPS), :] for key in currents_headers])
    else:
        pass
        # if sindex is None:
        #     currents = np.array([Data.jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] for key in selected_keys])
        # else:
        #     currents = np.array([Data.jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] - Data.jdata[key][START_STEP:(
        #                 START_STEP + NSTEPS), sindex] for key in selected_keys])
    data.currents = currents
    # create HeatCurrent object
    emsgs = []
    if volume is not -1:
        if temperature is not -1:
            data.j = tc.heatcurrent.HeatCurrent(currents, units, DT_FS, temperature, volume, psd_filter_w)
            gm.initialize(data.j)
        else:
            emsgs.append('Invalid temperature!')
            return -1, emsgs
    else:
        emsgs.append('Invalid volume!')
        return -1, emsgs


# -------- UTILITY SECTION --------


"""
This section contains some utility functions.
"""


def secure_exit(main_window):
    '''
    This function allow a secure exit from the execution
    of the software.
    '''

    # todo: Add multithread processes safe exit and a log file
    for process in main_window.open_windows:
        try:
            process.destroy()
        except AttributeError:
            process.master.destroy()
        except:
            raise SystemExit('Unable to exit the program safely!')

    exit()


def get_cor_index(arr1, arr2, corr_key):
    i = arr2.index(corr_key)
    return arr1[i], i


def get_temp(jdata, selected_key):
    if selected_key in jdata:
        temperature = np.mean(jdata[selected_key])
        temperature_std = np.std(jdata[selected_key])   # this is wrong (needs block average)
        print(' Mean Temperature (computed):  {} K  +/-  {}\n'.format(temperature, temperature_std))
    # elif 'Temp_ave' in jdata:
    #     temperature = jdata['Temp_ave']
    #     if 'Temp_std' in jdata:
    #         temperature_std = jdata['Temp_std']
    #         print(' Mean Temperature (file):      {} K  +/-  {}'.format(temperature, temperature_std))
    #         # logfile.write(' Mean Temperature (file):      {} K  +/-  {}\n'.format(temperature, temperature_std))
    #     else:
    #         print(' Mean Temperature (file):      {} K'.format(temperature))
    #         # logfile.write(' Mean Temperature (file):      {} K\n'.format(temperature))
    else:
        temperature = -1
        # raise RuntimeError('No Temp key found. Please provide Temperature (-T).')

    return temperature


def get_volume(jdata, structurefile):
    if structurefile is not None:
        _, volume = tc.i_o.read_lammps_datafile.get_box(structurefile)
        print(' Volume (structure file):    {} A^3'.format(volume))
        # logfile.write(' Volume (structure file):    {} A^3'.format(volume))
    elif 'Volume' in jdata:
        volume = jdata['Volume']
        print(' Volume (file):    {} A^3'.format(volume))
        # logfile.write(' Volume (file):    {} A^3\n'.format(volume))
    else:
        volume = -1
        # raise RuntimeError('No Volume key found. Please provide Volume (-V) of structure file (--structure).')

    return volume


def update_info(frame):
    """
    This function is used to write the info in
    the info screen.

    :param frame: the info screen.
    """

    frame.clear()
    frame.write('file name:        {}'.format(data.jfile.filename))
    frame.write('file size:        {}'.format(get_file_size(data.jfile.filename)))
    frame.write('input format:     {}'.format(data.inputformat))
    frame.write('data length:      {}'.format(data.jfile.MAX_NSTEPS))
    frame.write('selected ckeys:   {}'.format(data.jfile.select_ckeys))
    frame.write('------------------------------------')
    frame.write('Temperature:      {}'.format(data.temperature))
    frame.write('Volume:           {}'.format(data.volume))
    frame.write('DT_FS:            {}'.format(data.DT_FS))
    frame.write('psd filter width: {}'.format(data.psd_filter_width))
    frame.write('F*:               {}'.format(data.fstar))
    # if Data.xf.dct is not None:
    #     frame.write('aic type:         {}'.format(Data.xf.dct.aic_type))
    #     frame.write('aic min:          {}'.format(Data.xf.dct.aic_min))
    #     frame.write('P*:               {}'.format(Data.xf.dct.aic_Kmin + 1))
