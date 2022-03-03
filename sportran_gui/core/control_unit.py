# -*- coding: utf-8 -*-
"""
===========================================
    Sportran graphic user interface
===========================================
--------------------------
    Control unit file
--------------------------

This file contains the functions that the GUI will directly call to start new
multithread operations and to make requests to the sportran calculus unit.
"""

import os

import sportran as st
import numpy as np

from . import settings
try:
    from sportran.utils.utils import PrintMethod
except ImportError:
    from sportran_gui.utils.utils import PrintMethod

try:
    from sportran.plotter import CurrentPlotter
except ImportError:
    raise ImportError('Couldn\'t find sportran.utils plotter.py. The GUI needs this import to work.')

# Init print method
log = PrintMethod()
info = None

# -------- DATA SECTION --------


class Jfile:

    def __init__(self, name, keys, nsteps):
        self.filename = name
        self.select_ckeys = keys
        self.MAX_NSTEPS = nsteps


class Data:
    """
    This class is used as a structure to
    contains all the variables that the software use.
    """

    loaded = False
    options = ['None', 'Energy current', 'Other current', 'Temperature (K)', 'Volume (A^3)', 'DT (fs)']

    def __init__(self):
        self._CURRENT_FILE = ''
        self.changes = True
        self.recalc_pstar = True
        self._first_fstar = True

        self.inputformat = None
        self.jfile = None
        self.jdata = None
        self.j = None
        self.xf = None
        self.axis = None

        self.keys = None
        self._description = None

        self._temperature = 0.0
        self.temperature_std = 0.0
        self._volume = 0
        self._DT_FS = 0.0
        self.currents = None

        self._fstar = 0.0
        self._pstar = 0
        self.old_fstar = 0.0
        self._psd_filter_width = 0.1
        self._units = 'metal'
        self._current_type = 'heat'

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
    def current_type(self):
        return self._current_type

    @current_type.setter
    def current_type(self, value):
        self._verify_changes(value, self.current_type)
        self._current_type = value
        print(self.changes, ' current_type')

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
        print(self.changes, ' fstar')

    @property
    def pstar(self):
        return self._pstar

    @pstar.setter
    def pstar(self, value):
        self._verify_changes(value, self._pstar, set_changes=False)
        self._pstar = value
        print(self.changes, ' pstar')

    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, value):
        self._verify_changes(self.description, value)
        self._description = value
        print(self.changes, ' descriptions')

    @property
    def CURRENT_FILE(self):
        return self._CURRENT_FILE

    @CURRENT_FILE.setter
    def CURRENT_FILE(self, value):
        self._verify_changes(self.CURRENT_FILE, value)
        if self.CURRENT_FILE != value:
            Data.loaded = False
        self._CURRENT_FILE = value
        self._first_fstar = True
        print(self.changes, ' file')

    def _verify_changes(self, val1, val2, set_changes=True):
        if val1 != val2:
            if set_changes:
                self.changes = True
            self.recalc_pstar = True

    @property
    def first_fstar(self):
        res = False
        res, self._first_fstar = self._first_fstar, res
        return res

    @first_fstar.setter
    def first_fstar(self, value):
        self._first_fstar = value


# todo: add header property
# todo: check if the behaviour is correct

try:
    data
except:
    data = Data()


def new(main):
    global data
    global log
    data = Data()
    main.show_frame(main.home)
    log = PrintMethod()


# -------- GRAPH SECTION --------
"""
This section contains the functions that deal with
the graph manager.
"""

gm = CurrentPlotter


class FakeCurrent:

    def __init__(self):
        raise RuntimeError('you must set the current type by calling select_current method')


Current = FakeCurrent


def select_current(current):
    global Current
    if current in st.current.all_currents.keys():
        Current = st.current.all_currents[current][0]
    else:
        raise KeyError(f'{current} is not a valid type of current {list(st.current.all_currents.keys())}')
    Current.set_plotter(CurrentPlotter)


def set_graph(axis_, func, **kwargs):
    """
    This function request to the graph manager to generate
    a new graph to be displayed in the axis_ screen.

    :param axis_: the location to draw the new graph.
    :param func: the function to generate the graph.
    :param kwargs: other arguments to pass to the file manager.
    :return axis: the new graph.
    """
    axis = func(axes=[axis_], **kwargs)
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
        return f'{file_size} MB'
    else:
        file_size //= 1000
        return f'{file_size} KB'


def load_settings():
    """
    This function loads from the .ini file the paths
    where the file will be stored.
    """
    global settings
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
                            settings.DATA_PATH = './'
                            #os.mkdir(settings.DATA_PATH)
                    elif var == 'OP':
                        settings.OUTPUT_PATH = val
                        if not os.path.exists(settings.OUTPUT_PATH):
                            settings.OUTPUT_PATH = './'
                            #os.mkdir(settings.OUTPUT_PATH)
                    elif var == 'LP':
                        settings.LOG_PATH = val
                        if not os.path.exists(settings.LOG_PATH):
                            settings.LOG_PATH = './'
                            #os.mkdir(settings.LOG_PATH)
                    elif var == 'FS':
                        settings.FONT_SIZE = val
                    elif var == 'PL':
                        settings.PREVIEW_LINES = val
                    elif var == 'LANG':
                        settings.LANGUAGE = val
        else:
            set_defaults()


def set_defaults():
    """
    This function is used to set the default values of the .ini file
    """
    global settings
    # Create the file
    with open('thcp.ini', 'w') as file:
        defaults = (('DP', 'data'), ('OP', 'outputs'), ('LP', 'logs'))

        # Write the defaults paths
        for setting in defaults:
            try:
                #os.mkdir(os.path.join(settings.BASE_PATH, setting[1]))
                pass
            except:
                pass
            file.write(f'{setting[0]}:{os.path.join(settings.BASE_PATH, setting[1])}\n')

        # Write default preferences
        defaults = (('FS', 11), ('PL', 10), ('LANG', 'en-EN'))
        for setting in defaults:
            file.write(f'{setting[0]}:{setting[1]}\n')

    load_settings()


def load_keys(inputfile):
    """
    This function is used to load the header keys of the selected file.

    :param inputfile: the path of the selected file.
    :return:
    """
    global data
    if data.inputformat == 'table':
        jfile = st.i_o.TableFile(inputfile, group_vectors=True)
        return jfile.all_ckeys
    elif data.inputformat == 'dict':
        try:
            data.jdata = np.load(inputfile, allow_pickle=True).tolist()
        except:   # to allow loading of python2 pickle files
            data.jdata = np.load(inputfile, allow_pickle=True, encoding='latin1').tolist()
        return {key: i for i, key in enumerate(data.jdata) if key[0] != '_'}
    elif data.inputformat == 'lammps':
        jfile = st.i_o.LAMMPSLogFile(inputfile)
        return jfile.all_ckeys
    else:
        raise RuntimeError('inputformat {} not handled'.format(data.inputformat))


def load_data(inputfile, input_format, _selected_keys, temperature=None, NSTEPS=0, START_STEP=0, run_keyword='',
              units=None, DT_FS=None, volume=None, psd_filter_w=None, axis_=None, structurefile=None, _descriptions=[]):
    global data
    selected_keys = _selected_keys.copy()
    descriptions = _descriptions.copy()
    if temperature is not None:
        data.temperature = temperature
    if volume is not None:
        data.volume = volume
    if DT_FS is not None:
        data.DT_FS = DT_FS
    data.inputformat = input_format
    data.psd_filter_width = psd_filter_w
    if units is not None:
        data.units = units

    if input_format == 'table':
        jfile = st.i_o.TableFile(inputfile, group_vectors=True)
        data.jfile = jfile
        jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
        data.jdata = jfile.data
    elif input_format == 'dict':
        data.jfile = Jfile(inputfile, _selected_keys.copy(), data.jdata[selected_keys[0]].shape[0])

        # data.jdata = np.load(inputfile) #already loaded at the header selector section
    elif input_format == 'lammps':
        jfile = st.i_o.LAMMPSLogFile(inputfile, run_keyword=run_keyword)
        jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
        data.jdata = jfile.data
    else:
        raise NotImplemented('input format not implemented.')

    ## Define currents
    # print(selected_keys, jindex)

    if descriptions.count(Data.options[3]) == 1:
        if data.inputformat != 'dict':
            temperature = get_temp(data.jdata, get_cor_index(selected_keys, descriptions, Data.options[3])[0])
            data.temperature = temperature
    # if volume is None:
    #     volume = get_volume(data.jdata, structurefile)

    # if True:    # jindex is None:
    heat_current, i = get_cor_index(selected_keys, descriptions, Data.options[1])
    currents_headers = [heat_current]
    for i, other in enumerate(selected_keys):
        if descriptions[i] == Data.options[2]:
            currents_headers.append(other)
    if NSTEPS == 0:
        NSTEPS = data.jdata[heat_current].shape[0]
    currents = np.array([data.jdata[key][START_STEP:(START_STEP + NSTEPS), :] for key in currents_headers])
    # else:
    #     pass
    #     # if sindex is None:
    #     #     currents = np.array([Data.jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] for key in selected_keys])
    #     # else:
    #     #     currents = np.array([Data.jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] - Data.jdata[key][START_STEP:(
    #     #                 START_STEP + NSTEPS), sindex] for key in selected_keys])
    data.currents = currents
    # create Current object
    emsgs = []
    if volume != -1:
        if temperature != -1:
            data.j = Current(currents, UNITS=data.units, DT_FS=DT_FS, TEMPERATURE=temperature, VOLUME=volume,
                             PSD_FILTER_W=psd_filter_w)
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


def export_data(fileout):
    global data
    if data.jdata != None:
        if Data.loaded:
            if (not 'Temperature' in data.jdata.keys() or settings.OVERWRITE) and data.temperature is not None:
                data.jdata['Temperature'] = data.temperature
            if (not 'Volume_A' in data.jdata.keys() or settings.OVERWRITE) and data.volume is not None:
                data.jdata['Volume_A'] = data.volume
            if (not 'DT_FS' in data.jdata.keys() or settings.OVERWRITE) and data.DT_FS is not None:
                data.jdata['DT_FS'] = data.DT_FS
            if (not '_UNITS' in data.jdata.keys() or settings.OVERWRITE) and data.units is not None:
                data.jdata['_UNITS'] = data.units
            if (not '_CURRENT' in data.jdata.keys() or settings.OVERWRITE) and data.current_type is not None:
                data.jdata['_CURRENT'] = data.current_type
            if (not '_HEADERS' in data.jdata.keys() or
                    settings.OVERWRITE) and data.keys is not None and data.description is not None:
                data.jdata['_HEADERS'] = {}
                data.jdata['_HEADERS']['keys'] = data.keys
                data.jdata['_HEADERS']['description'] = data.description
            if (not '_FSTAR' in data.jdata.keys() or settings.OVERWRITE) and data.fstar != 0.0:
                data.jdata['_FSTAR'] = data.fstar
            if (not '_PSTAR' in data.jdata.keys() or settings.OVERWRITE) and data.pstar != 0:
                data.jdata['_PSTAR'] = data.pstar

        np.save(fileout, data.jdata)
        return True
    return False


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
        _, volume = st.i_o.read_lammps_datafile.get_box(structurefile)
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

    print('Info')
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

    if data.xf is not None and data.xf.cepf is not None:
        frame.write('aic type:         {}'.format(data.xf.cepf.aic_type))
        frame.write('aic min:          {}'.format(data.xf.cepf.aic_min))
        frame.write('P*:               {}'.format(data.xf.cepf.cutoffK + 1))
