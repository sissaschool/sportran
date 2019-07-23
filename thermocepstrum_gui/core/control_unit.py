'''
--------------------------------------------
    Thermocepstrum graphic user interface

    Control unit file
--------------------------------------------

This file contains the functions that the GUI will directly call to start new
multithread operations and to make requests to the thermocepstrum calculus unit.
'''

import os
from . import settings

import thermocepstrum as tc
import numpy as np
from thermocepstrum_gui.utils import Graph


class Data:
    CURRENT_FILE = ''
    loaded = False

    inputformat = None
    jfile = None
    jdata = None
    j = None
    xf = None
    axis = None

    keys = None
    description = None

    temperature = 0.0
    temperature_std = 0.0
    volume = 0
    DT_FS = 0.0
    currents = None

    fstar = 0.0
    old_fstar = 0.0
    psd_filter_width = 0.1


gm = Graph.GraphManager()


def set_graph(axis_, func, **kwargs):
    axis = func(axis=axis_, external_object=Data, **kwargs)
    return axis


def get_file_size(path):
    file_size = os.path.getsize(path)
    if file_size >= 1000000:
        file_size //= 1000000
        return f"{file_size} KB"
    else:
        file_size //= 1000
        return f"{file_size} MB"


def secure_exit(main_window):
    '''
    This function allow a secure exit from the execution
    of the software.
    '''

    # todo: Add multithread processes safe exit and a log file
    for process in main_window.open_windows:
        process.destroy()

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

# -------- LOAD SECTION --------

'''
In this section there are the functions that are used to load the data
that the software will use.

The starter settings are loaded from the .ini file that contains the
values that need to be stored.
'''


def load_path():
    '''
    This function loads the paths where the files are stored
    '''
    if not os.path.exists('thcp.ini'):
        set_defaults()

    with open('thcp.ini', 'r') as file:
        settings_file = file.read().split('\n')

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
    '''
    This function is used to set the default values of the .ini file
    '''
    with open('thcp.ini', 'w') as file:
        defaults = (('DP', 'data'),
                    ('OP', 'outputs'),
                    ('LP', 'logs'))

        for setting in defaults:
            try:
                os.mkdir(os.path.join(settings.BASE_PATH, setting[1]))
            except:
                pass
            file.write(f"{setting[0]}:{os.path.join(settings.BASE_PATH, setting[1])}\n")

    load_path()


def load_keys(inputfile):
    jfile = tc.i_o.TableFile(inputfile, group_vectors=True)
    return jfile.all_ckeys


def load_data(inputfile,input_format,selected_keys,temperature=None,NSTEPS=0,START_STEP=0,run_keyword='',units=None,DT_FS=None,volume=None,psd_filter_w=None,axis_=None, structurefile=None, descriptions=[]):

    Data.temperature = temperature
    Data.volume = volume
    Data.DT_FS = DT_FS
    Data.inputformat = input_format

    if input_format == 'table':
        jfile = tc.i_o.TableFile(inputfile, group_vectors=True)
        Data.jfile = jfile
        jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
        Data.jdata = jfile.data
    elif input_format == 'dict':
        Data.jdata = np.load(inputfile)
    elif input_format == 'lammps':
        jfile = tc.i_o.LAMMPSLogFile(inputfile, run_keyword=run_keyword)
        jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
        Data.jdata = jfile.data
    else:
        raise NotImplemented('input format not implemented.')

        ## Define currents
#    print(selected_keys, jindex)

    if descriptions.count('Temperature') == 1:
        temperature = get_temp(Data.jdata, get_cor_index(selected_keys, descriptions, 'Temperature')[0])
        Data.temperature = temperature
        i = get_cor_index(selected_keys, descriptions, 'Temperature')[1]
        del descriptions[i]
        del selected_keys[i]
    if volume is None:
        volume = get_volume(Data.jdata, structurefile)

    if NSTEPS == 0:
        NSTEPS = Data.jdata[list(Data.jdata.keys())[0]].shape[0]
    if True: #jindex is None:
        heat_current, i = get_cor_index(selected_keys, descriptions, 'Energy current')
        del descriptions[i]
        del selected_keys[i]
        currents_headers = [heat_current]
        for other in selected_keys:
            currents_headers.append(other)

        currents = np.array([Data.jdata[key][START_STEP:(START_STEP + NSTEPS), :] for key in currents_headers])
    else:
        pass
        # if sindex is None:
        #     currents = np.array([Data.jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] for key in selected_keys])
        # else:
        #     currents = np.array([Data.jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] - Data.jdata[key][START_STEP:(
        #                 START_STEP + NSTEPS), sindex] for key in selected_keys])
    Data.currents = currents
    # create HeatCurrent object
    emsgs = []
    if volume is not -1:
        if temperature is not -1:
            Data.j = tc.heatcurrent.HeatCurrent(currents, units, DT_FS, temperature, volume, psd_filter_w)
            gm.initialize(Data.j)
        else:
            emsgs.append('Invalid temperature!')
            return -1, emsgs
    else:
        emsgs.append('Invalid volume!')
        return -1, emsgs


def update_info(frame):
    frame.clear()
    frame.write('file name:        {}'.format(Data.jfile.filename))
    frame.write('file size:        {}'.format(get_file_size(Data.jfile.filename)))
    frame.write('input format:     {}'.format(Data.inputformat))
    frame.write('data length:      {}'.format(Data.jfile.MAX_NSTEPS))
    frame.write('selected ckeys:   {}'.format(Data.jfile.select_ckeys))
    frame.write('------------------------------------')
    frame.write('Temperature:      {}'.format(Data.temperature))
    frame.write('Volume:           {}'.format(Data.volume))
    frame.write('DT_FS:            {}'.format(Data.DT_FS))
    frame.write('psd filter width: {}'.format(Data.psd_filter_width))
    frame.write('F*:               {}'.format(Data.fstar))
    # if Data.xf.dct is not None:
    #     frame.write('aic type:         {}'.format(Data.xf.dct.aic_type))
    #     frame.write('aic min:          {}'.format(Data.xf.dct.aic_min))
    #     frame.write('P*:               {}'.format(Data.xf.dct.aic_Kmin + 1))