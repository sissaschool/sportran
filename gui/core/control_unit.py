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
from utils import Graph


class Data:
    CURRENT_FILE = ''
    loaded = False

    jdata = None
    j = None
    axis = None

    fstar = 0.0

gm = Graph.GraphManager()


def set_graph(axis_, func, **kwargs):
    gm.initialize(Data.j)
    axis = func(axis=axis_, **kwargs)
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


def load_data(inputfile,input_format,selected_keys,temperature=None,NSTEPS=0,START_STEP=0,run_keyword='',units=None,DT_FS=None,volume=None,psd_filter_w=None,axis_=None, logs=None):

    if input_format == 'table':
        if temperature is None:
            selected_keys.append('Temp')
        #      if 'Press' in jfile.ckey:
        #         selected_keys.append('Press')
        jfile = tc.i_o.TableFile(inputfile, group_vectors=True)
        jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
        Data.jdata = jfile.data
    elif input_format == 'dict':
        Data.jdata = np.load(inputfile)
    elif input_format == 'lammps':
        jfile = tc.i_o.LAMMPSLogFile(inputfile, run_keyword=run_keyword)
        if temperature is None:
            selected_keys.append('Temp')
        #      if 'Press' in jfile.ckey:
        #         selected_keys.append('Press')
        jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
        Data.jdata = jfile.data
    else:
        raise NotImplemented('input format not implemented.')

        ## Define currents
#    print(selected_keys, jindex)

    if NSTEPS == 0:
        NSTEPS = Data.jdata[list(Data.jdata.keys())[0]].shape[0]
    if True: #jindex is None:
        currents = np.array([Data.jdata[key][START_STEP:(START_STEP + NSTEPS), :] for key in selected_keys])
    else:
        pass
        # if sindex is None:
        #     currents = np.array([Data.jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] for key in selected_keys])
        # else:
        #     currents = np.array([Data.jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] - Data.jdata[key][START_STEP:(
        #                 START_STEP + NSTEPS), sindex] for key in selected_keys])

    if logs:
        logs.write('currents shape is {}'.format(currents.shape))
        logs.write('snippets')
        logs.write(currents)
    else:
        print('  currents shape is {}'.format(currents.shape))
        # logfile.write('  currents shape is {}\n'.format(currents.shape))
        print('snippet:')
        print(currents)

    # create HeatCurrent object
    Data.j = tc.heatcurrent.HeatCurrent(currents, units, DT_FS, temperature, volume, psd_filter_w)