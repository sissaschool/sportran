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


CURRENT_FILE = ''


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
