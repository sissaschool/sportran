import pytest
import os

pytest_plugins = 'pytester'


@pytest.fixture(scope='session')
def filepath_tests():
    """Return the absolute filepath of the `tests` folder.

    .. warning:: if this file moves with respect to the `tests` folder, the implementation should change.

    :return: absolute filepath of `tests` folder which is the basepath for all test resources.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope='session')
def data_NaCl_path(filepath_tests):
    return filepath_tests + '/data/NaCl.dat'


@pytest.fixture(scope='session')
def data_SiO2_path(filepath_tests):
    return filepath_tests + '/data/Silica.dat'


@pytest.fixture(scope='session')
def data_NaCl(data_NaCl_path):
    import thermocepstrum as tc
    import numpy as np

    jfile = tc.i_o.TableFile(data_NaCl_path, group_vectors=True)
    jfile.read_datalines(start_step=0, NSTEPS=0, select_ckeys=['Temp', 'flux', 'vcm[1]'])
    DT_FS = 5.0   # time step [fs]
    TEMPERATURE = np.mean(jfile.data['Temp'])   # temperature [K]
    VOLUME = 40.21**3   # volume [A^3]
    print('T = {:f} K'.format(TEMPERATURE))
    print('V = {:f} A^3'.format(VOLUME))
    return jfile, DT_FS, TEMPERATURE, VOLUME


@pytest.fixture(scope='session')
def data_SiO2(data_SiO2_path):
    import thermocepstrum as tc

    jfile = tc.i_o.TableFile(data_SiO2_path, group_vectors=True)
    jfile.read_datalines(start_step=0, NSTEPS=0, select_ckeys=['flux1'])
    DT_FS = 1.0   # time step [fs]
    TEMPERATURE = 1065.705630   # temperature [K]
    VOLUME = 3130.431110818   # volume [A^3]
    return jfile, DT_FS, TEMPERATURE, VOLUME
