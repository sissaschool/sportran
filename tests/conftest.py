# -*- coding: utf-8 -*-

import pytest
import os
import numpy as np

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
    return filepath_tests + '/data/NaCl/NaCl.dat'


@pytest.fixture(scope='session')
def data_SiO2_path(filepath_tests):
    return filepath_tests + '/data/Silica/Silica.dat'


@pytest.fixture(scope='session')
def data_NaCl(data_NaCl_path):
    import sportran as st
    import numpy as np

    jfile = st.i_o.TableFile(data_NaCl_path, group_vectors=True)
    jfile.read_datalines(start_step=0, NSTEPS=0, select_ckeys=['Temp', 'flux', 'vcm[1]'])
    DT_FS = 5.0   # time step [fs]
    TEMPERATURE = np.mean(jfile.data['Temp'])   # temperature [K]
    VOLUME = 40.21**3   # volume [A^3]
    print('T = {:f} K'.format(TEMPERATURE))
    print('V = {:f} A^3'.format(VOLUME))
    return jfile, DT_FS, TEMPERATURE, VOLUME


@pytest.fixture(scope='session')
def data_SiO2(data_SiO2_path):
    import sportran as st
    import numpy as np

    jfile = st.i_o.TableFile(data_SiO2_path, group_vectors=True)
    jfile.read_datalines(start_step=0, NSTEPS=0, select_ckeys=['flux1', 'Temp'])
    DT_FS = 1.0   # time step [fs]
    TEMPERATURE = np.mean(jfile.data['Temp'])   # temperature [K]
    VOLUME = 3130.431110818   # volume [A^3]
    print('T = {:f} K'.format(TEMPERATURE))
    return jfile, DT_FS, TEMPERATURE, VOLUME


@pytest.fixture
def check_reg(num_regression, data_regression):

    def _check(jf):
        num_regression.check({
            'psd': jf.psd,
            'logpsdK': jf.cepf.logpsdK,
            'logpsdK_THEORY_std': jf.cepf.logpsdK_THEORY_std,
            'logtau': jf.cepf.logtau,
            'logtau_THEORY_std': jf.cepf.logtau_THEORY_std,
            'KAPPA_SCALE': np.array([float(jf.KAPPA_SCALE)]),
            'Nyquist_f_THz': np.array([float(jf.Nyquist_f_THz)]),
            'kappa_Kmin': np.array([float(jf.kappa)]),
            'kappa_Kmin_std': np.array([float(jf.kappa_std)]),
            'aic_min': np.array([float(jf.cepf.aic_min)]),
        })
        data_regression.check({'aic_Kmin': int(jf.cepf.aic_Kmin)})

    return _check
