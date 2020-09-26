# -*- coding: utf-8 -*-
import pytest
import numpy as np


@pytest.fixture
def check_reg(num_regression, data_regression):

    def _check(jf):
        num_regression.check({
            'psd': jf.psd,
            'logpsdK': jf.dct.logpsdK,
            'logpsdK_THEORY_std': jf.dct.logpsdK_THEORY_std,
            'logtau': jf.dct.logtau,
            'logtau_THEORY_std': jf.dct.logtau_THEORY_std,
        })
        data_regression.check({
            'KAPPA_SCALE': float(jf.KAPPA_SCALE),
            'Nyquist_f_THz': float(jf.Nyquist_f_THz),
            'kappa_Kmin': float(jf.kappa_Kmin),
            'kappa_Kmin_std': float(jf.kappa_Kmin_std),
            'aic_Kmin': int(jf.dct.aic_Kmin),
            'aic_min': float(jf.dct.aic_min)
        })

    return _check


def test_example_NaCl(data_NaCl, check_reg):
    import thermocepstrum as tc

    jfile, DT_FS, TEMPERATURE, VOLUME = data_NaCl
    j = tc.HeatCurrent([jfile.data['flux'], jfile.data['vcm[1]']], DT_FS=DT_FS, UNITS='metal', TEMPERATURE=TEMPERATURE,
                       VOLUME=VOLUME)
    print(j.Nyquist_f_THz)
    FSTAR_THZ = 14.0
    jf = j.resample(fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
    jf.cepstral_analysis()
    print(jf.cepstral_log)
    check_reg(jf)


def test_example_SiO2(data_SiO2, check_reg):
    import thermocepstrum as tc

    jfile, DT_FS, TEMPERATURE, VOLUME = data_SiO2
    j = tc.HeatCurrent(jfile.data['flux1'], DT_FS=DT_FS, UNITS='metal', TEMPERATURE=TEMPERATURE, VOLUME=VOLUME)
    print(j.Nyquist_f_THz)
    FSTAR_THZ = 28.0
    jf = j.resample(fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
    jf.cepstral_analysis()
    print(jf.cepstral_log)
    check_reg(jf)


if __name__ == '__main__':
    test_example_NaCl()
    test_example_SiO2()
