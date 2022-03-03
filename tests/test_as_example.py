# -*- coding: utf-8 -*-
import pytest
import numpy as np


def test_example_NaCl(data_NaCl, check_reg):
    import sportran as st

    jfile, DT_FS, TEMPERATURE, VOLUME = data_NaCl
    j = st.HeatCurrent([jfile.data['flux'], jfile.data['vcm[1]']], DT_FS=DT_FS, UNITS='metal', TEMPERATURE=TEMPERATURE,
                       VOLUME=VOLUME)
    print(j.Nyquist_f_THz)
    FSTAR_THZ = 14.0
    jf = j.resample(fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
    jf.cepstral_analysis()
    print(jf.cepstral_log)
    check_reg(jf)


def test_example_SiO2(data_SiO2, check_reg):
    import sportran as st

    jfile, DT_FS, TEMPERATURE, VOLUME = data_SiO2
    j = st.HeatCurrent(jfile.data['flux1'], DT_FS=DT_FS, UNITS='metal', TEMPERATURE=TEMPERATURE, VOLUME=VOLUME)
    print(j.Nyquist_f_THz)
    FSTAR_THZ = 28.0
    jf = j.resample(fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
    jf.cepstral_analysis()
    print(jf.cepstral_log)
    check_reg(jf)


def test_example_SiO2_fixed_K(data_SiO2, check_reg):
    import sportran as st

    jfile, DT_FS, TEMPERATURE, VOLUME = data_SiO2
    j = st.HeatCurrent(jfile.data['flux1'], DT_FS=DT_FS, UNITS='metal', TEMPERATURE=TEMPERATURE, VOLUME=VOLUME)
    print(j.Nyquist_f_THz)
    FSTAR_THZ = 28.0
    jf = j.resample(fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
    jf.cepstral_analysis(manual_cutoffK=42)
    print(jf.cepstral_log)
    check_reg(jf)


if __name__ == '__main__':
    test_example_NaCl()
    test_example_SiO2()
