def test_example_NaCl():
    import numpy as np
    import os
    import thermocepstrum as tc

    print(os.environ)

    jfile = tc.i_o.TableFile('./data/NaCl.dat', group_vectors=True)
    jfile.read_datalines(start_step=0, NSTEPS=0, select_ckeys=['Temp', 'flux', 'vcm[1]'])
    DT_FS = 5.0   # time step [fs]
    TEMPERATURE = np.mean(jfile.data['Temp'])   # temperature [K]
    VOLUME = 40.21**3   # volume [A^3]
    print('T = {:f} K'.format(TEMPERATURE))
    print('V = {:f} A^3'.format(VOLUME))

    j = tc.HeatCurrent([jfile.data['flux'], jfile.data['vcm[1]']], 'metal', DT_FS, TEMPERATURE, VOLUME)
    print(j.Nyquist_f_THz)

    FSTAR_THZ = 14.0
    jf = tc.heatcurrent.resample_current(j, fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
    jf.cepstral_analysis()
    print('K of AIC_min = {:d}'.format(jf.dct.aic_Kmin))
    print('AIC_min = {:f}'.format(jf.dct.aic_min))
    print(jf.cepstral_log)

    assert abs(jf.dct.aic_min - 1410.555757) < 1.0e-6
    assert jf.dct.aic_Kmin == 3
    assert abs(jf.kappa_Kmin - 0.498310) < 1.0e-6
    assert abs(jf.kappa_Kmin_std - 0.028018) < 1.0e-6
    print('*********************\n   TEST:  passed.\n*********************\n')


def test_example_SiO2():

    import os
    import thermocepstrum as tc

    print(os.environ)

    jfile = tc.i_o.TableFile('./data/Silica.dat', group_vectors=True)
    jfile.read_datalines(start_step=0, NSTEPS=0, select_ckeys=['flux1'])
    DT_FS = 1.0   # time step [fs]
    TEMPERATURE = 1065.705630   # temperature [K]
    VOLUME = 3130.431110818   # volume [A^3]

    j = tc.HeatCurrent(jfile.data['flux1'], 'metal', DT_FS, TEMPERATURE, VOLUME)
    print(j.Nyquist_f_THz)

    FSTAR_THZ = 28.0
    jf = tc.heatcurrent.resample_current(j, fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
    jf.cepstral_analysis()
    print('K of AIC_min = {:d}'.format(jf.dct.aic_Kmin))
    print('AIC_min = {:f}'.format(jf.dct.aic_min))
    print(jf.cepstral_log)

    assert abs(jf.dct.aic_min - 2820.804281) < 1.0e-6
    assert jf.dct.aic_Kmin == 37
    assert abs(jf.kappa_Kmin - 2.484534) < 1.0e-6
    assert abs(jf.kappa_Kmin_std - 0.256596) < 1.0e-6
    print('*********************\n   TEST:  passed.\n*********************\n')


if __name__ == '__main__':
    test_example_NaCl()
    test_example_SiO2()
