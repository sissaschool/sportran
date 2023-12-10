# -*- coding: utf-8 -*-
import sportran as st
import numpy as np
import yaml, pandas as pd
import argparse


def write_results(jf, filename):
    d = {
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
    }
    df = pd.DataFrame({k: pd.Series(v) for k, v in d.items()})
    df.to_csv(filename + '.csv')

    d = {'aic_Kmin': int(jf.cepf.aic_Kmin)}
    with open(filename + '.yml', 'w') as f:
        f.write(yaml.dump(d))

def main(inputs):
    if (inputs == 'all' or inputs == 'NaCl' ):
        # NaCl
        jfile = st.i_o.TableFile('../data/NaCl/NaCl.dat', group_vectors=True)
        jfile.read_datalines(start_step=0, NSTEPS=0, select_ckeys=['Temp', 'flux', 'vcm[1]'])
        DT_FS = 5.0   # time step [fs]
        TEMPERATURE = np.mean(jfile.data['Temp'])   # temperature [K]
        VOLUME = 40.21**3   # volume [A^3]
        print('T = {:f} K'.format(TEMPERATURE))
        print('V = {:f} A^3'.format(VOLUME))

        j = st.HeatCurrent([jfile.data['flux'], jfile.data['vcm[1]']], DT_FS=DT_FS, UNITS='metal', TEMPERATURE=TEMPERATURE,
                        VOLUME=VOLUME)
        print(j.Nyquist_f_THz)
        FSTAR_THZ = 14.0
        jf = j.resample(fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
        jf.cepstral_analysis()
        print(jf.cepstral_log)

        write_results(jf, 'test_example_NaCl')
    
    if (inputs == 'all' or inputs == 'SiO2' ):
        # SiO2
        jfile = st.i_o.TableFile('../data/Silica/Silica.dat', group_vectors=True)
        jfile.read_datalines(start_step=0, NSTEPS=0, select_ckeys=['flux1', 'Temp'])
        DT_FS = 1.0   # time step [fs]
        TEMPERATURE = np.mean(jfile.data['Temp'])   # temperature [K]
        VOLUME = 3130.431110818   # volume [A^3]
        print('T = {:f} K'.format(TEMPERATURE))

        j = st.HeatCurrent(jfile.data['flux1'], DT_FS=DT_FS, UNITS='metal', TEMPERATURE=TEMPERATURE, VOLUME=VOLUME)
        print(j.Nyquist_f_THz)
        FSTAR_THZ = 28.0
        jf = j.resample(fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
        jf.cepstral_analysis()
        print(jf.cepstral_log)

        write_results(jf, 'test_example_SiO2')

        j = st.HeatCurrent(jfile.data['flux1'], DT_FS=DT_FS, UNITS='metal', TEMPERATURE=TEMPERATURE, VOLUME=VOLUME)
        print(j.Nyquist_f_THz)
        FSTAR_THZ = 28.0
        jf = j.resample(fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
        jf.cepstral_analysis(manual_cutoffK=42)
        print(jf.cepstral_log)

        write_results(jf, 'test_example_SiO2_fixed_K')
        
    if (inputs == 'all' or inputs == 'h2o' ):
        # Viscosity H2O
        jfile = st.i_o.TableFile('../data/h2o/h2o.dat', group_vectors=True)
        jfile.read_datalines()
        DT_FS = 1                                  # time step [fs]
        TEMPERATURE = np.mean(jfile.data['Temp'])  # temperature [K] (983.173 K)
        VOLUME = 1918.4149658                      # volume [A^3]
        print('T = {:f} K'.format(TEMPERATURE))
        print('V = {:f} A^3'.format(VOLUME))

        j = st.StressCurrent(np.column_stack([jfile.data['s_xy'], jfile.data['s_xz'], jfile.data['s_yz']]),
                            UNITS='real', DT_FS=DT_FS, TEMPERATURE=TEMPERATURE, VOLUME=VOLUME)
        print(j.Nyquist_f_THz)
        FSTAR_THZ = 25.0
        jf = j.resample(fstar_THz=FSTAR_THZ, plot=False, freq_units='thz')
        jf.cepstral_analysis()
        print(jf.cepstral_log)

        write_results(jf, 'test_example_viscosity_h2o')


if  __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Regenerate data for the test.')
    parser.add_argument('-i', type=str,
                    help='''what data do you want to regenerate? 
                    Valid options "all", "NaCl", "SiO2", "h2o".
                    Combinations can be done adding "_", e.g. "NaCl_SiO2"''')

    inputs = parser.parse_args().i.split('_')
    for i in inputs:
        if i not in ["all", "NaCl", "SiO2", "h2o"]:
            print('input values not valid use "-h" to get help')
        else:
            main(i)
