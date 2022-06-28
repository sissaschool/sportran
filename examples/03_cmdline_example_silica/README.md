# Example 3: command line (silica)

In this example we perform the same analysis that is performed in the correspondent [example notebook](../01_example_cepstrum_singlecomp_silica.ipynb), but using the command-line interface.

Simply run, after installing the package, the following command (see [`run_example.sh`](run_example.sh)):

```bash
sportran-analysis ../data/Silica/Silica.dat --input-format table -k flux1 -C heat -u metal -t 1.0 --VOLUME 3130.431110818 --param-from-input-file-column Temp TEMPERATURE -w 0.1 --FSTAR 28.0 -r
```

The options have the following meaning:

| command | explanation |
| --- | --- |
| `--input-format table` | input is a plain-text table-formatted file, column headers are used as keys |
| `-k flux1` | use the columns with header `flux1` as the main (energy) flux |
| `-C heat` | the flux is a heat current |
| `-u metal` | use LAMMPS metal units |
| `-t 1.0` | set the timestep to 1.0fs |
| `--VOLUME 3130.431110818` | set the volume of the system |
| `--param-from-input-file-column Temp TEMPERATURE` | use the column with header `Temp` as the temperature of the system |
| `-w 0.1` | set the width of the moving average filter to 0.1THz, that is used only to visualize the spectrum |
| `--FSTAR 28.0` | set the $f^*$ cutoff frequency to 28.0THz |
| `-r` | resample the time-series according to the value of $f^*$ specified with `--FSTAR` |

The [output](output_ref.log) of the program in the terminal is:

```text
 Input file (table):      ../data/Silica/Silica.dat
 Units:      metal
 Time step:      1.0 fs
# Solid Silica - BKS potential, melted and quenched
# 216 atoms, T~1000K, dens~2.295g/cm^3
# NVE, dt = 1.0 fs, 100 ps, print_step = 1.0 fs
# Temperature = 983.172635 K, Volume = 3130.431110818 A^3
# LAMMPS metal units
Temp c_flux1[1] c_flux1[2] c_flux1[3]
 #####################################
  all_ckeys = [('Temp', [0]), ('flux1', array([1, 2, 3]))]
 #####################################
Data length = 100001
  ckey = [('Temp', [0]), ('flux1', array([1, 2, 3]))]
    step =    100000 - 100.00% completed
  ( 100000 ) steps read.
DONE.  Elapsed time: 0.6583120822906494seconds
VOLUME (input): 3130.431110818
Mean TEMPERATURE (computed): 983.1726353043 +/- 39.36090003953625
 Time step (input):  1.0 fs
  currents shape is (1, 100000, 3)
snippet:
[[[ -265.30586   1520.6107      67.461829]
  [ -168.68352   1377.4459     101.82146 ]
  [  -93.688306  1180.375      117.20939 ]
  ...
  [ 1226.9778     212.0939   -1126.4643  ]
  [ 1223.8753     186.93836   -881.39541 ]
  [ 1232.7723     141.30647   -620.41895 ]]]
Using single component code.
 Number of currents = 1
 Number of equivalent components = 3
 KAPPA_SCALE = 6.144312221539281e-06
 Nyquist_f   = 500.0  THz
Using single component code.
-----------------------------------------------------
  RESAMPLE TIME SERIES
-----------------------------------------------------
 Original Nyquist freq  f_Ny =     500.00000 THz
 Resampling freq          f* =      27.77778 THz
 Sampling time         TSKIP =            18 steps
                             =        18.000 fs
 Original  n. of frequencies =         50001
 Resampled n. of frequencies =          2778
 PSD      @cutoff  (pre-filter&sample) ~ 2802468.65938
                  (post-filter&sample) ~ 2455132.46201
 log(PSD) @cutoff  (pre-filter&sample) ~     14.59530
                  (post-filter&sample) ~     14.38333
 min(PSD)          (pre-filter&sample) =      4.03008
 min(PSD)         (post-filter&sample) =  60168.84968
 % of original PSD Power f<f* (pre-filter&sample)  = 77.164 %
-----------------------------------------------------

-----------------------------------------------------
  CEPSTRAL ANALYSIS
-----------------------------------------------------
  cutoffK = (P*-1) = 33  (auto, AIC_Kmin = 33, corr_factor =  1.0)
  L_0*   =          13.114928 +/-   0.097614
  S_0*   =      717769.108506 +/- 70064.257396
-----------------------------------------------------
  kappa* =           2.205099 +/-   0.215248  W/m/K
-----------------------------------------------------
```

The program outputs raw data and some [PDF plots](output_ref.plots.pdf).

In this example the output files are called `"output.*"`.
