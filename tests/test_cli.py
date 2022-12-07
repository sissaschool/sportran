# -*- coding: utf-8 -*-

import pytest
import numpy as np


@pytest.fixture
def run_cli(testdir, filepath_tests):

    def do_run(*args):
        args = [filepath_tests + '/../sportran/analysis.py'] + list(args)
        return testdir.run(*args)

    return do_run


def test_cli_NaCl(tmpdir, run_cli, data_NaCl_path, num_regression, file_regression):
    fileout_pref = tmpdir.join('output')
    output = run_cli(data_NaCl_path, '-k', 'flux', '-j', 'vcm[1]', '-t', '5.0', '--VOLUME', '65013.301261',
                     '--param-from-input-file-column', 'Temp', 'TEMPERATURE', '-w', '0.1', '--FSTAR', '14.0', '-r',
                     '--test-suite-run', '-o', fileout_pref)
    file_list = [
        'output.cepstral.dat', 'output.cepstrumfiltered_psd.dat', 'output.cospectrum.dat', 'output.cospectrum.filt.dat',
        'output.psd.dat', 'output.resampled_psd.dat'
    ]
    with open(str(tmpdir.join('output.log'))) as l:
        file_regression.check(l.read(), basename='output.log')
    file_regression.check(output.stdout.str(), basename='stdout')
    file_regression.check(output.stderr.str(), basename='stderr')
    readed = {}
    for f in file_list:
        file_out = tmpdir.join(f)
        readed[f] = np.loadtxt(file_out).flatten()

    num_regression.check(readed)


def test_cli_NaCl_no_w(tmpdir, run_cli, data_NaCl_path, num_regression, file_regression):
    fileout_pref = tmpdir.join('output_no_w')
    output = run_cli(data_NaCl_path, '-k', 'flux', '-j', 'vcm[1]', '-t', '5.0', '--VOLUME', '65013.301261',
                     '--param-from-input-file-column', 'Temp', 'TEMPERATURE', '--FSTAR', '14.0', '-r',
                     '--test-suite-run', '-o', fileout_pref)
    file_list = [
        'output_no_w.cepstral.dat', 'output_no_w.cepstrumfiltered_psd.dat', 'output_no_w.cospectrum.dat',
        'output_no_w.psd.dat', 'output_no_w.resampled_psd.dat'
    ]
    with open(str(tmpdir.join('output_no_w.log'))) as l:
        file_regression.check(l.read(), basename='output_no_w.log')
    file_regression.check(output.stdout.str(), basename='stdout_no_w')
    file_regression.check(output.stderr.str(), basename='stderr_no_w')
    readed = {}
    for f in file_list:
        file_out = tmpdir.join(f)
        readed[f] = np.loadtxt(file_out).flatten()

    num_regression.check(readed)


def test_cli_bin_output_NaCl(tmpdir, run_cli, data_NaCl_path, num_regression, file_regression):
    fileout_pref = tmpdir.join('output')
    output = run_cli(data_NaCl_path, '-k', 'flux', '-j', 'vcm[1]', '-t', '5.0', '--VOLUME', '65013.301261',
                     '--param-from-input-file-column', 'Temp', 'TEMPERATURE', '-w', '0.1', '--FSTAR', '14.0', '-r',
                     '--test-suite-run', '-O', '--bin-output-old', '-o', fileout_pref)
    file_list = [
        'output.cepstral.npy', 'output.cepstrumfiltered_psd.npy', 'output.cospectrum.npy', 'output.cospectrum.filt.npy',
        'output.psd.npy', 'output.resampled_psd.npy'
    ]
    with open(str(tmpdir.join('output.log'))) as l:
        file_regression.check(l.read(), basename='output.log')
    file_regression.check(output.stdout.str(), basename='stdout')
    file_regression.check(output.stderr.str(), basename='stderr')
    readed = {}
    for f in file_list:
        file_out = tmpdir.join(f)
        readed[f] = np.load(str(file_out), allow_pickle=True).flatten().real
        # NOTE: TODO: imaginary part of cospectrum is not checked!
    num_regression.check(readed)
