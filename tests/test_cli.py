import pytest


@pytest.fixture
def run_cli(testdir, filepath_tests):

    def do_run(*args):
        args = [filepath_tests + '/../thermocepstrum/analysis.py'] + list(args)
        return testdir.run(*args)

    return do_run


def test_cli_NaCl(tmpdir, run_cli, data_NaCl_path, file_regression):
    #FIXME: test in a branch where cli works
    fileout_pref = tmpdir.join('output')
    output = run_cli(data_NaCl_path, '-k', 'flux', '-j', 'vcm[1]', '-t', '5.0', '-V', '65013.301261', '-w', '0.1',
                     '--FSTAR', '14.0', '-r', '--test-suite-run', '--savetxt-format', '%.14e', '-o', fileout_pref)
    file_list = [
        'output.cepstral.dat', 'output.cepstrumfiltered_psd.dat', 'output.cospectrum.dat', 'output.cospectrum.filt.dat',
        'output.log', 'output.psd.dat', 'output.resampled_psd.dat'
    ]
    file_regression.check(output.stdout.str(), basename='stdout')
    file_regression.check(output.stderr.str(), basename='stderr')
    for f in file_list:
        file_out = tmpdir.join(f)
        with file_out.open('r') as outp:
            print(f)
            file_regression.check(outp.read(), basename=f)
