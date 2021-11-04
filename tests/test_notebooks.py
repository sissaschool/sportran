# -*- coding: utf-8 -*-
import pytest
from testbook import testbook
import os

test_folder_path = os.path.dirname(os.path.abspath(__file__))

import functools


def local_dir(func):

    @functools.wraps(func)
    def wrapper_ld(*args, **kwargs):
        oldwd = os.getcwd()
        os.chdir(test_folder_path)
        ret = func(*args, **kwargs)
        os.chdir(oldwd)
        return ret

    return wrapper_ld


@local_dir
@testbook(test_folder_path + '/../examples/example_cepstrum_doublecomp_NaCl.ipynb', execute=True)
def test_doublecomp_NaCl(tb, file_regression):
    res = tb.cell_output_text('results_cell')
    file_regression.check(res, basename='final_result1')


@local_dir
@testbook(test_folder_path + '/../examples/example_cepstrum_singlecomp_silica.ipynb', execute=True)
def test_singlecomp_silica(tb, file_regression):
    res = tb.cell_output_text('results_cell')
    file_regression.check(res, basename='final_result2')
