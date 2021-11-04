# -*- coding: utf-8 -*-
import pytest
from testbook import testbook
import os

test_folder_path = os.path.dirname(os.path.abspath(__file__))


@testbook(test_folder_path + '/../examples/example_cepstrum_doublecomp_NaCl.ipynb', execute=True)
def test_doublecomp_NaCl(tb, file_regression):
    res = tb.cell_output_text('results_cell')
    file_regression.check(res, basename='final_result1')


@testbook(test_folder_path + '/../examples/example_cepstrum_singlecomp_silica.ipynb', execute=True)
def test_singlecomp_silica(tb, file_regression):
    res = tb.cell_output_text('results_cell')
    file_regression.check(res, basename='final_result2')
