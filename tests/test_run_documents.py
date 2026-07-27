#!/usr/bin/env python
"""Regression checks for run-local SUMMARY.md and EVAL_REPORT.md generation."""

import os
import sys
import tempfile

import h5py
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_REPO, 'src'))

from run_documents import save_run_documents


def _fixture(path):
    with h5py.File(path, 'w') as h5:
        metadata = h5.create_group('metadata')
        metadata.attrs['wavelength_convention'] = 'vacuum'
        metadata.attrs['photons_per_g_effective'] = 1.0e7
        metadata.attrs['Nrun'] = 3
        band = h5.create_group('o2a')
        band.create_dataset('wvl', data=[757.0, 757.001, 757.002])
        band.create_dataset('sza', data=[0.0])
        band.create_dataset('albedo', data=[0.0, 0.1])
        rho = np.array([[[0.01, 0.02, 0.03], [0.10, 0.11, 0.12]]])
        band.create_dataset('reflectance', data=rho)
        band.create_dataset('reflectance_stderr', data=rho * 0.002)
        optical = band.create_group('optical_thickness')
        optical.create_dataset('o2_column', data=[0.1, 1.0, 10.0])
        optical.create_dataset('h2o_column', data=[0.0, 0.0, 0.0])
        optical.create_dataset('rayleigh_column', data=[0.03, 0.03, 0.03])


def main():
    with tempfile.TemporaryDirectory() as directory:
        h5_path = os.path.join(directory, 'o2band_benchmark.h5')
        _fixture(h5_path)

        result = save_run_documents(h5_path)
        summary = result['SUMMARY.md']['path']
        report = result['EVAL_REPORT.md']['path']
        assert os.path.dirname(summary) == directory
        assert os.path.dirname(report) == directory
        assert result['SUMMARY.md']['status'] == 'written'
        assert result['EVAL_REPORT.md']['status'] == 'written'
        assert 'vacuum' in open(summary).read()
        assert 'PASS' in open(report).read()

        # Generated files are refreshed on a later assembly.
        result = save_run_documents(h5_path)
        assert result['SUMMARY.md']['status'] == 'written'

        # A curated file with no generator marker is left untouched.
        with open(summary, 'w') as stream:
            stream.write('curated summary\n')
        result = save_run_documents(h5_path)
        assert result['SUMMARY.md']['status'] == 'preserved'
        assert open(summary).read() == 'curated summary\n'

        result = save_run_documents(h5_path, overwrite=True)
        assert result['SUMMARY.md']['status'] == 'written'
        assert 'run summary' in open(summary).read()

    print('run-document checks: PASS')


if __name__ == '__main__':
    main()
