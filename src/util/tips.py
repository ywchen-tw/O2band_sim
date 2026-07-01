"""
HITRAN TIPS-2021 total internal partition sums (Gamache et al. 2021), used for
the line-intensity temperature scaling Q(296)/Q(T).

TIPS-2021 is the partition-sum set produced for the HITRAN2020 line list, so it
is internally consistent with the S(296 K) intensities in that database (the
S(296) values were normalised with these same Q(296)).  Using it therefore
avoids a Q(296) mismatch in the temperature correction.

Data source: ``src/TIPS_2021_PYTHON/QTpy/{mol}_{iso}.QTpy`` -- one pickle per
(HITRAN molecule id, isotopologue), a dict mapping integer temperature (K, as a
string) to Q(T).  We read Q(296) exactly and linearly interpolate Q(T).
"""

import os
import pickle
import functools
import numpy as np


__all__ = ['tips2021']

# default location of the TIPS-2021 QTpy tables (src/TIPS_2021_PYTHON/QTpy)
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_QTPY_DIR = os.path.normpath(os.path.join(_HERE, '..', 'TIPS_2021_PYTHON', 'QTpy'))

T_REF = 296.0  # HITRAN reference temperature (K)


class tips2021:

    """
    Q(296)/Q(T) from TIPS-2021.

    Input:
        qtpy_dir : directory holding the {mol}_{iso}.QTpy pickle files
                   (default: src/TIPS_2021_PYTHON/QTpy)

    Methods:
        Q(mol, iso, T)      : partition sum Q at temperature T (K), scalar or array
        ratio(mol, iso, T)  : Q(296)/Q(T)  (the factor used in S(T) scaling)
    """

    def __init__(self, qtpy_dir=DEFAULT_QTPY_DIR):
        self.qtpy_dir = qtpy_dir
        if not os.path.isdir(qtpy_dir):
            raise OSError('Error [tips2021]: QTpy directory not found: %s' % qtpy_dir)

    @functools.lru_cache(maxsize=None)
    def _table(self, mol, iso):
        fname = os.path.join(self.qtpy_dir, '%d_%d.QTpy' % (mol, iso))
        if not os.path.isfile(fname):
            raise OSError('Error [tips2021]: missing partition file %s' % fname)
        with open(fname, 'rb') as f:
            d = pickle.loads(f.read())
        # convert {str(T): Q} -> sorted arrays for interpolation
        temps = np.array(sorted(int(k) for k in d), dtype=np.float64)
        qvals = np.array([float(d[str(int(t))]) for t in temps], dtype=np.float64)
        return temps, qvals

    def Q(self, mol, iso, T):
        temps, qvals = self._table(mol, iso)
        return np.interp(T, temps, qvals)

    def ratio(self, mol, iso, T):
        """Q(296)/Q(T) -- the partition factor in the HITRAN intensity scaling."""
        q296 = self.Q(mol, iso, T_REF)
        return q296 / self.Q(mol, iso, T)


if __name__ == '__main__':

    tips = tips2021()
    print('QTpy dir:', tips.qtpy_dir)
    for mol, iso, name in [(7, 1, 'O2 (66)'), (1, 1, 'H2O (161)')]:
        print('%-10s Q(296)=%.4f  Q(250)=%.4f  Q(296)/Q(250)=%.5f'
              % (name, tips.Q(mol, iso, 296.0), tips.Q(mol, iso, 250.0),
                 tips.ratio(mol, iso, 250.0)))
