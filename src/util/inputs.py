"""
Locate the prescribed input files, and fail with one clear message when they
are missing.

The prescribed inputs (HITRAN 2020, the AFGL profile, the solar spectra) and the
TIPS-2021 QTpy tables are deliberately NOT in git -- ``data/`` is gitignored, so
a fresh clone has no ``data/`` at all.  On CURC they live on scratch and
``setup_env.sh`` points ``O2BAND_DATA_DIR`` / ``O2BAND_QTPY_DIR`` at them.

Forgetting to source it is the common failure, and without a guard it surfaces
as a stack of FileNotFoundError tracebacks from deep inside numpy/pickle, one
per check, which buries the actual cause.  ``require_inputs`` turns that into a
single line naming the missing files, the directory searched, and the fix.

Note ``util.tips`` reads ``O2BAND_QTPY_DIR`` at *import* time, so the variable
must be exported before Python starts -- setting ``os.environ`` inside a running
session is too late.  That is why the remedy is "source setup_env.sh", not "set
the variable".
"""

import os

__all__ = ['data_dir', 'qtpy_dir', 'require_inputs', 'MissingInputs']

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, '..', '..'))


class MissingInputs(OSError):
    """Prescribed input files could not be found."""


def data_dir():
    """Directory holding the prescribed inputs: $O2BAND_DATA_DIR, else <repo>/data."""
    return os.environ.get('O2BAND_DATA_DIR', os.path.join(_REPO, 'data'))


def qtpy_dir():
    """TIPS-2021 QTpy directory: $O2BAND_QTPY_DIR, else the in-repo copy."""
    return os.environ.get(
        'O2BAND_QTPY_DIR',
        os.path.join(_REPO, 'src', 'TIPS_2021_PYTHON', 'QTpy'))


def require_inputs(*names, **kwargs):
    """
    Check that every file in ``names`` exists under :func:`data_dir`, and
    optionally that the QTpy directory exists (``qtpy=True``).

    Raises :class:`MissingInputs` naming what is missing and how to fix it.
    Returns the data directory so callers can use it directly.
    """
    qtpy = kwargs.pop('qtpy', False)
    if kwargs:
        raise TypeError('unexpected keyword arguments: %s' % sorted(kwargs))

    ddir = data_dir()
    missing = [n for n in names if not os.path.isfile(os.path.join(ddir, n))]
    qdir = qtpy_dir()
    qtpy_missing = qtpy and not os.path.isdir(qdir)

    if missing or qtpy_missing:
        lines = ['prescribed input(s) not found:']
        for n in missing:
            lines.append('    %s' % os.path.join(ddir, n))
        if qtpy_missing:
            lines.append('    %s   (TIPS-2021 QTpy tables)' % qdir)
        lines.append('')
        lines.append('  data dir searched: %s  (from %s)'
                     % (ddir, 'O2BAND_DATA_DIR'
                        if os.environ.get('O2BAND_DATA_DIR') else 'default <repo>/data'))
        lines.append('  data/ is gitignored, so a fresh clone has none; on CURC the')
        lines.append('  inputs live on scratch.  Fix:')
        lines.append('      source setup_env.sh        # exports O2BAND_DATA_DIR / _QTPY_DIR')
        lines.append('  (must be sourced BEFORE python starts -- util.tips reads')
        lines.append('   O2BAND_QTPY_DIR at import time)')
        raise MissingInputs('\n'.join(lines))

    return ddir
