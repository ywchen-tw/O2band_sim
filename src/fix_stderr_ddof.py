#!/usr/bin/env python
"""
One-off, idempotent MC-stderr correction for existing output files.

Chunk and assembled HDF5s written before the ddof fix computed the run-to-run
standard deviation with the population convention (``np.std``, ddof=0) over
Nrun independent MCARaTS runs.  The unbiased sample std (ddof=1) is larger by
exactly sqrt(Nrun/(Nrun-1)) (~1.2247 at Nrun=3), so legacy std/stderr datasets
can be corrected in place by that factor -- no re-simulation, and every other
dataset stays bit-identical.

Marker + idempotency
--------------------
Every corrected chunk file (root attrs) and assembled band (root for a
per-band file, band group for a merged file) is stamped ``stderr_ddof = 1`` --
the same marker the fixed code writes -- and marked objects are skipped.  The
script is therefore safe to re-run at any time, including while still-running
array tasks (old code in memory) keep producing unmarked ddof=0 chunks: just
re-run it after they drain.  ``assemble`` also applies the identical rescale
on the fly to unmarked chunks, so assembled output is correct either way.

Writes are atomic (copy to a temp file, modify, ``os.replace``), matching the
tmp+rename convention of the driver, so a crash cannot leave a half-scaled file.

    python src/fix_stderr_ddof.py /scratch/.../z120_p1e6_n3 [more roots/files]
    python src/fix_stderr_ddof.py --dry-run OUT_DIR
"""

import os
import sys
import glob
import shutil
import argparse
import numpy as np
import h5py


CHUNK_DSETS = ['ref_std', 'ref_stderr', 'rad_std', 'rad_stderr']
BAND_DSETS = ['reflectance_stderr', 'radiance_stderr']


def _factor(nrun):
    return float(np.sqrt(nrun / (nrun - 1.0))) if nrun > 1 else 1.0


def _plan(path):
    """Inspect one HDF5; return a list of (group_path, dsets, factor) still to fix.

    Empty list -> nothing to do (already marked, or not a recognised file)."""
    todo = []
    with h5py.File(path, 'r') as f:
        if 'ref_std' in f and 'Nrun' in f.attrs:                 # chunk file
            if int(f.attrs.get('stderr_ddof', 0)) != 1:
                todo.append(('/', CHUNK_DSETS, _factor(int(f.attrs['Nrun']))))
        elif 'reflectance_stderr' in f:                          # per-band file
            if int(f.attrs.get('stderr_ddof', 0)) != 1:
                nrun = int(f['metadata'].attrs['Nrun'])
                todo.append(('/', BAND_DSETS, _factor(nrun)))
        else:                                                    # merged file?
            for name in f:
                g = f[name]
                if isinstance(g, h5py.Group) and 'reflectance_stderr' in g \
                        and int(g.attrs.get('stderr_ddof', 0)) != 1:
                    nrun = int(f['metadata'].attrs['Nrun'])
                    todo.append((name, BAND_DSETS, _factor(nrun)))
    return todo


def fix_file(path, dry_run=False):
    """Rescale + mark one file atomically.  Returns 'fixed'/'skip'/'ignored'."""
    try:
        todo = _plan(path)
    except (OSError, KeyError) as e:
        print('  [ignored] %s (%s)' % (path, e))
        return 'ignored'
    if not todo:
        return 'skip'
    if dry_run:
        for grp, dsets, fac in todo:
            print('  [dry-run] %s :: %s x%.6f (%s)' % (path, grp, fac, ','.join(dsets)))
        return 'fixed'

    tmp = path + '.fixtmp'
    shutil.copy2(path, tmp)
    try:
        with h5py.File(tmp, 'r+') as f:
            for grp, dsets, fac in todo:
                g = f if grp == '/' else f[grp]
                for d in dsets:
                    g[d][...] = g[d][...] * fac
                g.attrs['stderr_ddof'] = 1
        os.replace(tmp, path)
    except BaseException:
        if os.path.isfile(tmp):
            os.remove(tmp)
        raise
    for grp, dsets, fac in todo:
        print('  [fixed] %s :: %s x%.6f' % (path, grp, fac))
    return 'fixed'


def collect(roots):
    """Expand file/dir arguments into a sorted list of candidate HDF5 files."""
    files = []
    for root in roots:
        if os.path.isfile(root):
            files.append(root)
        elif os.path.isdir(root):
            files += glob.glob(os.path.join(root, '**', 'chunk_*.h5'), recursive=True)
            files += glob.glob(os.path.join(root, '*.h5'))
        else:
            sys.exit('Error [fix_stderr_ddof]: no such path: %s' % root)
    return sorted(set(f for f in files if not f.endswith(('.tmp', '.fixtmp'))))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('roots', nargs='+',
                   help='output tree(s) (chunk + assembled HDF5s found recursively) '
                        'and/or individual .h5 files')
    p.add_argument('--dry-run', action='store_true',
                   help='report what would be rescaled without writing')
    args = p.parse_args()

    files = collect(args.roots)
    n = {'fixed': 0, 'skip': 0, 'ignored': 0}
    for fn in files:
        n[fix_file(fn, dry_run=args.dry_run)] += 1
    print('done: %d fixed, %d already-marked/skipped, %d ignored (of %d files)'
          % (n['fixed'], n['skip'], n['ignored'], len(files)))
