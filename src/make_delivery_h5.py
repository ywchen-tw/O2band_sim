"""
Make a delivery copy of an assembled benchmark HDF5 file with the
solar-flux-dependent products removed.

Reflectance (I/F) is independent of the absolute solar spectrum, so the
delivery file keeps reflectance, its Monte-Carlo stderr, optical thickness,
and atmosphere/metadata, and drops ``radiance`` / ``radiance_stderr`` (which
fold in the coarse 0.5 nm solar flux).  Works on both the merged
``o2band_benchmark.h5`` (bands as groups) and the per-band ``o2a.h5`` /
``o2b.h5`` files.

Usage:
    python src/make_delivery_h5.py IN.h5 [-o OUT.h5] [--drop NAME ...]

Default output: ``IN_noradiance.h5`` next to the input.
"""

import os
import sys
import argparse

import h5py


DEFAULT_DROP = ('radiance', 'radiance_stderr')


def _copy_attrs(src, dst):
    for k, v in src.attrs.items():
        dst.attrs[k] = v


def copy_without(src, dst, drop, prefix=''):
    """Recursively copy group ``src`` into ``dst``, skipping datasets whose
    basename is in ``drop``.  Returns the list of skipped paths."""
    skipped = []
    _copy_attrs(src, dst)
    for name, obj in src.items():
        path = '%s/%s' % (prefix, name)
        if isinstance(obj, h5py.Group):
            skipped += copy_without(obj, dst.create_group(name), drop, path)
        elif name in drop:
            skipped.append(path)
        else:
            d = dst.create_dataset(name, data=obj[()],
                                   compression=obj.compression,
                                   compression_opts=obj.compression_opts)
            _copy_attrs(obj, d)
    return skipped


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    p.add_argument('input', help='assembled benchmark or per-band HDF5 file')
    p.add_argument('-o', '--output', default=None,
                   help='output path (default: <input>_noradiance.h5)')
    p.add_argument('--drop', nargs='+', default=list(DEFAULT_DROP),
                   metavar='NAME',
                   help='dataset basenames to drop (default: %(default)s)')
    args = p.parse_args(argv)

    out = args.output
    if out is None:
        stem, ext = os.path.splitext(args.input)
        out = stem + '_noradiance' + (ext or '.h5')
    if os.path.abspath(out) == os.path.abspath(args.input):
        p.error('output path equals input path')

    with h5py.File(args.input, 'r') as fi, h5py.File(out, 'w') as fo:
        skipped = copy_without(fi, fo, set(args.drop))

    for s in skipped:
        print('dropped %s' % s)
    if not skipped:
        print('warning: nothing dropped (none of %s found)' % (args.drop,))
    print('wrote %s' % out)
    return out


if __name__ == '__main__':
    main()
