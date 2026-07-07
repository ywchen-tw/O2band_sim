"""
Driver for the Phase-1 O2 A/B-band RT intercomparison benchmark.

Loops over bands x SZA x albedo, runs MCARaTS (v0.10.4, IPA) at 0.001 nm
line-by-line resolution, and writes TOA reflectance + separable O2 / H2O /
Rayleigh optical thickness with full reproducibility metadata (CLAUDE.md sec.3-4).

Design (see PLAN.md sec.10)
---------------------------
- Absorption is geometry/surface independent -> computed once per band and cached
  to disk (pickle keyed by a config hash).  The full-band 0.001 nm grid is the
  canonical lattice; an optional ``wvl_range`` selects a contiguous subset of grid
  indices (validated in-band by o2band_absorption.subrange_indices).
- Rayleigh is made exact at every 0.001 nm g-point via set_per_g_rayleigh, so
  chunking is purely operational (batch/file granularity), not a physics knob.
- Work unit = (band, SZA, albedo, chunk).  Each writes a chunk HDF5 at a
  deterministic, lattice-anchored path and is skip-if-done with atomic writes, so
  a failed run resumes by rerunning and sub-range runs compose with the full run.
- Solar source: MCARaTS runs unit-source (Src_flx=1); the CU composite solar
  spectrum F0(lambda) is folded in post to give absolute radiance, while
  reflectance rho = pi*R_raw/mu0 is F0-independent (PLAN.md sec.7.5).

Reproducibility: MCARaTS v0.10.4 is pinned explicitly via MCARATS_V010_EXE and
recorded in metadata alongside every setting, input-file identity and git commit.
"""

import os
import sys
import copy
import glob
import json
import time
import hashlib
import pickle
import subprocess
import numpy as np
import h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                   # src/  (for util)

# --- locate er3t ------------------------------------------------------------ #
# Prefer an already-importable er3t (e.g. an editable `pip install -e` in the
# active conda env -- the CURC `er3t` env has this).  Otherwise add ER3T_HOME --
# the directory that CONTAINS the `er3t` package -- to sys.path.  ER3T_HOME
# defaults to the original dev checkout for backward compatibility; on CURC set
# it (or `source setup_env.sh`) to /projects/yuch8913/wen_soft/er3t.
try:
    import er3t                                             # noqa: F401
except ImportError:
    _ER3T_HOME = os.environ.get('ER3T_HOME', '/Users/yuch8913/programming/er3t')
    if os.path.isdir(_ER3T_HOME) and _ER3T_HOME not in sys.path:
        sys.path.insert(0, _ER3T_HOME)

from util.atmosphere import afgl_atmosphere
from util.solar import solar_cu
from util.tips import tips2021
from util.absorption import hitran_lines, o2band_absorption, cal_rayleigh_od, BANDS
from util.er3t_abs import mca_atm_lbl, mca_abs_lbl, set_per_g_rayleigh
from util.mca_out_lbl import mca_out_lbl

from er3t.rtm.mca import mca_atm_1d, mcarats_ng


__all__ = ['O2BandConfig', 'O2BandSim']

_SCHEMA_VERSION = '1.0'

# CURC (Alpine) fast scratch base for simulation inputs + result files
_CURC_SCRATCH = '/scratch/alpine/yuch8913/O2band_sim'


def _default_out_dir():
    """Where simulation + result files go when ``out_dir`` is not given.

    Priority: explicit ``O2BAND_OUT_DIR`` env var > CURC/Linux scratch
    (``/scratch/alpine/yuch8913/O2band_sim``) > a local ``../out`` beside src/.
    """
    env = os.environ.get('O2BAND_OUT_DIR')
    if env:
        return env
    if sys.platform.startswith('linux'):
        return _CURC_SCRATCH
    return os.path.normpath(os.path.join(_HERE, '..', 'out'))


def _default_data_dir():
    """Where prescribed inputs (HITRAN/AFGL/solar/QTpy) live when not given.

    Priority: ``O2BAND_DATA_DIR`` env var > CURC/Linux scratch
    (``/scratch/alpine/yuch8913/O2band_sim/data``) > a local ``../data``.
    """
    env = os.environ.get('O2BAND_DATA_DIR')
    if env:
        return env
    if sys.platform.startswith('linux'):
        return os.path.join(_CURC_SCRATCH, 'data')
    return os.path.normpath(os.path.join(_HERE, '..', 'data'))


def _default_qtpy_dir(data_dir):
    """Locate the TIPS-2021 QTpy tables.

    Priority: ``O2BAND_QTPY_DIR`` env var > ``<data_dir>/TIPS_2021_PYTHON/QTpy``
    (how they are staged on CURC) > the in-repo ``src/TIPS_2021_PYTHON/QTpy`` >
    the first candidate as a best guess (tips2021 raises a clear error if absent).
    """
    env = os.environ.get('O2BAND_QTPY_DIR')
    if env:
        return env
    candidates = [os.path.join(data_dir, 'TIPS_2021_PYTHON', 'QTpy'),
                  os.path.join(_HERE, 'TIPS_2021_PYTHON', 'QTpy')]
    for c in candidates:
        if os.path.isdir(c):
            return c
    return candidates[0]


# ----------------------------------------------------------------------------- #
#  Configuration (frozen Phase-1 defaults; CLAUDE.md sec.3)
# ----------------------------------------------------------------------------- #
class O2BandConfig:

    """All settings for a run.  Defaults are the frozen Phase-1 benchmark; only
    operational knobs (paths, photons, chunking, z_top, wvl_range) are meant to
    be overridden.  Any physics deviation is recorded in metadata."""

    def __init__(self,
                 bands=('o2a', 'o2b'),
                 szas=(0.0, 30.0, 60.0),
                 albedos=(0.0, 0.1),
                 dwvl=0.001, R=20000.0, ncut=3.0,
                 include_h2o=True,
                 z_top=70.0,                       # default 70 km; production run: 120
                 photons=1.0e6, Nrun=3,
                 min_photons_per_g=1.0e4,          # floor on per-g-point photons
                 chunk_size=1000,
                 wvl_range=None,                   # (a,b) air nm, subset of a band
                 solver='ipa',
                 Ncpu='auto',
                 data_dir=None,
                 out_dir=None,
                 qtpy_dir=None,
                 mcarats_exe_env='MCARATS_V010_EXE'):

        self.bands = tuple(bands)
        self.szas = tuple(float(s) for s in szas)
        self.albedos = tuple(float(a) for a in albedos)
        self.dwvl = float(dwvl)
        self.R = float(R)
        self.ncut = float(ncut)
        self.include_h2o = bool(include_h2o)
        self.z_top = float(z_top)
        self.photons = float(photons)
        self.min_photons_per_g = float(min_photons_per_g)
        self.Nrun = int(Nrun)
        self.chunk_size = int(chunk_size)
        self.wvl_range = None if wvl_range is None else (float(wvl_range[0]), float(wvl_range[1]))
        self.solver = str(solver)
        self.Ncpu = Ncpu

        self.data_dir = data_dir or _default_data_dir()
        self.out_dir = out_dir or _default_out_dir()
        self.qtpy_dir = qtpy_dir or _default_qtpy_dir(self.data_dir)

        self.fname_hitran = os.path.join(self.data_dir, 'hitran2020_lines.txt')
        self.fname_afgl = os.path.join(self.data_dir, 'afglms.dat')
        self.fname_solar = os.path.join(self.data_dir, 'CU_composite_solar.dat')

        self.mcarats_exe_env = mcarats_exe_env
        self.mcarats_exe = os.environ.get(mcarats_exe_env, None)

    # ---- absorption cache identity: everything that changes the OD spectrum ---
    def absorption_key(self, band):
        payload = dict(band=band, dwvl=self.dwvl, R=self.R, ncut=self.ncut,
                       include_h2o=self.include_h2o, z_top=self.z_top,
                       hitran=_file_id(self.fname_hitran),
                       afgl=_file_id(self.fname_afgl))
        blob = json.dumps(payload, sort_keys=True).encode()
        return hashlib.sha1(blob).hexdigest()[:12], payload

    def metadata(self):
        return dict(
            schema_version=_SCHEMA_VERSION,
            bands=list(self.bands), szas=list(self.szas), albedos=list(self.albedos),
            dwvl_nm=self.dwvl, resolving_power_R=self.R, wing_cutoff_N=self.ncut,
            include_h2o=self.include_h2o, z_top_km=self.z_top,
            photons=self.photons, min_photons_per_g=self.min_photons_per_g,
            photons_per_g_effective=max(self.photons, self.min_photons_per_g),
            Nrun=self.Nrun, chunk_size=self.chunk_size,
            wvl_range_nm=('full' if self.wvl_range is None else list(self.wvl_range)),
            solver=self.solver,
            line_shape='Voigt', hitran_version='HITRAN 2020',
            partition_sums='TIPS-2021', wavelength_convention='air',
            o2_cia='excluded', rayleigh='Bodhaine 1999',
            reflectance_def='rho = pi*I/(mu0*F0) = pi*R_raw/mu0',
            solar_source='CU_composite_solar.dat (F0 folded post-RT)',
            mcarats_exe=self.mcarats_exe, mcarats_version='v0.10.4',
            er3t_home=os.environ.get('ER3T_HOME', ''),
            data_dir=os.path.abspath(self.data_dir),
            qtpy_dir=os.path.abspath(self.qtpy_dir),
            output_dir=os.path.abspath(self.out_dir),
            platform=sys.platform,
            git_commit=_git_commit(_HERE),
            inputs=dict(hitran=_file_id(self.fname_hitran),
                        afgl=_file_id(self.fname_afgl),
                        solar=_file_id(self.fname_solar)),
        )


# ----------------------------------------------------------------------------- #
#  Driver
# ----------------------------------------------------------------------------- #
class O2BandSim:

    def __init__(self, cfg):
        self.cfg = cfg
        if cfg.mcarats_exe is None or not os.path.isfile(cfg.mcarats_exe):
            raise OSError('Error [O2BandSim]: MCARaTS v0.10.4 executable not found '
                          'via env <%s>=%r.' % (cfg.mcarats_exe_env, cfg.mcarats_exe))
        os.makedirs(cfg.out_dir, exist_ok=True)
        self._cache_dir = os.path.join(cfg.out_dir, '_absorption_cache')
        os.makedirs(self._cache_dir, exist_ok=True)

        self.atm = afgl_atmosphere(cfg.fname_afgl, z_top=cfg.z_top)
        self.solar = solar_cu(cfg.fname_solar)
        self.tips = tips2021(cfg.qtpy_dir)
        self.mca_atm = mca_atm_lbl(self.atm)

    # --- absorption: build once per band, cache to disk ----------------------
    def absorption(self, band):
        key, payload = self.cfg.absorption_key(band)
        fcache = os.path.join(self._cache_dir, 'abs_%s_%s.pkl' % (band, key))
        if os.path.isfile(fcache):
            with open(fcache, 'rb') as f:
                return pickle.load(f)
        wl_min, wl_max = BANDS[band]
        lines = hitran_lines(self.cfg.fname_hitran, wl_range=(wl_min, wl_max), margin_cm=5.0)
        absb = o2band_absorption(self.atm, lines, band=band, dwvl=self.cfg.dwvl,
                                 R=self.cfg.R, ncut=self.cfg.ncut,
                                 include_h2o=self.cfg.include_h2o, tips=self.tips)
        _atomic_pickle(fcache, absb)
        return absb

    # --- chunk index ranges over the selected (sub)grid ----------------------
    def _chunks(self, absb):
        idx = absb.subrange_indices(self.cfg.wvl_range)     # validates in-band
        cs = self.cfg.chunk_size
        for i in range(0, idx.size, cs):
            yield idx[i:i + cs]

    def _chunk_path(self, band, sza, alb, idx_chunk):
        i0, i1 = int(idx_chunk[0]), int(idx_chunk[-1])       # inclusive, band-relative
        d = os.path.join(self.cfg.out_dir, band, 'sza%02d_alb%.2f' % (int(round(sza)), alb))
        os.makedirs(d, exist_ok=True)
        return os.path.join(d, 'chunk_%05d_%05d.h5' % (i0, i1))

    # --- one work unit: (band, sza, alb, chunk) ------------------------------
    def run_chunk(self, band, absb, idx_chunk, sza, alb, overwrite=False):
        fout = self._chunk_path(band, sza, alb, idx_chunk)
        if (not overwrite) and _chunk_valid(fout, idx_chunk.size):
            return fout, 'skip'

        abs_lbl = mca_abs_lbl(absb, idx=idx_chunk, solar=self.solar)
        atm1d = mca_atm_1d(atm_obj=self.mca_atm, abs_obj=abs_lbl)
        set_per_g_rayleigh(atm1d, self.mca_atm, abs_lbl)

        fdir = os.path.join(self.cfg.out_dir, band,
                            'sza%02d_alb%.2f' % (int(round(sza)), alb),
                            '_mca_%05d_%05d' % (int(idx_chunk[0]), int(idx_chunk[-1])))
        # Each g-point is an INDEPENDENT monochromatic solve that must get the
        # full per-g photon budget.  er3t's distribute_photon() expects weights
        # summing to 1 and splits the TOTAL photons across g-points; so pass
        # ones/Ng and scale the total by Ng -> every g-point gets p_per_g.
        # Floor the per-g count so no g-point is starved of photons (which would
        # spike MC noise or trip MCARaTS's ptot>=1 check).
        Ng = abs_lbl.Ng
        p_per_g = max(self.cfg.photons, self.cfg.min_photons_per_g)
        weights = np.full(Ng, 1.0 / Ng)
        photons_total = p_per_g * Ng
        mca0 = mcarats_ng(
            atm_1ds=[atm1d], Ng=Ng, target='radiance',
            surface_albedo=alb, solar_zenith_angle=sza, solar_azimuth_angle=0.0,
            sensor_zenith_angle=0.0, sensor_azimuth_angle=0.0,
            fdir=fdir, Nrun=self.cfg.Nrun,
            weights=weights,
            photons=photons_total, solver=self.cfg.solver,
            Ncpu=self.cfg.Ncpu, mp_mode='py', overwrite=True, quiet=True)
        out = mca_out_lbl(mca0, abs_lbl)
        if not out.valid():
            raise RuntimeError('Error [run_chunk]: non-finite radiance in %s' % fout)

        self._write_chunk(fout, band, sza, alb, idx_chunk, abs_lbl, out)
        return fout, 'run'

    def _write_chunk(self, fout, band, sza, alb, idx_chunk, abs_lbl, out):
        tmp = fout + '.tmp'
        with h5py.File(tmp, 'w') as f:
            f.attrs['band'] = band
            f.attrs['sza'] = sza
            f.attrs['albedo'] = alb
            f.attrs['idx0'] = int(idx_chunk[0])
            f.attrs['idx1'] = int(idx_chunk[-1])
            f.attrs['Ng'] = int(idx_chunk.size)
            f.attrs['photons'] = self.cfg.photons
            f.attrs['Nrun'] = self.cfg.Nrun
            f.attrs['stderr_ddof'] = 1     # std/stderr use sample std (ddof=1)
            f.attrs['mu0'] = out.mu0
            f.create_dataset('idx', data=idx_chunk.astype(np.int32))
            f.create_dataset('wvl', data=out.wvl)
            f.create_dataset('f0', data=out.f0)
            f.create_dataset('ref', data=out.ref)
            f.create_dataset('ref_std', data=out.ref_std)
            f.create_dataset('ref_stderr', data=out.ref_stderr)
            f.create_dataset('rad', data=out.rad)
            f.create_dataset('rad_std', data=out.rad_std)
            f.create_dataset('rad_stderr', data=out.rad_stderr)
            f.create_dataset('r_raw', data=out.r_raw)
        os.replace(tmp, fout)

    # --- full loop -----------------------------------------------------------
    # --- deterministic flat list of work units (for sharding) ----------------
    def _work_units(self):
        """Ordered list of (band, absb, idx_chunk, ic, nchunk, sza, alb) units.

        The order is fixed and identical in every process (band -> chunk -> sza ->
        albedo), so a stride ``global_index % ntasks == task_id`` partitions the
        units into non-overlapping shards without any inter-task coordination.
        Absorption is loaded once per band from the (cached) pickle.
        """
        units = []
        for band in self.cfg.bands:
            absb = self.absorption(band)
            chunks = list(self._chunks(absb))
            for ic, idx_chunk in enumerate(chunks):
                for sza in self.cfg.szas:
                    for alb in self.cfg.albedos:
                        units.append((band, absb, idx_chunk, ic, len(chunks), sza, alb))
        # Order by a stable digest so a plain index stride (run's shard filter)
        # is BOTH size-balanced (stride property, differs by <=1 across shards)
        # AND geometry/wavelength-interleaved (digest scrambles the nested loop
        # order, so no shard ends up all-slow-SZA).  md5, not the salted builtin
        # hash, so every array task derives the identical order.
        units.sort(key=lambda u: hashlib.md5(
            ('%s|%d|%d|%.4f|%.4f' % (u[0], int(u[2][0]), int(u[2][-1]), u[5], u[6])
             ).encode()).hexdigest())
        return units

    def prep(self, verbose=True):
        """Build + cache the absorption object for every band (no RT).  Run once
        before a job array so tasks load the cache instead of racing to write it."""
        out = []
        for band in self.cfg.bands:
            t0 = time.time()
            self.absorption(band)
            out.append(band)
            if verbose:
                print('[prep] absorption cached for %s  %.1fs' % (band, time.time() - t0))
        return out

    def run(self, overwrite=False, shard=None, verbose=True):
        """Run work units.  ``shard=(task_id, ntasks)`` runs only the units this
        task owns -- a stride over the digest-ordered unit list, so shards are
        size-balanced and geometry-interleaved (job-array mode); ``None`` runs
        all of them."""
        if shard is not None:
            task_id, ntasks = int(shard[0]), int(shard[1])
            if not (0 <= task_id < ntasks):
                raise ValueError('shard task_id %d out of range for ntasks %d'
                                 % (task_id, ntasks))
        units = self._work_units()
        done = []
        for gidx, (band, absb, idx_chunk, ic, nchunk, sza, alb) in enumerate(units):
            if shard is not None and gidx % ntasks != task_id:
                continue
            t0 = time.time()
            fout, status = self.run_chunk(band, absb, idx_chunk, sza, alb,
                                          overwrite=overwrite)
            if verbose:
                tag = '' if shard is None else '[shard %d/%d gu%d] ' % (task_id, ntasks, gidx)
                print('%s[%s sza%2.0f alb%.2f] chunk %d/%d %s  %.1fs  %s'
                      % (tag, band, sza, alb, ic + 1, nchunk, status,
                         time.time() - t0, os.path.basename(fout)))
            done.append(fout)
        return done

    # --- stitch chunk files into per-band + merged HDF5 ----------------------
    def assemble(self, verbose=True):
        merged_path = os.path.join(self.cfg.out_dir, 'o2band_benchmark.h5')
        meta = self.cfg.metadata()
        band_files = {}
        with h5py.File(merged_path, 'w') as fm:
            _write_meta(fm, meta)
            for band in self.cfg.bands:
                absb = self.absorption(band)
                idx_all = absb.subrange_indices(self.cfg.wvl_range)
                band_path = os.path.join(self.cfg.out_dir, '%s.h5' % band)
                grp = fm.create_group(band)
                self._assemble_band(band, absb, idx_all, band_path, grp, meta,
                                    verbose=verbose)
                band_files[band] = band_path
        if verbose:
            print('assembled merged file: %s' % merged_path)
        return merged_path, band_files

    def _assemble_band(self, band, absb, idx_all, band_path, merged_grp, meta, verbose=True):
        wvl = absb.wvl[idx_all]
        nsza, nalb, nw = len(self.cfg.szas), len(self.cfg.albedos), idx_all.size
        ref = np.full((nsza, nalb, nw), np.nan)
        ref_err = np.full((nsza, nalb, nw), np.nan)
        rad = np.full((nsza, nalb, nw), np.nan)
        rad_err = np.full((nsza, nalb, nw), np.nan)
        pos = {int(k): i for i, k in enumerate(idx_all)}     # band-idx -> array col

        for isza, sza in enumerate(self.cfg.szas):
            for ialb, alb in enumerate(self.cfg.albedos):
                for idx_chunk in self._chunks(absb):
                    fchunk = self._chunk_path(band, sza, alb, idx_chunk)
                    if not _chunk_valid(fchunk, idx_chunk.size):
                        raise OSError('Error [assemble]: missing/invalid chunk %s' % fchunk)
                    with h5py.File(fchunk, 'r') as f:
                        cidx = f['idx'][:]
                        cols = [pos[int(k)] for k in cidx]
                        # chunks written before the ddof fix carry population-std
                        # (ddof=0) stderr; rescale to sample-std on the fly
                        fac = 1.0
                        if int(f.attrs.get('stderr_ddof', 0)) != 1:
                            nrun = int(f.attrs.get('Nrun', self.cfg.Nrun))
                            if nrun > 1:
                                fac = float(np.sqrt(nrun / (nrun - 1.0)))
                        ref[isza, ialb, cols] = f['ref'][:]
                        ref_err[isza, ialb, cols] = f['ref_stderr'][:] * fac
                        rad[isza, ialb, cols] = f['rad'][:]
                        rad_err[isza, ialb, cols] = f['rad_stderr'][:] * fac

        # geometry-independent optical thickness on the same grid
        od_o2 = absb.od['o2'][:, idx_all]
        od_h2o = absb.od['h2o'][:, idx_all] if 'h2o' in absb.od else np.zeros_like(od_o2)
        od_ray = cal_rayleigh_od(absb.atm, wvl)               # (Nlay, Nw)
        f0 = self.solar.interp(wvl)

        # write both the standalone per-band file and the merged-file group
        for target in (band_path, merged_grp):
            _write_band_payload(target, meta, band, self.cfg, wvl, f0,
                                ref, ref_err, rad, rad_err,
                                od_o2, od_h2o, od_ray, absb)
        if verbose:
            print('assembled band %s -> %s' % (band, band_path))


# ----------------------------------------------------------------------------- #
#  helpers
# ----------------------------------------------------------------------------- #
def _write_band_payload(target, meta, band, cfg, wvl, f0,
                        ref, ref_err, rad, rad_err, od_o2, od_h2o, od_ray, absb):
    """Write a band's datasets into either a fresh HDF5 file (path) or a group."""
    close = False
    if isinstance(target, str):
        f = h5py.File(target, 'w')
        _write_meta(f, meta)
        close = True
    else:
        f = target

    f.attrs['band'] = band
    f.attrs['band_range_nm'] = list(BANDS[band])
    f.attrs['stderr_ddof'] = 1             # stderr datasets use sample std (ddof=1)
    f.create_dataset('wvl', data=wvl)                        # air nm
    f.create_dataset('sza', data=np.array(cfg.szas))
    f.create_dataset('albedo', data=np.array(cfg.albedos))
    f.create_dataset('f0', data=f0)                          # W m-2 nm-1

    d = f.create_dataset('reflectance', data=ref)
    d.attrs['dims'] = '(sza, albedo, wvl)'
    f.create_dataset('reflectance_stderr', data=ref_err).attrs['dims'] = '(sza, albedo, wvl)'
    r = f.create_dataset('radiance', data=rad)
    r.attrs['dims'] = '(sza, albedo, wvl)'; r.attrs['units'] = 'W m-2 nm-1 sr-1'
    f.create_dataset('radiance_stderr', data=rad_err).attrs['dims'] = '(sza, albedo, wvl)'

    g = f.create_group('optical_thickness')
    g.create_dataset('o2_layer', data=od_o2).attrs['dims'] = '(layer, wvl)'
    g.create_dataset('h2o_layer', data=od_h2o).attrs['dims'] = '(layer, wvl)'
    g.create_dataset('rayleigh_layer', data=od_ray).attrs['dims'] = '(layer, wvl)'
    g.create_dataset('o2_column', data=od_o2.sum(axis=0))
    g.create_dataset('h2o_column', data=od_h2o.sum(axis=0))
    g.create_dataset('rayleigh_column', data=od_ray.sum(axis=0))

    geo = f.create_group('atmosphere')
    geo.create_dataset('z_lay_km', data=absb.atm.lay['z'])
    geo.create_dataset('p_lay_hPa', data=absb.atm.lay['p'])
    geo.create_dataset('t_lay_K', data=absb.atm.lay['temperature'])
    geo.create_dataset('z_lev_km', data=absb.atm.lev['z'])
    geo.create_dataset('p_lev_hPa', data=absb.atm.lev['p'])

    if close:
        f.close()


def _write_meta(f, meta):
    g = f.create_group('metadata') if 'metadata' not in f else f['metadata']
    for k, v in meta.items():
        g.attrs[k] = json.dumps(v) if isinstance(v, (dict, list)) else v


def _chunk_valid(fout, n_expected):
    if not os.path.isfile(fout):
        return False
    try:
        with h5py.File(fout, 'r') as f:
            if f['ref'].shape[0] != n_expected:
                return False
            return bool(np.all(np.isfinite(f['ref'][:])))
    except (OSError, KeyError):
        return False


def _atomic_pickle(fout, obj):
    tmp = fout + '.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, fout)


def _file_id(fname):
    """Cheap content identity: (size, mtime). Full hash is overkill for large data."""
    try:
        st = os.stat(fname)
        return dict(path=os.path.abspath(fname), size=st.st_size,
                    mtime=int(st.st_mtime))
    except OSError:
        return dict(path=os.path.abspath(fname), size=-1, mtime=-1)


def _git_commit(path):
    try:
        return subprocess.check_output(
            ['git', '-C', path, 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return 'unknown'


def _build_cli():
    import argparse
    p = argparse.ArgumentParser(
        description='Phase-1 O2 A/B-band line-by-line RT benchmark driver.')
    p.add_argument('--bands', nargs='+', default=None,
                   help="bands to run (default frozen: o2a o2b)")
    p.add_argument('--szas', nargs='+', type=float, default=None,
                   help='solar zenith angles, deg (default 0 30 60)')
    p.add_argument('--albedos', nargs='+', type=float, default=None,
                   help='Lambertian albedos (default 0.0 0.1)')
    p.add_argument('--wvl-range', nargs=2, type=float, default=None,
                   metavar=('A', 'B'), help='air-nm sub-window (default: full band)')
    p.add_argument('--photons', type=float, default=None, help='photons per g-point')
    p.add_argument('--nrun', type=int, default=None, help='independent MC runs')
    p.add_argument('--chunk-size', type=int, default=None, help='g-points per chunk')
    p.add_argument('--z-top', type=float, default=None, help='top altitude, km')
    p.add_argument('--ncpu', default=None,
                   help="cores for MCARaTS ('auto' or an int; default auto)")
    p.add_argument('--out-dir', default=None,
                   help='output base (default: $O2BAND_OUT_DIR or scratch on Linux)')
    p.add_argument('--test', action='store_true',
                   help='quick prototype: narrow window, few photons (sanity run)')
    p.add_argument('--overwrite', action='store_true',
                   help='recompute chunks even if a valid file exists')
    p.add_argument('--no-assemble', action='store_true',
                   help='run chunks only; skip the per-band/merged HDF5 assembly')
    p.add_argument('--stage', choices=('all', 'prep', 'run', 'assemble'),
                   default='all',
                   help="pipeline stage (job-array mode): 'prep' caches absorption, "
                        "'run' executes work units (use with --shard), 'assemble' "
                        "stitches chunks; default 'all' = prep+run+assemble locally")
    p.add_argument('--shard', nargs=2, type=int, default=None,
                   metavar=('TASK_ID', 'NTASKS'),
                   help='run only work units with global_index %% NTASKS == TASK_ID '
                        '(SLURM job-array sharding); implies stage run, no assemble')
    return p


if __name__ == '__main__':
    args = _build_cli().parse_args()

    # --test presets: tiny window / small photon count for a fast end-to-end check
    kw = {}
    if args.test:
        kw.update(bands=('o2a',), szas=(30.0,), albedos=(0.0, 0.1),
                  wvl_range=(763.20, 763.30), photons=1.0e4, Nrun=2, chunk_size=50)

    # explicit CLI flags override the presets / frozen defaults
    if args.bands is not None:       kw['bands'] = tuple(args.bands)
    if args.szas is not None:        kw['szas'] = tuple(args.szas)
    if args.albedos is not None:     kw['albedos'] = tuple(args.albedos)
    if args.wvl_range is not None:   kw['wvl_range'] = tuple(args.wvl_range)
    if args.photons is not None:     kw['photons'] = args.photons
    if args.nrun is not None:        kw['Nrun'] = args.nrun
    if args.chunk_size is not None:  kw['chunk_size'] = args.chunk_size
    if args.z_top is not None:       kw['z_top'] = args.z_top
    if args.out_dir is not None:     kw['out_dir'] = args.out_dir
    if args.ncpu is not None:
        kw['Ncpu'] = args.ncpu if args.ncpu == 'auto' else int(args.ncpu)

    # --shard implies the 'run' stage of a job array (no assemble in the task)
    stage = args.stage
    shard = tuple(args.shard) if args.shard is not None else None
    if shard is not None and stage == 'all':
        stage = 'run'

    cfg = O2BandConfig(**kw)
    print('[sim_o2band] stage=%s shard=%s out_dir=%s bands=%s szas=%s albedos=%s '
          'wvl_range=%s photons=%g Nrun=%d Ncpu=%s'
          % (stage, shard, cfg.out_dir, cfg.bands, cfg.szas, cfg.albedos,
             cfg.wvl_range, cfg.photons, cfg.Nrun, cfg.Ncpu))
    sim = O2BandSim(cfg)

    if stage == 'prep':
        sim.prep()
    elif stage == 'run':
        sim.run(overwrite=args.overwrite, shard=shard)
    elif stage == 'assemble':
        sim.assemble()
    else:  # 'all' -- local end-to-end
        sim.run(overwrite=args.overwrite)
        if not args.no_assemble:
            sim.assemble()
