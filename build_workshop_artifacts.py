"""Pre-build the artifacts consumed by ``Euclid_AGN_Workshop.ipynb``.

This script runs offline before the AAS workshop and produces a small,
self-contained set of files so the 30-minute hands-on notebook does no
SOM training and no batched IRSA queries at workshop time.

What it builds
--------------
- ``data/cache/workshop_som.pkl``           : trained MiniSom + (mu, sigma, qe_thresh)
- ``data/cache/workshop_agn_edfn.fits``     : EDF-N AGN that pass photometry gates,
                                              with bmu_row, bmu_col, qe, colors, mH, z attached
- ``data/cache/workshop_train_edfn.fits``   : north-system training galaxies with bmu attached
- ``data/cache/workshop_picks.json``        : the three pre-picked AGN with notes

It reads the existing tutorial caches (``tutorial_*.fits``, ``tutorial_*.ecsv``,
``irsa_spectra_agn/``) and the EDF-N AGN catalog, so it does not call IRSA.
"""

from __future__ import annotations

import json
import os
import pickle
import sys

import numpy as np
from astropy.table import Table, join
from minisom import MiniSom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from agn_tutorial_utils import (
    COLOR_NAMES, MSZ, abmag, as_float_array, coalesce, project, snr, valid_mag,
)

ROOT = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(ROOT, 'data', 'input')
CACHE_DIR = os.path.join(ROOT, 'data', 'cache')
SPECTRA_CACHE_DIR = os.path.join(CACHE_DIR, 'irsa_spectra_agn')

Z_WIN = (1.30, 1.80)
SNR_MIN = 5.0
MAG_MIN, MAG_MAX = 12, 30
N_ITER = 100_000
RANDOM_SEED = 0
QE_PCTILE = 95

COMMON_FLAGS = [
    'PRF_qso_candidate', 'B24a_qso_candidate',
    'C75_agn_candidate', 'R90_agn_candidate',
    'GDR3_qso_candidate',
    'JH_IeY_qso_candidate', 'IeH_gz_qso_candidate',
]
EDFN_EXTRA = [
    'B24b_qso_candidate',
    'DESI_qso_candidate', 'DESI_broadline_galaxy_candidate', 'DESI_broadline_qso_candidate',
    'DESI_niibpt_agn_candidate', 'DESI_siibpt_agn_candidate', 'DESI_oibpt_agn_candidate',
    'DESI_whan_agn_candidate', 'DESI_blue_agn_candidate', 'DESI_kex_agn_candidate',
    'AGN_sed_candidate',
]
ALL_FLAGS = COMMON_FLAGS + EDFN_EXTRA


def load_edfn_agn():
    t = Table.read(os.path.join(INPUT_DIR, 'edfn.fits'))
    bright = (t['bright_vis_mag_bin'] == 1) | (t['medium_vis_mag_bin'] == 1) | (t['faint_vis_mag_bin'] == 1)
    clean = (t['good_flags'] == 1) & bright & (t['star_candidate_all'] == 0)
    t = t[clean]

    any_agn = np.zeros(len(t), dtype=bool)
    for f in ALL_FLAGS:
        any_agn |= (t[f] == 1)
    t = t[any_agn]

    keep = ['object_id_euclid', 'ra_euclid', 'dec_euclid', 'Z_desi', 'Z_sdss'] + ALL_FLAGS
    out = t[keep]
    print(f'EDF-N AGN candidates after catalog gates: {len(out):,}')
    return out


def attach_z_and_phot(agn):
    """Join the EDF-N AGN with the cached photo-z table and MER photometry,
    then derive z_final, z_source, six adjacent colors, and the S/N gate.
    """
    z = Table.read(os.path.join(CACHE_DIR, 'tutorial_agn_z.fits'))
    if 'object_id' in z.colnames:
        z.rename_column('object_id', 'object_id_euclid')
    agn_z = join(agn, z, keys='object_id_euclid', join_type='left')

    z_spec = np.array(agn_z['Z_desi'], dtype=float)
    bad = ~np.isfinite(z_spec) | (z_spec <= 0)
    z_spec[bad] = np.array(agn_z['Z_sdss'], dtype=float)[bad]
    z_phot = np.array(agn_z['phz_median'], dtype=float)
    agn_z['z_final']  = np.where(np.isfinite(z_spec) & (z_spec > 0), z_spec, z_phot)
    agn_z['z_source'] = np.where(np.isfinite(z_spec) & (z_spec > 0), 'spec', 'phot')

    in_z = np.isfinite(agn_z['z_final']) & (agn_z['z_final'] > Z_WIN[0]) & (agn_z['z_final'] < Z_WIN[1])
    agn_in_z = agn_z[in_z]
    n_spec = int(np.sum(agn_in_z['z_source'] == 'spec'))
    print(f'In z window {Z_WIN}: {len(agn_in_z):,} AGN ({n_spec:,} with spec-z)')

    phot = Table.read(os.path.join(CACHE_DIR, 'tutorial_agn_phot.fits'))
    if 'object_id' in phot.colnames:
        phot.rename_column('object_id', 'object_id_euclid')
    merged = join(agn_in_z, phot, keys='object_id_euclid', join_type='inner')

    # EDF-N native: HSC g, MegaCam r, PanSTARRS i, HSC z. coalesce() picks the
    # north-system fluxes (which exist) and tags survival as 'north'.
    gf, ge, gsys = coalesce(merged, 'flux_g_ext_hsc_templfit',       'flux_g_ext_decam_templfit')
    rf, re_, _   = coalesce(merged, 'flux_r_ext_megacam_templfit',   'flux_r_ext_decam_templfit')
    if_, ie, _   = coalesce(merged, 'flux_i_ext_panstarrs_templfit', 'flux_i_ext_decam_templfit')
    zf, ze, _    = coalesce(merged, 'flux_z_ext_hsc_templfit',       'flux_z_ext_decam_templfit')
    yf = np.asarray(merged['flux_y_templfit'], float); ye = np.asarray(merged['fluxerr_y_templfit'], float)
    jf = np.asarray(merged['flux_j_templfit'], float); je = np.asarray(merged['fluxerr_j_templfit'], float)
    hf = np.asarray(merged['flux_h_templfit'], float); he = np.asarray(merged['fluxerr_h_templfit'], float)

    g, r, i, z = abmag(gf), abmag(rf), abmag(if_), abmag(zf)
    y, j, h    = abmag(yf), abmag(jf), abmag(hf)
    snr_ok = (
        (snr(gf, ge)  > SNR_MIN) & (snr(rf, re_) > SNR_MIN) & (snr(if_, ie) > SNR_MIN) &
        (snr(zf, ze)  > SNR_MIN) & (snr(yf, ye)  > SNR_MIN) & (snr(jf, je) > SNR_MIN) &
        (snr(hf, he)  > SNR_MIN)
    )
    mag_ok = (
        (g > MAG_MIN) & (r > MAG_MIN) & (i > MAG_MIN) & (z > MAG_MIN) &
        (y > MAG_MIN) & (j > MAG_MIN) & (h > MAG_MIN) &
        (g < MAG_MAX) & (r < MAG_MAX) & (i < MAG_MAX) & (z < MAG_MAX) &
        (y < MAG_MAX) & (j < MAG_MAX) & (h < MAG_MAX)
    )
    merged['gr'], merged['ri'], merged['iz'] = g - r, r - i, i - z
    merged['zy'], merged['yj'], merged['jh'] = z - y, y - j, j - h
    merged['mH_obs_mer'] = h
    merged['phot_ok'] = snr_ok & mag_ok
    merged['phot_system'] = gsys

    # EDF-N AGN should be tagged 'north'. Restrict to north only as a safety net.
    north = (merged['phot_system'] == 'north') & merged['phot_ok']
    out = merged[north]
    print(f'EDF-N AGN through north-system + 7-band S/N>{SNR_MIN} gate: {len(out):,}')
    return out


def attach_mH(tab, id_col):
    """Merge the PHYS H-band flux, fall back to MER H, store mH_obs."""
    phys = Table.read(os.path.join(CACHE_DIR, 'tutorial_phys.fits'))
    phys_h = as_float_array(phys['flux_h_total_corrected'])
    phys = phys[np.isfinite(phys_h) & (phys_h > 0)].copy()
    phys['mH_obs_phys'] = abmag(as_float_array(phys['flux_h_total_corrected']))
    if id_col != 'object_id':
        phys.rename_column('object_id', id_col)
    tab = join(tab, phys[[id_col, 'mH_obs_phys']], keys=id_col, join_type='left')
    mH_phys = as_float_array(tab['mH_obs_phys'])
    mH_mer  = as_float_array(tab['mH_obs_mer'])
    mH = mH_phys.copy()
    fallback = ~np.isfinite(mH) & valid_mag(mH_mer)
    mH[fallback] = mH_mer[fallback]
    tab['mH_obs'] = mH
    return tab


def load_training():
    """Load the cached 50k-row Q1 training sample and apply the same photometric
    quality gates as the original tutorial. The existing cache happens to be
    fully south-system (DECam) because the TAP ORDER BY object_id pulled a
    contiguous southern block; EDF-N AGN are projected onto this SOM with the
    same cross-bandpass compromise the original notebook accepts. Cleaner
    long-term fix is a north-system-only training query, but it is not needed
    for the workshop pedagogical signal."""
    t = Table.read(os.path.join(CACHE_DIR, 'tutorial_mer_training.fits'))

    gf, ge, gsys = coalesce(t, 'flux_g_ext_hsc_templfit',       'flux_g_ext_decam_templfit')
    rf, re_, _   = coalesce(t, 'flux_r_ext_megacam_templfit',   'flux_r_ext_decam_templfit')
    if_, ie, _   = coalesce(t, 'flux_i_ext_panstarrs_templfit', 'flux_i_ext_decam_templfit')
    zf, ze, _    = coalesce(t, 'flux_z_ext_hsc_templfit',       'flux_z_ext_decam_templfit')
    yf = np.asarray(t['flux_y_templfit'], float); ye = np.asarray(t['fluxerr_y_templfit'], float)
    jf = np.asarray(t['flux_j_templfit'], float); je = np.asarray(t['fluxerr_j_templfit'], float)
    hf = np.asarray(t['flux_h_templfit'], float); he = np.asarray(t['fluxerr_h_templfit'], float)

    g, r, i, z = abmag(gf), abmag(rf), abmag(if_), abmag(zf)
    y, j, h    = abmag(yf), abmag(jf), abmag(hf)
    snr_ok = (
        (snr(gf, ge)  > SNR_MIN) & (snr(rf, re_) > SNR_MIN) & (snr(if_, ie) > SNR_MIN) &
        (snr(zf, ze)  > SNR_MIN) & (snr(yf, ye)  > SNR_MIN) & (snr(jf, je) > SNR_MIN) &
        (snr(hf, he)  > SNR_MIN)
    )
    mag_ok = (
        (g > MAG_MIN) & (r > MAG_MIN) & (i > MAG_MIN) & (z > MAG_MIN) &
        (y > MAG_MIN) & (j > MAG_MIN) & (h > MAG_MIN) &
        (g < MAG_MAX) & (r < MAG_MAX) & (i < MAG_MAX) & (z < MAG_MAX) &
        (y < MAG_MAX) & (j < MAG_MAX) & (h < MAG_MAX)
    )
    t['gr'], t['ri'], t['iz'] = g - r, r - i, i - z
    t['zy'], t['yj'], t['jh'] = z - y, y - j, j - h
    t['mH_obs_mer'] = h
    t['phot_ok'] = snr_ok & mag_ok
    t['phot_system'] = gsys
    out = t[t['phot_ok']]
    n_north = int((out['phot_system'] == 'north').sum())
    n_south = int((out['phot_system'] == 'south').sum())
    print(f'Training galaxies through gates: {len(out):,} (north {n_north:,}, south {n_south:,})')
    return out


def train_som(train_tab):
    X = np.vstack([train_tab[c] for c in COLOR_NAMES]).T
    keep = np.all(np.isfinite(X), axis=1)
    X = X[keep]
    train_tab = train_tab[keep]

    mu = np.nanmean(X, axis=0)
    sigma = np.nanstd(X, axis=0)
    sigma[sigma == 0] = 1.0
    Xn = (X - mu) / sigma

    som = MiniSom(MSZ[0], MSZ[1], len(COLOR_NAMES),
                  sigma=10.0, learning_rate=0.1, random_seed=RANDOM_SEED)
    som.random_weights_init(Xn)
    print(f'Training SOM on {len(Xn):,} galaxies, {N_ITER:,} iterations...')
    som.train_random(Xn, N_ITER)

    weights = som.get_weights()
    bmu = np.array([som.winner(v) for v in Xn])
    qe  = np.linalg.norm(Xn - weights[bmu[:, 0], bmu[:, 1]], axis=1)
    qe_thresh = float(np.percentile(qe, QE_PCTILE))
    print(f'Training QE (mean / median / {QE_PCTILE}th pct): '
          f'{qe.mean():.3f} / {np.median(qe):.3f} / {qe_thresh:.3f}')

    train_tab['bmu_row'] = bmu[:, 0]
    train_tab['bmu_col'] = bmu[:, 1]
    train_tab['qe'] = qe
    return som, mu, sigma, qe_thresh, train_tab


def project_agn(agn, som, mu, sigma, qe_thresh):
    ok, bmu, qe = project(agn, som, mu, sigma)
    out = agn[ok].copy()
    out['bmu_row'] = bmu[:, 0]
    out['bmu_col'] = bmu[:, 1]
    out['qe'] = qe
    out['qe_inlier'] = qe <= qe_thresh
    print(f'AGN projected: {len(out):,}; inliers (QE<= {qe_thresh:.2f}): '
          f'{int(out["qe_inlier"].sum()):,}')
    return out


def select_picks(agn_proj):
    """Pick three AGN whose spectra are already in the local cache, spanning
    distinct selection-flag combinations and SOM neighborhoods. Returns a
    list of dicts ready for JSON serialisation.
    """
    cached_ids = []
    for fn in os.listdir(SPECTRA_CACHE_DIR):
        if fn.endswith('.npz'):
            try:
                cached_ids.append(int(fn[:-4]))
            except ValueError:
                continue
    cached_ids = np.array(cached_ids, dtype=np.int64)
    mask = np.isin(np.asarray(agn_proj['object_id_euclid'], dtype=np.int64), cached_ids)
    candidates = agn_proj[mask & agn_proj['qe_inlier']]
    print(f'AGN with cached spectra also projected as inliers: {len(candidates)}')

    # Score each candidate by how many flags it has, prefer spec-z + diverse BMU.
    if len(candidates) == 0:
        raise RuntimeError('No cached spectra match the projected AGN catalog.')

    flag_count = np.zeros(len(candidates), dtype=int)
    for f in ALL_FLAGS:
        flag_count = flag_count + (np.asarray(candidates[f]) == 1).astype(int)
    candidates['n_flags'] = flag_count

    picks = []
    # Order: spec-z first, then by flag count descending, for reproducibility.
    order = np.lexsort((-flag_count, np.asarray(candidates['z_source']) != 'spec'))
    used_bmus = set()
    for idx in order:
        row = candidates[int(idx)]
        bmu = (int(row['bmu_row']), int(row['bmu_col']))
        # Skip if too close to a previously-picked BMU (keeps the three diverse).
        too_close = any(abs(bmu[0] - b[0]) <= 3 and abs(bmu[1] - b[1]) <= 3 for b in used_bmus)
        if too_close:
            continue
        used_bmus.add(bmu)
        flagged = [f for f in ALL_FLAGS if int(row[f]) == 1]
        picks.append({
            'object_id_euclid': int(row['object_id_euclid']),
            'ra': float(row['ra_euclid']),
            'dec': float(row['dec_euclid']),
            'z_final': float(row['z_final']),
            'z_source': str(row['z_source']),
            'bmu_row': int(row['bmu_row']),
            'bmu_col': int(row['bmu_col']),
            'flags_set': flagged,
        })
        if len(picks) >= 3:
            break
    if len(picks) < 3:
        # Fallback: drop the spacing constraint to fill the slate.
        for idx in order:
            row = candidates[int(idx)]
            oid = int(row['object_id_euclid'])
            if any(p['object_id_euclid'] == oid for p in picks):
                continue
            flagged = [f for f in ALL_FLAGS if int(row[f]) == 1]
            picks.append({
                'object_id_euclid': oid,
                'ra': float(row['ra_euclid']),
                'dec': float(row['dec_euclid']),
                'z_final': float(row['z_final']),
                'z_source': str(row['z_source']),
                'bmu_row': int(row['bmu_row']),
                'bmu_col': int(row['bmu_col']),
                'flags_set': flagged,
            })
            if len(picks) >= 3:
                break
    return picks


def main():
    print('=' * 72)
    print('Build workshop artifacts (EDF-N only)')
    print('=' * 72)

    agn = load_edfn_agn()
    agn = attach_z_and_phot(agn)
    agn = attach_mH(agn, id_col='object_id_euclid')

    train = load_training()
    train = attach_mH(train, id_col='object_id')

    som, mu, sigma, qe_thresh, train_proj = train_som(train)
    agn_proj = project_agn(agn, som, mu, sigma, qe_thresh)

    picks = select_picks(agn_proj)
    for p in picks:
        print(f'  pick: {p["object_id_euclid"]}  bmu=({p["bmu_row"]},{p["bmu_col"]})  '
              f'z={p["z_final"]:.3f} ({p["z_source"]})  '
              f'flags={p["flags_set"][:3]}{"..." if len(p["flags_set"]) > 3 else ""}')

    # --- write outputs ------------------------------------------------------
    agn_out = os.path.join(CACHE_DIR, 'workshop_agn_edfn.fits')
    train_out = os.path.join(CACHE_DIR, 'workshop_train_edfn.fits')
    som_out = os.path.join(CACHE_DIR, 'workshop_som.pkl')
    picks_out = os.path.join(CACHE_DIR, 'workshop_picks.json')

    agn_proj.write(agn_out, overwrite=True)
    train_proj.write(train_out, overwrite=True)
    with open(som_out, 'wb') as fh:
        pickle.dump({
            'som': som,
            'mu': mu,
            'sigma': sigma,
            'qe_thresh': qe_thresh,
            'msz': MSZ,
            'color_names': list(COLOR_NAMES),
            'z_window': list(Z_WIN),
            'snr_min': SNR_MIN,
        }, fh, protocol=4)
    with open(picks_out, 'w') as fh:
        json.dump(picks, fh, indent=2)

    for path in [agn_out, train_out, som_out, picks_out]:
        size_kb = os.path.getsize(path) / 1024
        print(f'  wrote {path}  ({size_kb:,.1f} KB)')


if __name__ == '__main__':
    main()
