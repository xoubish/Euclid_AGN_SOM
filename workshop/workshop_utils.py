"""Routine helpers for ``Euclid_AGN_Tutorial.ipynb``.

The notebook imports everything via ``from agn_tutorial_utils import *`` and
keeps the scientific narrative (config knobs, ``plot_*`` functions, analysis
cells) inline. Boilerplate that you can name without explanation -- data IO,
table utilities, SOM grid arithmetic, plotting helpers, cache management,
SAS spectrum fetching, SPE Halpha querying -- lives here.

Module-level constants below act as defaults; the notebook may override them
via kwargs at call sites or by reassigning the attribute on this module.
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import scipy.interpolate
from scipy.ndimage import gaussian_filter, median_filter

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize, PowerNorm

from astropy import units as u
from astropy.cosmology import Planck18 as COSMO
from astropy.io import fits
from astropy.table import QTable, Table, join, vstack
from astropy.units import UnitsWarning

try:
    from tqdm.auto import tqdm
except ImportError:
    def tqdm(x, **kwargs):
        return x

from astroquery.ipac.irsa import Irsa


# ---------------------------------------------------------------------------
# Module-level constants. Notebook config mirrors these for visibility; the
# helpers use the values defined here.
# ---------------------------------------------------------------------------
MSZ = (30, 30)
RANDOM_SEED = 0
BATCH = 1000
MAG_MIN, MAG_MAX = 12, 30
COLOR_NAMES = ('gr', 'ri', 'iz', 'zy', 'yj', 'jh')

# AGN density overlay
SMOOTH_SIGMA = 1.35
DENSITY_FLOOR = 0.08
DENSITY_PCTILE = 18
DENSITY_VMAX_PCTILE = 98.0
DENSITY_N_CONTOURS = 8
AGN_JITTER = 0.34
AGN_OVERLAY_COLOR = '#064d2c'

# Method-excess overlay
METHOD_EXCESS_FLOOR = 1.25
METHOD_EXCESS_VMAX_PCTILE = 98
METHOD_MIN_EXPECTED = 0.03
METHOD_CONTOUR_COLOR = '#064d2c'

# SPE Halpha query
SPE_HALPHA_BATCH = 900
HALPHA_SNR_MIN = 3.0
SPE_LINES_TABLE = 'euclid_q1_spe_lines_line_features'
SPE_HALPHA_QUERY_COLS = [
    'spe_rank',
    'spe_line_name',
    'spe_line_snr_gf',
    'spe_line_flux_gf',
    'spe_line_flux_err_gf',
    'spe_line_central_wl_gf',
    'spe_line_ew_gf',
]

# Spectra cache / SAS
SPECTRA_CACHE_DIR = 'data/cache/irsa_spectra_agn'
SPECTRA_BUCKET = 'nasa-irsa-euclid-q1'
SPECTRA_ASSOC_TABLE = 'euclid.objectid_spectrafile_association_q1'

# Halpha-complex rest-frame wavelengths (Angstrom)
HALPHA_REST = 6562.8
NII_REST = [6548.05, 6583.45]
SII_REST = [6716.44, 6730.82]


# ---------------------------------------------------------------------------
# Colormaps used by overlays
# ---------------------------------------------------------------------------
AGN_CMAP = LinearSegmentedColormap.from_list('agn_density_green', [
    (0.82, 0.96, 0.70, 0.34),
    (0.66, 0.88, 0.58, 0.58),
    (0.28, 0.68, 0.35, 0.78),
    (0.00, 0.33, 0.18, 0.96),
], N=256)

METHOD_CMAP = LinearSegmentedColormap.from_list('method_excess_green', [
    (1.0, 1.0, 1.0, 0.00),
    (0.86, 0.96, 0.78, 0.34),
    (0.36, 0.74, 0.42, 0.72),
    (0.00, 0.33, 0.18, 0.95),
], N=256)

HALPHA_DERIVED_CMAP = plt.get_cmap('magma').copy()
HALPHA_DERIVED_CMAP.set_bad(color='white', alpha=1.0)


# Some Euclid Q1 spectrum FITS tables use the metadata unit string 'Number'.
# It is not a standard FITS unit but is harmless for the columns used here.
try:
    _NUMBER_UNIT = u.def_unit('Number')
    u.add_enabled_units([_NUMBER_UNIT])
except ValueError:
    pass
warnings.filterwarnings(
    'ignore',
    message=r".*'Number' did not parse as fits unit.*",
    category=UnitsWarning,
)

os.makedirs(SPECTRA_CACHE_DIR, exist_ok=True)


# ===========================================================================
# Photometry helpers
# ===========================================================================
def abmag(flux_uJy):
    """Convert microJansky fluxes to AB magnitudes.

    Non-finite or non-positive values return -99.
    """
    f = np.asarray(flux_uJy, dtype=float)
    bad = ~np.isfinite(f) | (f <= 0)
    out = np.full_like(f, -99.0)
    out[~bad] = -2.5 * np.log10(f[~bad]) + 23.9
    return out


def snr(flux, err):
    """Per-element signal-to-noise ratio; zero / non-finite errors return 0."""
    e = np.where(np.isfinite(err) & (err > 0), err, np.inf)
    return flux / e


def as_float_array(col, fill=np.nan):
    """Convert Astropy/Numpy columns to float arrays, preserving masks."""
    return np.asarray(np.ma.filled(np.ma.asarray(col), fill), dtype=float)


def valid_mag(m, mag_min=MAG_MIN, mag_max=MAG_MAX):
    return np.isfinite(m) & (m > mag_min) & (m < mag_max)


def coalesce(t, primary, fallback):
    """Take `primary` flux/err where finite & positive, else `fallback`.

    Returns (flux, err, system) where system labels each source 'north' or
    'south' based on which photometric system was used.
    """
    pf = np.asarray(t[primary], float)
    pe = np.asarray(t[primary.replace('flux_', 'fluxerr_')], float)
    ff = np.asarray(t[fallback], float)
    fe = np.asarray(t[fallback.replace('flux_', 'fluxerr_')], float)
    ok = np.isfinite(pf) & (pf > 0)
    return np.where(ok, pf, ff), np.where(ok, pe, fe), np.where(ok, 'north', 'south')


def fill_nan_nearest(m):
    """Nearest-neighbour fill of NaN cells in a 2D SOM map (display only)."""
    arr = np.asarray(m, dtype=float)
    mask = np.isfinite(arr)
    if not mask.any():
        return arr.copy()
    xx, yy = np.meshgrid(np.arange(arr.shape[1]), np.arange(arr.shape[0]))
    xym = np.vstack((np.ravel(xx[mask]), np.ravel(yy[mask]))).T
    return scipy.interpolate.NearestNDInterpolator(
        xym, np.ravel(arr[mask])
    )(np.ravel(xx), np.ravel(yy)).reshape(xx.shape)


# ===========================================================================
# IRSA TAP batched query
# ===========================================================================
def tap_query_with_retry(adql, retries=5, base_delay=2.0, async_job=False):
    """Run a TAP query with exponential backoff for transient failures.

    IRSA's TAP endpoint occasionally returns 502/503/504 or VOTable parse
    errors under load; these are not query bugs and clear on retry. Set
    ``async_job=True`` for long queries that need the async TAP endpoint.
    """
    import time
    last_exc = None
    for attempt in range(retries):
        try:
            return Irsa.query_tap(adql, async_job=async_job).to_table()
        except Exception as exc:
            last_exc = exc
            msg = str(exc)
            transient = any(code in msg for code in ('502', '503', '504', 'Bad Gateway',
                                                     'Gateway Time-out', 'Service Unavailable',
                                                     'VOTABLE'))
            if not transient or attempt == retries - 1:
                raise
            delay = base_delay * (2 ** attempt)
            warnings.warn(f'TAP transient failure (attempt {attempt+1}/{retries}): {exc!r}; '
                          f'retrying in {delay:.1f}s', RuntimeWarning)
            time.sleep(delay)
    raise last_exc  # unreachable, kept for static analysis


def batched_query(table, cols, ids, batch=BATCH, id_col='object_id', desc=None):
    """Run an ADQL ``WHERE id IN (...)`` query in batches.

    IRSA's TAP service caps SQL length per request, so large ID lists must
    be split. One query per batch, results vstacked, progress via tqdm.
    Each batch is wrapped with ``tap_query_with_retry`` so transient
    5xx responses self-heal without losing earlier batches' progress.
    """
    if len(ids) == 0:
        return Table(names=[id_col] + cols)
    cols_sql = ', '.join([id_col] + cols)
    out, n_batches = [], (len(ids) + batch - 1) // batch
    iterator = tqdm(range(0, len(ids), batch), total=n_batches,
                    desc=desc or table, unit='batch')
    for k in iterator:
        in_list = ','.join(str(int(x)) for x in ids[k:k + batch])
        adql = f'SELECT {cols_sql} FROM {table} WHERE {id_col} IN ({in_list})'
        out.append(tap_query_with_retry(adql))
    return vstack(out)


# ===========================================================================
# Astropy Table utilities
# ===========================================================================
def table_float(tab, col):
    return np.asarray(np.ma.filled(tab[col], np.nan), dtype=float)


def table_int(tab, col, fill=-1):
    return np.asarray(np.ma.filled(tab[col], fill), dtype=np.int64)


def ordered_unique_int(values):
    out, seen = [], set()
    for value in values:
        oid = int(value)
        if oid not in seen:
            seen.add(oid)
            out.append(oid)
    return out


def row_lookup_by_object_id(tab):
    return {int(row['object_id_euclid']): row for row in tab}


def group_mask(tab, flags):
    """Boolean OR across a list of {0,1} flag columns; missing columns are False."""
    m = np.zeros(len(tab), dtype=bool)
    for f in flags:
        if f in tab.colnames:
            m |= (tab[f] == 1)
    return m


# ===========================================================================
# SOM grid arithmetic
# ===========================================================================
def median_map(values, bmu_xy, msz=MSZ):
    """Median statistic per SOM cell. Empty cells return NaN."""
    m = np.full(msz, np.nan)
    for ix in range(msz[0]):
        for iy in range(msz[1]):
            sel = (bmu_xy[:, 0] == ix) & (bmu_xy[:, 1] == iy)
            if sel.any():
                v = as_float_array(values[sel])
                v = v[np.isfinite(v)]
                if v.size:
                    m[ix, iy] = np.median(v)
    return m


def project(table, som, mu, sigma, color_names=COLOR_NAMES):
    """Project a source table onto a trained SOM.

    Each source's colors are z-score-normalised with the *training* mean
    and standard deviation, then mapped to the best-matching unit (BMU).
    Sources with non-finite colors in any band are dropped.

    Returns ``(ok, bmu, qe)`` where ``ok`` masks the input rows, ``bmu`` is
    an array of (row, col) coords, and ``qe`` is the per-source quantization
    error (distance to the BMU weight in normalised color space).
    """
    X = np.vstack([table[c] for c in color_names]).T
    ok = np.all(np.isfinite(X), axis=1)
    Xn = (X[ok] - mu) / sigma
    bmu = np.array([som.winner(v) for v in Xn])
    weights = som.get_weights()
    qe = np.linalg.norm(Xn - weights[bmu[:, 0], bmu[:, 1]], axis=1)
    return ok, bmu, qe


def bmu_density(bmu_xy, mask=None, msz=MSZ):
    """2D histogram of BMU coordinates."""
    xs = bmu_xy[:, 0] if mask is None else bmu_xy[mask, 0]
    ys = bmu_xy[:, 1] if mask is None else bmu_xy[mask, 1]
    if len(xs) == 0:
        return None
    H, *_ = np.histogram2d(
        xs, ys, bins=[np.arange(msz[0] + 1), np.arange(msz[1] + 1)]
    )
    return H


def jitter_bmu(bmu_xy, jitter=AGN_JITTER, seed=RANDOM_SEED):
    """Return (x, y) display coords for BMU crosses, jittered within cells."""
    rng = np.random.default_rng(seed)
    y = bmu_xy[:, 0].astype(float) + rng.uniform(-jitter, jitter, len(bmu_xy))
    x = bmu_xy[:, 1].astype(float) + rng.uniform(-jitter, jitter, len(bmu_xy))
    return x, y


def cell_count_grid(bmu_xy, shape=MSZ):
    grid = np.zeros(shape, dtype=float)
    if len(bmu_xy) == 0:
        return grid
    for r, c in np.asarray(bmu_xy, dtype=int):
        if 0 <= r < shape[0] and 0 <= c < shape[1]:
            grid[r, c] += 1.0
    return grid


def detection_count_grid(det_table, shape=MSZ):
    grid = np.zeros(shape, dtype=float)
    if len(det_table) == 0:
        return grid
    for r, c in zip(det_table['som_row'], det_table['som_col']):
        grid[int(r), int(c)] += 1.0
    return grid


def mean_metric_grid_including_zeros(sample_bmu, det_table, value_col,
                                     value_transform=None, shape=MSZ):
    """Per-cell mean of `value_col` from `det_table`, dividing by the sample
    total (so non-detections enter as zeros). Empty cells return NaN."""
    total = cell_count_grid(sample_bmu, shape=shape)
    summed = np.zeros(shape, dtype=float)
    if len(det_table) > 0 and value_col in det_table.colnames:
        values = table_float(det_table, value_col)
        if value_transform is not None:
            values = value_transform(values)
        good = np.isfinite(values) & (values > 0)
        for r, c, value in zip(det_table['som_row'][good],
                               det_table['som_col'][good],
                               values[good]):
            summed[int(r), int(c)] += float(value)
    mean = np.full(shape, np.nan, dtype=float)
    occupied = total > 0
    mean[occupied] = summed[occupied] / total[occupied]
    return mean, total


def detection_fraction_grid(sample_bmu, det_table, shape=MSZ):
    total = cell_count_grid(sample_bmu, shape=shape)
    detected = detection_count_grid(det_table, shape=shape)
    fraction = np.full(shape, np.nan, dtype=float)
    occupied = total > 0
    fraction[occupied] = detected[occupied] / total[occupied]
    return fraction, total, detected


def ensure_halpha_luminosity(tab):
    """Add a 'halpha_luminosity' column (erg/s) from halpha_flux + z_final."""
    tab = tab.copy()
    if 'halpha_luminosity' in tab.colnames:
        return tab
    luminosity = np.full(len(tab), np.nan, dtype=float)
    flux = table_float(tab, 'halpha_flux')
    redshift = table_float(tab, 'z_final')
    good = np.isfinite(flux) & (flux > 0) & np.isfinite(redshift) & (redshift > 0)
    if np.any(good):
        dl_cm = COSMO.luminosity_distance(redshift[good]).to_value(u.cm)
        luminosity[good] = 4.0 * np.pi * dl_cm**2 * flux[good]
    tab['halpha_luminosity'] = luminosity
    return tab


# ===========================================================================
# Plotting helpers (overlays, colorbars, normalisations)
# ===========================================================================
def smooth_density(H, smooth=SMOOTH_SIGMA):
    if H is None:
        return None
    return gaussian_filter(H.astype(float), sigma=smooth) if smooth > 0 else H.astype(float)


def overlay(ax, H, smooth=SMOOTH_SIGMA, floor=DENSITY_FLOOR,
            pctile=DENSITY_PCTILE, vmax_pctile=DENSITY_VMAX_PCTILE):
    """Plot only high-density AGN regions over an existing SOM background."""
    if H is None:
        return None, np.nan, np.nan
    Hs = gaussian_filter(H.astype(float), sigma=smooth) if smooth > 0 else H.astype(float)
    positive = Hs[Hs > 0]
    if positive.size == 0:
        return None, np.nan, np.nan
    vmin = max(float(floor), float(np.nanpercentile(positive, pctile)))
    Hm = np.ma.masked_where(Hs < vmin, Hs)
    if Hm.count() == 0:
        vmin = 0.5 * float(np.nanmax(positive))
        Hm = np.ma.masked_where(Hs < vmin, Hs)
    vmax = max(float(np.nanpercentile(positive, vmax_pctile)), vmin * 1.35)
    im = ax.imshow(Hm, origin='lower', cmap=AGN_CMAP,
                   norm=LogNorm(vmin=vmin, vmax=vmax), interpolation='bilinear')
    ax.contour(Hs, levels=np.geomspace(vmin, vmax, DENSITY_N_CONTOURS),
               colors=AGN_OVERLAY_COLOR, linewidths=0.55, alpha=0.68,
               origin='lower')
    return im, vmin, vmax


def overlay_method_excess(ax, H_group, H_union, floor=METHOD_EXCESS_FLOOR,
                          min_expected=METHOD_MIN_EXPECTED,
                          vmax_pctile=METHOD_EXCESS_VMAX_PCTILE):
    """Overlay where a method is enriched relative to the union AGN distribution."""
    if H_group is None or H_union is None or H_group.sum() <= 0 or H_union.sum() <= 0:
        return None
    Hg = smooth_density(H_group)
    Hu = smooth_density(H_union)
    expected = Hu * (float(H_group.sum()) / float(H_union.sum()))
    ratio = np.full_like(Hg, np.nan, dtype=float)
    good = expected >= min_expected
    ratio[good] = (Hg[good] + 0.02) / (expected[good] + 0.02)
    shown = ratio[np.isfinite(ratio) & (ratio >= floor)]
    if shown.size == 0:
        return None
    vmax = max(float(np.nanpercentile(shown, vmax_pctile)), floor * 1.2)
    ratio_m = np.ma.masked_where(~np.isfinite(ratio) | (ratio < floor), ratio)
    im = ax.imshow(ratio_m, origin='lower', cmap=METHOD_CMAP,
                   norm=LogNorm(vmin=floor, vmax=vmax), interpolation='bilinear')
    levels = np.geomspace(floor, vmax, 4)
    ax.contour(ratio, levels=levels, colors=METHOD_CONTOUR_COLOR, linewidths=0.35,
               alpha=0.35, origin='lower')
    return im


def add_panel_colorbar(fig, ax, im):
    """Compact per-panel colorbar for independently scaled method maps."""
    cax = ax.inset_axes([1.02, 0.16, 0.035, 0.68])
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=7, length=2)
    cbar.set_label('excess', fontsize=7)
    return cbar


def format_empty_som_axis(ax, msz=MSZ):
    ax.set_facecolor('white')
    ax.set_xlim(-0.5, msz[1] - 0.5); ax.set_ylim(-0.5, msz[0] - 0.5)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xticks(np.arange(-0.5, msz[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, msz[0], 1), minor=True)
    ax.grid(which='minor', color='0.90', linewidth=0.25)
    ax.tick_params(which='minor', bottom=False, left=False)


def power_norm_from_grids(*grids, gamma=0.35, pct=98):
    values = np.concatenate([
        g[np.isfinite(g)].ravel() for g in grids if np.isfinite(g).any()
    ])
    positive = values[values > 0]
    if positive.size == 0:
        vmax = 1.0
    else:
        vmax = float(np.nanpercentile(positive, pct))
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = float(np.nanmax(positive))
    return PowerNorm(gamma=gamma, vmin=0.0, vmax=vmax)


def halpha_flux_scale_norm(*arrays, gamma=0.35, pct=98):
    vals = []
    for arr in arrays:
        arr = np.asarray(arr, dtype=float)
        vals.append(arr[np.isfinite(arr)])
    vals = np.concatenate([v for v in vals if v.size]) if any(v.size for v in vals) else np.array([])
    if vals.size == 0:
        return None
    positive = vals[vals > 0]
    if positive.size == 0:
        return Normalize(vmin=0.0, vmax=1.0)
    vmax = np.nanpercentile(positive, pct) if positive.size > 2 else np.nanmax(positive)
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(np.nanmax(positive))
    return PowerNorm(gamma=gamma, vmin=0.0, vmax=vmax)


def log_flux_values(tab):
    if len(tab) == 0 or 'halpha_flux' not in tab.colnames:
        return np.array([], dtype=float)
    flux = np.asarray(tab['halpha_flux'], dtype=float)
    good = np.isfinite(flux) & (flux > 0)
    return np.log10(flux[good])


# ===========================================================================
# Spectra: cache + IRSA Spectra-Association + on-the-fly FITS read
# ===========================================================================
def association_to_s3_uri(path):
    key = str(path).replace('api/spectrumdm/convert/euclid/', '').split('?')[0].lstrip('/')
    if key.startswith('s3://'):
        return key
    return f's3://{SPECTRA_BUCKET}/{key}'


def load_cached_spectrum(object_id):
    cache_npz = os.path.join(SPECTRA_CACHE_DIR, f'{int(object_id)}.npz')
    if not os.path.exists(cache_npz):
        return None
    data = np.load(cache_npz)
    return {
        'wave': data['wave'] * u.angstrom,
        'flux': data['flux'] * (u.erg / u.s / u.cm**2 / u.angstrom),
        'error': data['error'] * (u.erg / u.s / u.cm**2 / u.angstrom),
        'object_id': int(object_id),
    }


def save_cached_spectrum(object_id, wave, flux, error):
    np.savez(
        os.path.join(SPECTRA_CACHE_DIR, f'{int(object_id)}.npz'),
        wave=np.asarray(wave.to_value(u.angstrom), dtype=float),
        flux=np.asarray(flux.value, dtype=float),
        error=np.asarray(error.value, dtype=float),
    )


def query_spectrum_associations(object_ids):
    object_ids = ordered_unique_int(object_ids)
    if len(object_ids) == 0:
        return Table(names=['objectid', 'path', 'hdu'])
    id_list = ','.join(map(str, object_ids))
    query = f'''
    SELECT objectid, path, hdu
    FROM {SPECTRA_ASSOC_TABLE}
    WHERE objectid IN ({id_list})
    '''
    try:
        return tap_query_with_retry(query)
    except Exception as err:
        print(f'Spectrum association query failed: {err}')
        return Table(names=['objectid', 'path', 'hdu'])


def read_spectra_from_associations(assoc, wanted_ids, spectra=None, n_keep=8):
    """Open the per-pointing FITS files referenced by an association table
    and extract the requested object IDs. Successful reads are cached to
    ``SPECTRA_CACHE_DIR`` via ``save_cached_spectrum``."""
    spectra = {} if spectra is None else dict(spectra)
    wanted = set(ordered_unique_int(wanted_ids))
    groups = {}
    skipped_hdu = 0
    for row in assoc:
        try:
            object_id = int(row['objectid'])
            hdu_index = int(row['hdu'])
        except Exception:
            skipped_hdu += 1
            continue
        if object_id not in wanted or object_id in spectra:
            continue
        s3_uri = association_to_s3_uri(row['path'])
        groups.setdefault(s3_uri, []).append((object_id, hdu_index))
    if skipped_hdu:
        print(f'Skipped {skipped_hdu} association rows with missing/invalid HDU indices.')

    for s3_uri, items in tqdm(list(groups.items()),
                              desc='Opening spectrum FITS files', unit='file'):
        if len(spectra) >= n_keep:
            break
        try:
            with fits.open(s3_uri, fsspec_kwargs={'anon': True}, lazy_load_hdus=True) as hdul:
                for object_id, hdu_index in items:
                    if len(spectra) >= n_keep:
                        break
                    try:
                        with warnings.catch_warnings():
                            warnings.filterwarnings(
                                'ignore',
                                message=r".*'Number' did not parse as fits unit.*",
                                category=UnitsWarning,
                            )
                            spec = QTable.read(hdul[hdu_index], format='fits')
                        header = hdul[hdu_index].header
                        fscale = header.get('FSCALE', 1.0)
                        wave = np.asarray(spec['WAVELENGTH'], dtype=float) * u.angstrom
                        signal = np.asarray(spec['SIGNAL'], dtype=float)
                        var = np.asarray(spec['VAR'], dtype=float)
                        mask = np.asarray(spec['MASK'])
                        valid = (mask % 2 == 0) & (mask < 64) & np.isfinite(signal) & np.isfinite(var) & (var > 0)
                        if np.sum(valid) < 20:
                            continue
                        wave = wave[valid]
                        flux = signal[valid] * fscale * u.erg / u.s / u.cm**2 / u.angstrom
                        error = np.sqrt(var[valid]) * fscale * flux.unit
                        spectra[object_id] = {'wave': wave, 'flux': flux, 'error': error, 'object_id': object_id}
                        save_cached_spectrum(object_id, wave, flux, error)
                    except Exception:
                        continue
        except Exception as err:
            print(f'Failed to open {s3_uri}: {err}')
    return spectra


def get_Q1_sir_spectra(tab, n=8, max_ids=1200, prefer_spec_z=True,
                       batch=250, cache_check=True):
    """Fetch up to ``n`` Q1 SIR 1D spectra for objects in ``tab``.

    Workflow: order rows (spec-z first when available), look up each object_id
    in the on-disk cache (``SPECTRA_CACHE_DIR``), and for any misses run a
    batched ADQL query against ``SPECTRA_ASSOC_TABLE`` followed by SAS FITS
    reads. Returns a dict {object_id: spectrum_table}.

    Parameters
    ----------
    cache_check : bool, default True
        If True, reuse any per-object spectra already in the local cache
        before hitting IRSA. Set False to skip the cache and re-query.
    """
    order = np.arange(len(tab))
    if prefer_spec_z and 'z_source' in tab.colnames:
        spec_first = np.asarray(tab['z_source'], dtype=str) == 'spec'
        order = np.r_[np.where(spec_first)[0], np.where(~spec_first)[0]]
    object_ids = ordered_unique_int(np.asarray(tab['object_id_euclid'])[order])[:max_ids]

    spectra, remaining = {}, []
    if cache_check:
        for object_id in object_ids:
            cached = load_cached_spectrum(object_id)
            if cached is not None:
                spectra[object_id] = cached
                if len(spectra) >= n:
                    return {oid: spectra[oid] for oid in object_ids if oid in spectra}
            else:
                remaining.append(object_id)
    else:
        remaining = list(object_ids)

    for start in range(0, len(remaining), batch):
        if len(spectra) >= n:
            break
        batch_ids = remaining[start:start + batch]
        assoc = query_spectrum_associations(batch_ids)
        if len(assoc) == 0:
            continue
        spectra = read_spectra_from_associations(assoc, batch_ids, spectra=spectra, n_keep=n)

    return {oid: spectra[oid] for oid in object_ids if oid in spectra}


def continuum_normalized_rest_spectrum(spec, z, cont_window=101):
    """Subtract a running-median continuum and shift to rest-frame wavelength."""
    rest = spec['wave'].to_value(u.angstrom) / (1.0 + z)
    flux = spec['flux'].value
    good = np.isfinite(rest) & np.isfinite(flux)
    rest, flux = rest[good], flux[good]
    if rest.size < 20:
        return None
    order = np.argsort(rest)
    rest, flux = rest[order], flux[order]
    k = int(cont_window)
    if k % 2 == 0:
        k += 1
    if k >= flux.size:
        k = max(5, (flux.size // 2) * 2 - 1)
    continuum = median_filter(flux, size=k, mode='nearest')
    resid = flux - continuum
    scale = np.nanpercentile(np.abs(resid), 95)
    if not np.isfinite(scale) or scale <= 0:
        scale = np.nanstd(resid)
    if not np.isfinite(scale) or scale <= 0:
        scale = 1.0
    return rest, resid / scale


def mark_halpha_complex(ax, full_range=True):
    """Vertical guide lines for H-alpha, [NII], and (optionally) [SII]."""
    ax.axvline(HALPHA_REST, color='crimson', ls='--', lw=1.0)
    for line in NII_REST:
        ax.axvline(line, color='0.35', ls=':', lw=0.8)
    if full_range:
        for line in SII_REST:
            ax.axvline(line, color='0.55', ls=':', lw=0.7)


# ===========================================================================
# MER imaging cutouts via IRSA IBE (server-side cutout, returns small FITS)
# ===========================================================================
MER_SIA_COLLECTION = 'euclid_DpdMerBksMosaic'

# Map common band identifiers to the energy_bandpassname returned by SIA.
MER_BAND_ALIASES = {
    'VIS': 'VIS', 'Ie': 'VIS', 'I_E': 'VIS',
    'Y': 'Y', 'NISP-Y': 'Y',
    'J': 'J', 'NISP-J': 'J',
    'H': 'H', 'NISP-H': 'H',
    'G': 'G', 'g': 'G',
    'R': 'R', 'r': 'R',
    'I': 'I', 'i': 'I',
    'Z': 'Z', 'z': 'Z',
    'U': 'U', 'u': 'U',
}


def _mer_sia_lookup(ra_deg, dec_deg, search_radius_deg=0.005):
    """Run a single SIA query at (ra, dec) and return the row list once."""
    import astropy.units as u_
    rows = Irsa.query_sia(
        pos=(ra_deg * u_.deg, dec_deg * u_.deg, search_radius_deg * u_.deg),
        collection=MER_SIA_COLLECTION,
    )
    return rows


def get_Q1_mer_cutout(ra_deg, dec_deg, band='H', size_arcsec=10.0,
                     sia_rows=None, verbose=True):
    """Fetch a Q1 MER mosaic cutout via IRSA's IBE cutout service.

    Cloud-native pattern: SIA query returns the mosaic file URL on IRSA's
    cloud-hosted IBE; we append ``?center=ra,dec&size=Xarcsec`` to ask the
    server to return only the cutout (typically ~40 KB gzipped FITS).
    ``astropy.io.fits.open`` fetches and decompresses the URL in one step.

    Parameters
    ----------
    ra_deg, dec_deg : float
        Sky position in decimal degrees.
    band : str, default 'H'
        Filter to fetch. See ``MER_BAND_ALIASES`` for accepted names.
    size_arcsec : float, default 10
        Cutout side length in arcsec.
    sia_rows : astropy.table.Table or None
        Cached SIA result (for repeated calls at the same position). If
        ``None``, runs a fresh SIA query.
    verbose : bool
        Print the cutout URL being fetched.

    Returns
    -------
    dict with keys: ``data`` (2D array), ``header`` (FITS header),
    ``wcs`` (astropy.wcs.WCS), ``s3_uri`` (the underlying S3 path of the
    full mosaic), ``cutout_url`` (the HTTPS IBE cutout URL).
    """
    from astropy.io import fits as _fits
    from astropy.wcs import WCS as _WCS
    import json as _json

    target = MER_BAND_ALIASES.get(band, band)
    if sia_rows is None:
        sia_rows = _mer_sia_lookup(ra_deg, dec_deg)
    matches = [
        r for r in sia_rows
        if r['energy_bandpassname'] == target and r['dataproduct_subtype'] == 'science'
    ]
    if not matches:
        bands_seen = sorted({str(r['energy_bandpassname']) for r in sia_rows})
        raise ValueError(
            f'No {band!r} (= {target!r}) MER mosaic found at '
            f'({ra_deg:.4f}, {dec_deg:.4f}). Bands available here: {bands_seen}.'
        )
    row = matches[0]
    base_url = row['access_url']
    cloud = _json.loads(row['cloud_access'])['aws']
    s3_uri = f"s3://{cloud['bucket_name']}/{cloud['key']}"
    cutout_url = (
        f'{base_url}?center={ra_deg:.6f},{dec_deg:.6f}'
        f'&size={size_arcsec:.2f}arcsec'
    )
    if verbose:
        print(f'  band {target}: {cutout_url}')
    with _fits.open(cutout_url) as hdul:
        sci_idx = next(
            i for i, h in enumerate(hdul)
            if h.header.get('NAXIS', 0) >= 2 and h.header.get('NAXIS1', 0) > 0
        )
        data = np.asarray(hdul[sci_idx].data, dtype=float)
        header = hdul[sci_idx].header.copy()
    return {
        'band': target,
        'data': data,
        'header': header,
        'wcs': _WCS(header),
        's3_uri': s3_uri,
        'cutout_url': cutout_url,
    }


# ===========================================================================
# SPE Halpha line-feature catalog
# ===========================================================================
def empty_spe_halpha_table():
    return Table(names=['object_id'] + SPE_HALPHA_QUERY_COLS)


def get_Q1_spe_halpha(object_ids, cache_path,
                      batch=SPE_HALPHA_BATCH, max_ids=None,
                      cache_check=True):
    """Fetch Q1 SPE Halpha line-feature rows for ``object_ids``.

    Runs a batched ADQL query against ``SPE_LINES_TABLE`` filtered to
    ``spe_line_name = 'Halpha'``. Results are cached as ECSV at
    ``cache_path``; subsequent calls reuse the cache unless
    ``cache_check=False`` is passed.

    Parameters
    ----------
    cache_check : bool, default True
        If True (default) and ``cache_path`` exists, the cached table is
        returned and no TAP queries are issued. Set False to discard the
        cache and re-query.
    """
    object_ids = ordered_unique_int(object_ids)
    if max_ids is not None:
        object_ids = object_ids[:int(max_ids)]
    if not cache_check and os.path.exists(cache_path):
        os.remove(cache_path)
        print(f'cache_check=False: removed {cache_path}')
    if cache_check and os.path.exists(cache_path):
        tab = Table.read(cache_path, format='ascii.ecsv')
        print(f'Loaded cached SPE Halpha measurements: {len(tab):,} rows ({cache_path})')
        return tab
    if len(object_ids) == 0:
        return empty_spe_halpha_table()

    out = []
    cols_sql = ', '.join(['object_id'] + SPE_HALPHA_QUERY_COLS)
    n_batches = (len(object_ids) + batch - 1) // batch
    iterator = tqdm(range(0, len(object_ids), batch), total=n_batches,
                    desc=f'SPE Halpha rows: {cache_path}', unit='batch')
    for start in iterator:
        id_list = ','.join(str(int(x)) for x in object_ids[start:start + batch])
        adql = f'''
        SELECT {cols_sql}
        FROM {SPE_LINES_TABLE}
        WHERE object_id IN ({id_list})
          AND spe_line_name = 'Halpha'
        '''
        try:
            part = tap_query_with_retry(adql)
        except Exception as err:
            print(f'SPE Halpha query failed for batch {start // batch + 1}: {err}')
            continue
        if len(part) > 0:
            out.append(part)

    if out:
        tab = vstack(out, metadata_conflicts='silent')
        tab.write(cache_path, format='ascii.ecsv', overwrite=True)
        print(f'Cached SPE Halpha measurements: {len(tab):,} rows ({cache_path})')
    else:
        tab = empty_spe_halpha_table()
        print(f'No SPE Halpha line rows found for {cache_path}.')
    return tab


def make_spe_meta_table(tab, id_col, z_col=None, z_source='phot', z_source_col=None):
    """Build a small (object_id_euclid, z_final, z_source) meta-table for joining."""
    meta = Table()
    meta['object_id_euclid'] = np.asarray(tab[id_col], dtype=np.int64)
    if z_col is not None and z_col in tab.colnames:
        meta['z_final'] = table_float(tab, z_col)
    else:
        meta['z_final'] = np.full(len(tab), np.nan)
    if z_source_col is not None and z_source_col in tab.colnames:
        meta['z_source'] = np.asarray(tab[z_source_col], dtype=str)
    else:
        meta['z_source'] = np.full(len(tab), z_source)
    return meta


def empty_best_halpha_table():
    return Table(names=[
        'object_id_euclid', 'halpha_flux', 'halpha_flux_err', 'halpha_snr',
        'halpha_luminosity', 'halpha_obs_wavelength', 'halpha_ew',
        'spe_rank', 'z_final', 'z_source',
    ])


def strongest_spe_halpha_per_object(line_table, sample_meta):
    """Pick the highest-S/N Halpha row per object and attach z/luminosity."""
    if len(line_table) == 0:
        return empty_best_halpha_table()
    tab = line_table.copy()
    flux = table_float(tab, 'spe_line_flux_gf')
    snr_arr = table_float(tab, 'spe_line_snr_gf')
    good = np.isfinite(flux) & (flux > 0) & np.isfinite(snr_arr)
    tab = tab[good]
    if len(tab) == 0:
        return empty_best_halpha_table()

    ids = table_int(tab, 'object_id')
    snr_arr = table_float(tab, 'spe_line_snr_gf')
    order = np.lexsort((-snr_arr, ids))
    tab = tab[order]
    ids = table_int(tab, 'object_id')
    _, first = np.unique(ids, return_index=True)
    best = tab[first]

    out = Table()
    out['object_id_euclid'] = table_int(best, 'object_id')
    out['halpha_flux'] = table_float(best, 'spe_line_flux_gf')
    out['halpha_flux_err'] = table_float(best, 'spe_line_flux_err_gf')
    out['halpha_snr'] = table_float(best, 'spe_line_snr_gf')
    out['halpha_obs_wavelength'] = table_float(best, 'spe_line_central_wl_gf')
    out['halpha_ew'] = table_float(best, 'spe_line_ew_gf')
    out['spe_rank'] = table_int(best, 'spe_rank')

    sample_meta = sample_meta[['object_id_euclid', 'z_final', 'z_source']]
    out = join(out, sample_meta, keys='object_id_euclid', join_type='left')

    luminosity = np.full(len(out), np.nan, dtype=float)
    flux = table_float(out, 'halpha_flux')
    redshift = table_float(out, 'z_final')
    good_lum = np.isfinite(flux) & (flux > 0) & np.isfinite(redshift) & (redshift > 0)
    if np.any(good_lum):
        dl_cm = COSMO.luminosity_distance(redshift[good_lum]).to_value(u.cm)
        luminosity[good_lum] = 4.0 * np.pi * dl_cm**2 * flux[good_lum]
    out['halpha_luminosity'] = luminosity
    return out


def add_bmu_to_halpha(halpha_table, sample_table, bmu_xy, sample_id_col,
                      snr_min=HALPHA_SNR_MIN, label='sample'):
    """Join SOM BMU coords onto the SPE-Halpha measurement table, with an
    SNR cut (falls back to all positive fluxes if the cut yields none)."""
    if len(halpha_table) == 0:
        print(f'No SPE Halpha measurements to project for {label}.')
        return Table()
    bmu_tab = Table()
    bmu_tab['object_id_euclid'] = np.asarray(sample_table[sample_id_col], dtype=np.int64)
    bmu_tab['som_row'] = bmu_xy[:, 0]
    bmu_tab['som_col'] = bmu_xy[:, 1]
    joined = join(halpha_table, bmu_tab, keys='object_id_euclid', join_type='inner')
    good = (
        np.isfinite(joined['halpha_flux']) & (joined['halpha_flux'] > 0) &
        np.isfinite(joined['halpha_snr']) & (joined['halpha_snr'] >= snr_min)
    )
    if np.sum(good) == 0:
        print(f'No {label} SPE Halpha measurements pass S/N >= {snr_min}; '
              f'showing all positive catalog fluxes instead.')
        good = np.isfinite(joined['halpha_flux']) & (joined['halpha_flux'] > 0)
    return joined[good]
