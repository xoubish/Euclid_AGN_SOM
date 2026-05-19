"""Generate ``Euclid_AGN_Workshop.ipynb`` from the cell sequence below.

The workshop notebook is a 30-minute hands-on derived from the full
``Euclid_AGN_Tutorial.ipynb``. It loads the pre-built artifacts created by
``build_workshop_artifacts.py`` (EDF-N AGN only, trained SOM, three pre-picked
candidates) and walks participants through two short fill-in-the-blank tasks:

  Task 1 - Plot one AGN selection method's density and excess on the SOM
  Task 2 - Project a chosen candidate, then pull data live from the cloud:
           1D SIR spectrum (direct S3), MER multi-band cutouts (IBE), and a
           live SPE Halpha catalog query (TAP). All three are cloud-native
           reads against IRSA infrastructure.

Re-run this script whenever the cell sequence is edited.
"""

from __future__ import annotations

import json
import os
from typing import List


def md(source: str) -> dict:
    return {'cell_type': 'markdown', 'metadata': {}, 'source': source}


def code(source: str) -> dict:
    return {
        'cell_type': 'code',
        'execution_count': None,
        'metadata': {},
        'outputs': [],
        'source': source,
    }


CELLS: List[dict] = []

# ---------------------------------------------------------------------------
# Title and orientation
# ---------------------------------------------------------------------------
CELLS.append(md(
    "# Euclid Q1 hands-on: AGN selection on a color SOM\n"
    "\n"
    "**AAS workshop:** *Euclid Data in the Cloud - Access, Analysis, and Science Opportunities* (90 min)\n"
    "\n"
    "This 30-minute hands-on uses a pre-trained color SOM and a pre-projected EDF-N AGN catalog so\n"
    "we can focus on the **cloud-access** half of the workshop title. The two tasks are:\n"
    "\n"
    "1. **Task 1** - Pick one AGN selection method and visualize where it sits on the color manifold,\n"
    "   relative to the union of all methods.\n"
    "2. **Task 2** - Pick an individual AGN candidate, project it, and verify it against three\n"
    "   independent **live cloud reads** against IRSA: its 1D Q1 SIR spectrum (direct S3 FITS read),\n"
    "   its multi-band MER imaging (IBE server-side cutout), and its SPE H$\\alpha$ line measurement\n"
    "   (live TAP query). None of those data products are pre-cached for you.\n"
    "\n"
    "**The cloud point.** Everything in Task 2 is fetched at workshop time from\n"
    "`s3://nasa-irsa-euclid-q1/...` or `https://irsa.ipac.caltech.edu/...`. No downloads happened\n"
    "in advance and none happen to disk now. Running this on Fornax it's region-local and\n"
    "instant; from a laptop it works too over the public HTTPS/anonymous-S3 endpoints. The same\n"
    "code scales unchanged to the full Q1 release and onward to DR1.\n"
    "\n"
    "**The pre-built bits.** Only two things are pre-built for you, both to save ~2 minutes:\n"
    "the trained SOM weights (`workshop_som.pkl`) and the EDF-N AGN catalog with each source's\n"
    "best-matching unit (BMU) attached (`workshop_agn_edfn.fits`). Both are derived from queries\n"
    "against the same cloud-hosted Q1 catalogs you'll touch live in Task 2.\n"
    "\n"
    "**Scope.** EDF-N AGN candidates from Euclid Collaboration: Matamoro Zatarain et al. (2025),\n"
    "in the redshift window `1.30 < z < 1.80` (H$\\alpha$ in the NISP red grism), with the same\n"
    "photometric gates (per-band S/N > 5, AB in [12, 30], MER 7-band) as the full tutorial.\n"
))

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
CELLS.append(md(
    "## 0. Setup\n"
    "\n"
    "Load libraries and the pre-built local artifacts (SOM + projected catalog).\n"
    "No cloud reads yet - those start in Task 2.\n"
))

CELLS.append(code(
    "# Standard library and third-party imports.\n"
    "import os\n"
    "import json\n"
    "import pickle\n"
    "import warnings\n"
    "\n"
    "import numpy as np\n"
    "import matplotlib.pyplot as plt\n"
    "from astropy.table import Table\n"
    "from astropy.utils.metadata import MergeConflictWarning\n"
    "from astroquery.ipac.irsa import Irsa\n"
    "\n"
    "# Project helpers (data IO, SOM utilities, plot overlays, spectrum/cutout fetch).\n"
    "from agn_tutorial_utils import (\n"
    "    AGN_OVERLAY_COLOR, COLOR_NAMES, METHOD_CONTOUR_COLOR, SPE_LINES_TABLE,\n"
    "    SPECTRA_ASSOC_TABLE, SPECTRA_BUCKET,\n"
    "    bmu_density, continuum_normalized_rest_spectrum, fill_nan_nearest,\n"
    "    get_Q1_mer_cutout, get_Q1_sir_spectra, jitter_bmu, mark_halpha_complex,\n"
    "    median_map, overlay, overlay_method_excess,\n"
    ")\n"
    "\n"
    "warnings.filterwarnings('ignore', category=MergeConflictWarning)\n"
    "%matplotlib inline\n"
    "\n"
    "CACHE_DIR = 'data/cache'\n"
    "print('Setup complete.')\n"
))

CELLS.append(code(
    "# Load the pre-built SOM and the projected EDF-N catalogs.\n"
    "with open(os.path.join(CACHE_DIR, 'workshop_som.pkl'), 'rb') as fh:\n"
    "    pack = pickle.load(fh)\n"
    "som = pack['som']\n"
    "mu, sigma = pack['mu'], pack['sigma']\n"
    "qe_thresh = pack['qe_thresh']\n"
    "MSZ = tuple(pack['msz'])\n"
    "Z_WIN = tuple(pack['z_window'])\n"
    "\n"
    "agn = Table.read(os.path.join(CACHE_DIR, 'workshop_agn_edfn.fits'))\n"
    "train = Table.read(os.path.join(CACHE_DIR, 'workshop_train_edfn.fits'))\n"
    "with open(os.path.join(CACHE_DIR, 'workshop_picks.json')) as fh:\n"
    "    picks = json.load(fh)\n"
    "\n"
    "# Inlier subset: AGN whose colors land inside the trained galaxy manifold.\n"
    "agn_in = agn[agn['qe_inlier']]\n"
    "bmu_in = np.column_stack([agn_in['bmu_row'], agn_in['bmu_col']]).astype(int)\n"
    "bmu_train = np.column_stack([train['bmu_row'], train['bmu_col']]).astype(int)\n"
    "\n"
    "print(f'SOM grid: {MSZ[0]} x {MSZ[1]}')\n"
    "print(f'Redshift window: {Z_WIN[0]} < z < {Z_WIN[1]}')\n"
    "print(f'Training galaxies on the SOM:       {len(train):,}')\n"
    "print(f'EDF-N AGN projected:                {len(agn):,}')\n"
    "print(f'  inliers (QE <= {qe_thresh:.2f}):    {len(agn_in):,}')\n"
    "print(f'  with spec-z:                      {int(np.sum(agn_in[\"z_source\"] == \"spec\")):,}')\n"
    "print(f'Three pre-picked candidates loaded: {[p[\"object_id_euclid\"] for p in picks]}')\n"
))

# ---------------------------------------------------------------------------
# Warm-up plot - SOM background
# ---------------------------------------------------------------------------
CELLS.append(md(
    "## 1. The Q1 color SOM at a glance\n"
    "\n"
    "The SOM was trained on Q1 galaxies in `1.3 < z < 1.8` using six adjacent colors\n"
    "(`g-r`, `r-i`, `i-z`, `z-Y`, `Y-J`, `J-H`). The map below is colored by the\n"
    "median observed-frame H-band magnitude per cell - a simple host-brightness\n"
    "proxy (darker = brighter). The green overlay on the right shows where the\n"
    "**union** of AGN candidates sits on this manifold.\n"
))

CELLS.append(code(
    "# Median observed-frame mH per cell (the manifold's brightness background).\n"
    "med_mH = median_map(train['mH_obs'], bmu_train, msz=MSZ)\n"
    "mH_bg = fill_nan_nearest(med_mH)\n"
    "mH_lo, mH_hi = np.nanpercentile(med_mH, [5, 95])\n"
    "\n"
    "fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))\n"
    "for ax in axes:\n"
    "    im_bg = ax.imshow(mH_bg, origin='lower', cmap='cividis_r', vmin=mH_lo, vmax=mH_hi)\n"
    "    ax.set_xticks([]); ax.set_yticks([])\n"
    "fig.colorbar(im_bg, ax=axes[0], fraction=0.046, pad=0.04).set_label('median $m_H$ per cell', fontsize=10)\n"
    "axes[0].set_title(f'Q1 color SOM colored by median $m_H$\\n(n={len(train):,} galaxies)', fontsize=11)\n"
    "\n"
    "# Right panel: union AGN density overlay.\n"
    "H_union = bmu_density(bmu_in)\n"
    "im_d, vlo, vhi = overlay(axes[1], H_union)\n"
    "if im_d is not None:\n"
    "    fig.colorbar(im_d, ax=axes[1], fraction=0.046, pad=0.04).set_label('smoothed AGN / cell (log)', fontsize=10)\n"
    "axes[1].set_title(f'Union AGN density on the SOM\\n(n={len(agn_in):,} inlier AGN)', fontsize=11)\n"
    "plt.tight_layout(); plt.show()\n"
))

# ---------------------------------------------------------------------------
# Task 1
# ---------------------------------------------------------------------------
CELLS.append(md(
    "## Task 1 - Where does *your* selection method live? (~10 min)\n"
    "\n"
    "The union map above mixes every AGN selection technique into one footprint.\n"
    "Different methods sample different physics: WISE colors trace hot dust,\n"
    "Euclid NIR colors trace mid-IR-like features in extended bandpasses, Gaia DR3\n"
    "selects bright quasars by parallax/photometry, etc. On a common manifold we\n"
    "can ask: *do these methods land in the same place, or different ones?*\n"
    "\n"
    "Available flags in this EDF-N sample (each is a 0/1 column on `agn_in`):\n"
    "\n"
    "| Flag column | What it selects |\n"
    "| --- | --- |\n"
    "| `R90_agn_candidate` | WISE mid-IR color, 90% reliability (Assef+18) |\n"
    "| `C75_agn_candidate` | WISE mid-IR color, 75% completeness (Assef+18) |\n"
    "| `JH_IeY_qso_candidate` | Euclid NIR two-color box (Bisigello+24) |\n"
    "| `IeH_gz_qso_candidate` | Euclid + optical color box (Bisigello+24) |\n"
    "| `B24a_qso_candidate` | Optical + NIR color cut (Bisigello+24 a) |\n"
    "| `B24b_qso_candidate` | Optical + NIR color cut (Bisigello+24 b) |\n"
    "| `GDR3_qso_candidate` | Gaia DR3 quasar candidate |\n"
    "| `PRF_qso_candidate` | Probabilistic Random Forest |\n"
    "| `AGN_sed_candidate` | SED-template based (EDF-N only) |\n"
    "| `DESI_broadline_qso_candidate` | DESI spec-confirmed broad-line QSO |\n"
    "| `DESI_niibpt_agn_candidate` | DESI [NII]/H$\\alpha$ BPT AGN |\n"
    "\n"
    "**Your task:** pick *one* flag in the cell below. The cell will plot two panels:\n"
    "left = where that method's AGN sit on the SOM; right = where they are *enriched*\n"
    "relative to the union AGN distribution (green = preferred regions for that method).\n"
    "\n"
    "Compare with your neighbor: did you pick the same place on the manifold?\n"
))

CELLS.append(code(
    "# === TODO 1: pick one flag column from the table above ====================\n"
    "method = 'R90_agn_candidate'      # <-- change this to any flag column name\n"
    "# ==========================================================================\n"
    "\n"
    "assert method in agn_in.colnames, f'{method!r} is not a column in agn_in.'\n"
    "method_mask = np.asarray(agn_in[method] == 1)\n"
    "n_method = int(method_mask.sum())\n"
    "print(f'AGN flagged by {method}: {n_method:,} (out of {len(agn_in):,} inliers)')\n"
    "\n"
    "fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))\n"
    "for ax in axes:\n"
    "    ax.imshow(mH_bg, origin='lower', cmap='cividis_r', vmin=mH_lo, vmax=mH_hi)\n"
    "    ax.set_xticks([]); ax.set_yticks([])\n"
    "\n"
    "# Left: BMUs and smoothed density of the selected method.\n"
    "H_method = bmu_density(bmu_in[method_mask]) if n_method else None\n"
    "if H_method is not None:\n"
    "    im_d, _, _ = overlay(axes[0], H_method)\n"
    "    if im_d is not None:\n"
    "        fig.colorbar(im_d, ax=axes[0], fraction=0.046, pad=0.04).set_label('smoothed AGN / cell (log)', fontsize=10)\n"
    "if 0 < n_method <= 60:\n"
    "    x_pt, y_pt = jitter_bmu(bmu_in[method_mask], seed=42)\n"
    "    axes[0].scatter(x_pt, y_pt, marker='x', s=30, linewidths=0.9,\n"
    "                    color='white', alpha=0.9, zorder=5)\n"
    "    axes[0].scatter(x_pt, y_pt, marker='x', s=30, linewidths=0.55,\n"
    "                    color=AGN_OVERLAY_COLOR, alpha=0.9, zorder=6)\n"
    "axes[0].set_title(f'{method} density\\n(n={n_method:,})', fontsize=11)\n"
    "\n"
    "# Right: method-excess over the union distribution.\n"
    "im_x = overlay_method_excess(axes[1], H_method, H_union)\n"
    "if im_x is not None:\n"
    "    fig.colorbar(im_x, ax=axes[1], fraction=0.046, pad=0.04).set_label('excess (method / union)', fontsize=10)\n"
    "axes[1].set_title(f'{method} excess over union\\n(green = method-enriched cells)', fontsize=11)\n"
    "plt.tight_layout(); plt.show()\n"
))

CELLS.append(md(
    "**Discussion prompt.** Where on the SOM did your method land? Bright-host cells (light)\n"
    "or faint-host cells (dark)? Is the enrichment spread out across the manifold or\n"
    "concentrated in one corner? Compare with someone who picked a different flag - the\n"
    "*excess* maps usually disagree, even when the *density* maps look similar.\n"
))

# ---------------------------------------------------------------------------
# Task 2
# ---------------------------------------------------------------------------
CELLS.append(md(
    "## Task 2 - Project a candidate, then verify it with three live cloud reads (~15 min)\n"
    "\n"
    "We pre-selected three AGN candidates that each land in a distinct SOM neighborhood\n"
    "and are selected by a different combination of methods. Pick one, project it onto\n"
    "the SOM, then pull three independent measurements **live from IRSA** to check\n"
    "whether the colors-and-manifold prediction agrees with the actual photons:\n"
    "\n"
    "- **2b** SIR 1D spectrum - direct S3 FITS read from `s3://nasa-irsa-euclid-q1/...`\n"
    "- **2c** MER multi-band imaging - server-side cutout from IRSA IBE\n"
    "- **2d** SPE H$\\alpha$ line measurement - one-row TAP query against the SPE line catalog\n"
))

CELLS.append(code(
    "# Pre-picked candidates, with the flags that selected them and their SOM cell.\n"
    "for i, p in enumerate(picks):\n"
    "    flags_short = ', '.join(p['flags_set'][:4])\n"
    "    if len(p['flags_set']) > 4:\n"
    "        flags_short += f' ... (+{len(p[\"flags_set\"]) - 4} more)'\n"
    "    print(f'[{i}] object_id = {p[\"object_id_euclid\"]}')\n"
    "    print(f'    RA, Dec   = {p[\"ra\"]:.5f}, {p[\"dec\"]:.5f}')\n"
    "    print(f'    z = {p[\"z_final\"]:.3f}  ({p[\"z_source\"]})')\n"
    "    print(f'    SOM cell  = ({p[\"bmu_row\"]}, {p[\"bmu_col\"]})')\n"
    "    print(f'    flagged by: {flags_short}')\n"
    "    print()\n"
))

CELLS.append(code(
    "# === TODO 2: pick a candidate =============================================\n"
    "# Recommended first pass: one of the three curated picks above. Each has a\n"
    "# Q1 SIR spectrum, spec-z, and a distinct SOM neighborhood.\n"
    "my_id = picks[0]['object_id_euclid']     # also try picks[1] or picks[2]\n"
    "\n"
    "# Or pick any AGN from the inlier sample. SIR coverage in EDF-N is very\n"
    "# high (essentially every inlier AGN has a spectrum), but the cells below\n"
    "# handle the missing case gracefully if you happen to hit one. Uncomment:\n"
    "#   my_id = int(agn_in['object_id_euclid'][0])              # first inlier\n"
    "#   my_id = int(agn_in[agn_in['z_source']=='spec']['object_id_euclid'][0])  # first spec-z\n"
    "# ==========================================================================\n"
    "\n"
    "row = agn[agn['object_id_euclid'] == my_id][0]\n"
    "ra, dec = float(row['ra_euclid']), float(row['dec_euclid'])\n"
    "z = float(row['z_final'])\n"
    "print(f'Chose object_id {my_id}')\n"
    "print(f'  RA, Dec = {ra:.6f}, {dec:.6f}')\n"
    "print(f'  z = {z:.3f} ({row[\"z_source\"]})')\n"
))

CELLS.append(md(
    "### 2a. Project the candidate onto the SOM (local compute)\n"
    "\n"
    "To project a source onto the SOM you need its six adjacent colors normalized the\n"
    "same way as the training sample (z-score with the training `mu`, `sigma`), then\n"
    "the SOM's `winner()` method returns the best-matching unit (BMU). This step is\n"
    "pure local arithmetic - no cloud reads.\n"
))

CELLS.append(code(
    "# === TODO 3: normalize the candidate's colors and find its BMU ============\n"
    "colors = np.array([row[c] for c in COLOR_NAMES])   # six adjacent colors (mag)\n"
    "print('Colors:', dict(zip(COLOR_NAMES, np.round(colors, 3))))\n"
    "\n"
    "v_norm = (colors - mu) / sigma                      # z-score with training stats\n"
    "bmu = som.winner(v_norm)                            # (row, col) on the SOM grid\n"
    "# ==========================================================================\n"
    "print(f'Projected BMU: {bmu}  (cached value was ({int(row[\"bmu_row\"])}, {int(row[\"bmu_col\"])}))')\n"
    "\n"
    "# Plot the BMU on top of the union AGN density background.\n"
    "fig, ax = plt.subplots(figsize=(6.4, 5.4))\n"
    "ax.imshow(mH_bg, origin='lower', cmap='cividis_r', vmin=mH_lo, vmax=mH_hi)\n"
    "if H_union is not None:\n"
    "    overlay(ax, H_union)\n"
    "ax.scatter(bmu[1], bmu[0], marker='*', s=320, edgecolors='black', linewidths=1.3,\n"
    "           facecolor='gold', zorder=10, label=f'your AGN  ({bmu[0]}, {bmu[1]})')\n"
    "ax.set_xticks([]); ax.set_yticks([])\n"
    "ax.legend(loc='upper right', frameon=True, fontsize=9)\n"
    "ax.set_title(f'object_id {my_id} on the SOM\\n(z = {z:.3f}, {row[\"z_source\"]}-z)', fontsize=11)\n"
    "plt.tight_layout(); plt.show()\n"
))

# --- 2b: live SIR 1D spectrum direct from S3 -------------------------------
CELLS.append(md(
    "### 2b. Pull the SIR spectrum live from S3\n"
    "\n"
    "Two cloud calls happen here. First, a TAP query against IRSA's spectrum-association\n"
    "table to find which FITS file in the public bucket contains this object's spectrum\n"
    "(a sub-second SQL lookup). Second, an `astropy.io.fits.open(s3_uri, ...)` against\n"
    "`s3://nasa-irsa-euclid-q1/...` with `fsspec` - that's a **direct anonymous S3 read**\n"
    "of the FITS file, no download.\n"
))

CELLS.append(code(
    "# Step 1: TAP query for the spectrum-file association.\n"
    "adql_assoc = (\n"
    "    'SELECT objectid, path, hdu '\n"
    "    f'FROM {SPECTRA_ASSOC_TABLE} '\n"
    "    f'WHERE objectid = {my_id}'\n"
    ")\n"
    "print('TAP query:')\n"
    "print(adql_assoc)\n"
    "assoc = Irsa.query_tap(adql_assoc).to_table()\n"
    "print(f'Returned {len(assoc)} row(s).')\n"
    "\n"
    "have_spectrum = len(assoc) > 0\n"
    "if not have_spectrum:\n"
    "    print('  -- no Q1 SIR spectrum recorded for this object_id.')\n"
    "    print('     Imaging (2c) still works; pick a different object_id if you')\n"
    "    print('     want to see a spectrum here.')\n"
    "else:\n"
    "    for r in assoc:\n"
    "        print(f'  hdu={int(r[\"hdu\"]):>3d}  path={r[\"path\"]}')\n"
    "    # Step 2: build the S3 URI for the spectrum FITS and open it directly.\n"
    "    from agn_tutorial_utils import association_to_s3_uri\n"
    "    first = assoc[0]\n"
    "    s3_uri = association_to_s3_uri(first['path'])\n"
    "    print(f'\\nReading from: {s3_uri}')\n"
    "    print(f'  HDU index in that file: {int(first[\"hdu\"])}')\n"
))

CELLS.append(code(
    "# Skip cleanly if the source has no SIR association; otherwise read the FITS\n"
    "# from S3 and plot. The helper does the same fits.open(s3_uri, fsspec_kwargs=...)\n"
    "# under the hood, with a local .npz cache for repeat calls on the same id.\n"
    "if not have_spectrum:\n"
    "    print('No SIR spectrum to read. Move on to cell 2c (imaging).')\n"
    "else:\n"
    "    one_row = agn[agn['object_id_euclid'] == my_id]\n"
    "    spectra = get_Q1_sir_spectra(one_row, n=1, cache_check=True)\n"
    "    if my_id not in spectra:\n"
    "        print('Association returned a row but the FITS HDU was unreadable. Skipping.')\n"
    "    else:\n"
    "        spec = spectra[my_id]\n"
    "        print(f'wavelength samples: {len(spec[\"wave\"]):,}')\n"
    "        print(f'observed range:     {spec[\"wave\"][0]:.0f} - {spec[\"wave\"][-1]:.0f}')\n"
    "\n"
    "        rest, y = continuum_normalized_rest_spectrum(spec, z)\n"
    "        fig, (axw, axz) = plt.subplots(1, 2, figsize=(13.0, 4.2),\n"
    "                                        gridspec_kw={'width_ratios': [1.6, 1.0]})\n"
    "        axw.plot(rest, y, color='#064d2c', lw=0.8)\n"
    "        axw.axhline(0, color='0.5', lw=0.6, alpha=0.5)\n"
    "        axw.axvspan(6400, 6760, color='#82c77a', alpha=0.15, lw=0)\n"
    "        mark_halpha_complex(axw, full_range=True)\n"
    "        axw.set_xlim(np.nanmin(rest), np.nanmax(rest))\n"
    "        axw.set_xlabel('rest wavelength (Angstrom)')\n"
    "        axw.set_ylabel('continuum-subtracted flux / scale')\n"
    "        axw.set_title(f'SIR 1D spectrum  -  object_id {my_id}, z={z:.3f}')\n"
    "        zoom = (rest > 6400) & (rest < 6760)\n"
    "        axz.plot(rest[zoom], y[zoom], color='#064d2c', lw=1.1)\n"
    "        axz.axhline(0, color='0.5', lw=0.6, alpha=0.5)\n"
    "        mark_halpha_complex(axz, full_range=True)\n"
    "        axz.set_xlim(6400, 6760)\n"
    "        axz.set_xlabel('rest wavelength (Angstrom)')\n"
    "        axz.set_title(r'H$\\alpha$ zoom (dashed=H$\\alpha$, dotted=[NII], [SII])')\n"
    "        plt.tight_layout(); plt.show()\n"
))

# --- 2c: live MER imaging cutouts via IBE ----------------------------------
CELLS.append(md(
    "### 2c. Pull multi-band MER imaging cutouts live from IRSA IBE\n"
    "\n"
    "Now an imaging look at the same source. We use IRSA's IBE cutout service:\n"
    "an SIA query returns the URL of the full MER mosaic FITS that covers our position;\n"
    "we append `?center=ra,dec&size=Xarcsec` to ask the server to return only the\n"
    "small region we want. The server-side cutout comes back as a ~40 KB gzipped FITS\n"
    "in a fraction of a second - no need to stream the full 1.4 GB mosaic.\n"
    "\n"
    "Below we fetch VIS (high-resolution optical) and three NISP NIR bands (Y, J, H).\n"
))

CELLS.append(code(
    "# Fetch four bands - the helper prints each cutout URL so you can see the cloud calls.\n"
    "from agn_tutorial_utils import _mer_sia_lookup\n"
    "import time\n"
    "\n"
    "t0 = time.time()\n"
    "sia_rows = _mer_sia_lookup(ra, dec)\n"
    "print(f'SIA query: {time.time()-t0:.2f} s, {len(sia_rows)} mosaic rows at this position')\n"
    "print('IBE cutout URLs being fetched:')\n"
    "\n"
    "BANDS = ['VIS', 'Y', 'J', 'H']\n"
    "cutouts = {}\n"
    "t0 = time.time()\n"
    "for b in BANDS:\n"
    "    cutouts[b] = get_Q1_mer_cutout(ra, dec, band=b, size_arcsec=8.0, sia_rows=sia_rows, verbose=True)\n"
    "print(f'Fetched {len(BANDS)} cutouts in {time.time()-t0:.2f} s')\n"
    "\n"
    "# Print one S3 URI to make the cloud location explicit.\n"
    "print(f'\\n(For reference, the underlying S3 path of the VIS mosaic was:')\n"
    "print(f' {cutouts[\"VIS\"][\"s3_uri\"]})')\n"
))

CELLS.append(code(
    "# Display the four bands side-by-side with consistent arcsinh stretch per panel.\n"
    "from matplotlib.colors import AsinhNorm\n"
    "\n"
    "fig, axes = plt.subplots(1, len(BANDS), figsize=(3.0 * len(BANDS) + 0.5, 3.2))\n"
    "for ax, b in zip(axes, BANDS):\n"
    "    data = cutouts[b]['data']\n"
    "    vmin = float(np.nanpercentile(data, 5))\n"
    "    vmax = float(np.nanpercentile(data, 99.5))\n"
    "    span = max(vmax - vmin, 1e-6)\n"
    "    ax.imshow(data, origin='lower', cmap='gray',\n"
    "              norm=AsinhNorm(linear_width=0.1 * span, vmin=vmin, vmax=vmax))\n"
    "    ax.set_xticks([]); ax.set_yticks([])\n"
    "    ax.set_title(f'{b}  ({data.shape[1]} x {data.shape[0]} pix)', fontsize=10)\n"
    "    # Center crosshair on source position.\n"
    "    yc, xc = data.shape[0] / 2, data.shape[1] / 2\n"
    "    ax.plot(xc, yc, marker='+', color='red', markersize=14, markeredgewidth=1.5)\n"
    "fig.suptitle(f'object_id {my_id}  -  live MER cutouts streamed from IRSA IBE', fontsize=11)\n"
    "plt.tight_layout(); plt.show()\n"
))

# --- 2d: live SPE Halpha catalog lookup -----------------------------------
CELLS.append(md(
    "### 2d. Look up the catalog SPE H$\\alpha$ measurement (live TAP)\n"
    "\n"
    "Finally a quick cross-check against Euclid's own automated line-feature catalog.\n"
    "SPE measures Gaussian fits to spectral lines and stores flux, EW, and S/N per line.\n"
    "We run a one-row ADQL query for this object's H$\\alpha$ row(s) and compare with\n"
    "what we just saw in the 1D spectrum.\n"
))

CELLS.append(code(
    "# === Live TAP query - the SPE line-feature catalog is not pre-cached =====\n"
    "adql_spe = (\n"
    "    'SELECT spe_rank, spe_line_central_wl_gf, spe_line_flux_gf, '\n"
    "    'spe_line_flux_err_gf, spe_line_snr_gf, spe_line_ew_gf '\n"
    "    f'FROM {SPE_LINES_TABLE} '\n"
    "    f'WHERE object_id = {my_id} '\n"
    '    \"AND spe_line_name = \'Halpha\'\"\n'
    ")\n"
    "print('Live TAP query:')\n"
    "print(adql_spe)\n"
    "\n"
    "result = Irsa.query_tap(adql_spe).to_table()\n"
    "print(f'\\nReturned {len(result)} row(s):')\n"
    "if len(result) == 0:\n"
    "    print('  (no Halpha measurement in the SPE catalog for this source.)')\n"
    "else:\n"
    "    result.pprint(max_lines=10, max_width=120)\n"
    "    print()\n"
    "    best = result[np.argmax(np.asarray(result['spe_line_snr_gf'], dtype=float))]\n"
    "    print(f'Best-S/N row:  central wavelength = {float(best[\"spe_line_central_wl_gf\"]):.1f} Angstrom (observed)')\n"
    "    print(f'               flux  = {float(best[\"spe_line_flux_gf\"]):.3e} erg/s/cm2')\n"
    "    print(f'               S/N   = {float(best[\"spe_line_snr_gf\"]):.2f}')\n"
    "    print(f'               |EW|  = {abs(float(best[\"spe_line_ew_gf\"])):.1f} Angstrom')\n"
    "    print(f'Predicted Halpha at z={z:.3f}: observed wavelength = '\n"
    "          f'{6562.8 * (1+z):.1f} Angstrom (rest 6562.8 A)')\n"
))

CELLS.append(md(
    "**Three independent cloud reads, one consistent answer.** The 1D spectrum shows the\n"
    "line, the imaging shows the host, and the SPE catalog quantifies the line - all\n"
    "fetched live from `irsa.ipac.caltech.edu` or `nasa-irsa-euclid-q1` in this notebook,\n"
    "with no preloaded data products beyond the SOM artifacts.\n"
    "\n"
    "Try changing `my_id` to one of the other curated picks, or to any `object_id_euclid`\n"
    "from `agn_in`, and rerunning Task 2. Different SOM neighborhoods often correspond to\n"
    "different AGN flavors - compact point source vs extended host, broad-line vs narrow,\n"
    "bright vs faint - which you can read off directly from the imaging cell.\n"
))

# ---------------------------------------------------------------------------
# Wrap
# ---------------------------------------------------------------------------
CELLS.append(md(
    "## Wrap-up\n"
    "\n"
    "In ~30 minutes you have:\n"
    "\n"
    "1. Loaded a pre-trained Q1 color SOM and an EDF-N AGN catalog projected onto it.\n"
    "2. Picked one AGN selection method and visualized where it is enriched on the manifold.\n"
    "3. Picked an individual candidate and verified it with **three independent live cloud reads**:\n"
    "   - a 1D SIR spectrum read directly from `s3://nasa-irsa-euclid-q1/...`,\n"
    "   - multi-band MER imaging cutouts streamed from IRSA's IBE server-side cutout service,\n"
    "   - a one-row SPE H$\\alpha$ measurement from a live TAP query.\n"
    "\n"
    "**The cloud point.** None of those data products were downloaded in advance. The\n"
    "same code, run from a Fornax JupyterLab session in the same AWS region, would\n"
    "complete in a fraction of the time you saw here (the SIA + IBE round-trips dominate\n"
    "from a laptop). It scales unchanged to the full Q1 release (~30M sources, ~150 deg2)\n"
    "and onward to DR1 with no rewrite.\n"
    "\n"
    "**Where to go next.**\n"
    "\n"
    "- Full pipeline (queries, training, all 3 fields, north/south overlays, SPE H$\\alpha$\n"
    "  diagnostics): [`Euclid_AGN_Tutorial.ipynb`](Euclid_AGN_Tutorial.ipynb).\n"
    "- IRSA's reference Euclid Q1 cloud-access tutorial:\n"
    "  https://caltech-ipac.github.io/irsa-tutorials/euclid-cloud-access/\n"
    "- Cluster-finding example: [IRSA Euclid clusters tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-clusters-tutorial/).\n"
    "- SPE line-feature catalog: [IRSA Euclid SPE catalog tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-spe-catalog/).\n"
    "\n"
    "**Acknowledgments.** AGN candidate catalogs from Euclid Collaboration: Matamoro Zatarain\n"
    "et al. (2025); Q1 data products distributed by NASA/IPAC IRSA. The SOM methodology\n"
    "extends Sanjaripour et al. (2024) and a companion Euclid Q1 paper (in prep).\n"
))


def main():
    nb = {
        'cells': CELLS,
        'metadata': {
            'kernelspec': {
                'display_name': 'Python 3',
                'language': 'python',
                'name': 'python3',
            },
            'language_info': {
                'name': 'python',
                'pygments_lexer': 'ipython3',
            },
        },
        'nbformat': 4,
        'nbformat_minor': 5,
    }
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'Euclid_AGN_Workshop.ipynb')
    with open(out_path, 'w') as fh:
        json.dump(nb, fh, indent=1)
    print(f'Wrote {out_path}  ({len(CELLS)} cells)')


if __name__ == '__main__':
    main()
