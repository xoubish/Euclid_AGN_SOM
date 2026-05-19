"""Generate the workshop assignment notebook from the filled-in solution.

Reads `workshop/Euclid_AGN_Workshop.ipynb` (solution, source of truth) and
writes `workshop/Euclid_AGN_Workshop_assignment.ipynb` with the bodies of the
TODO blocks replaced by hints. Hints are limited to things participants can't
guess from context: IRSA table names, column-name conventions (which differ
between tables!), and helper-function signatures.

Re-run after editing the solution so the assignment stays in sync. The script
errors out if any solution snippet no longer matches, so drift is visible.
"""
import json
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, 'workshop', 'Euclid_AGN_Workshop.ipynb')
DST = os.path.join(ROOT, 'workshop', 'Euclid_AGN_Workshop_assignment.ipynb')


# (find_in_solution, replace_in_assignment) pairs. Each find_text must be
# unique to one cell.
TODO_REPLACEMENTS = [
    # ----------------------------------------------------------------------
    # Task 1: live TAP query for MER photometry
    # ----------------------------------------------------------------------
    (
        "    phot = batched_query('euclid_q1_mer_catalogue', FLUX_COLS, ids,\n"
        "                         id_col='object_id', desc='MER phot')",
        "    # Hints:\n"
        "    #   - MER catalog table name: 'euclid_q1_mer_catalogue'\n"
        "    #   - In this table the ID column is 'object_id' (NOT 'object_id_euclid')\n"
        "    #   - Use FLUX_COLS (defined above) for the columns to SELECT\n"
        "    #   - batched_query(table, cols, ids, id_col=..., desc=...) is imported\n"
        "    #     from workshop_utils and returns an astropy Table\n"
        "    phot = ...   # <-- batched_query(...)",
    ),

    # ----------------------------------------------------------------------
    # Task 2: pick a flag column from the table above
    # ----------------------------------------------------------------------
    (
        "method = 'R90_agn_candidate'      # <-- change this to any flag column name",
        "# Pick any flag name (as a string) from the table above. Try a few - some\n"
        "# methods land in the same neighborhood on the SOM, others don't.\n"
        "method = ...   # <-- e.g. 'R90_agn_candidate'",
    ),

    # ----------------------------------------------------------------------
    # Task 3: pick a candidate (random draw from spec-z inliers)
    # ----------------------------------------------------------------------
    (
        "# Change the seed (or set seed=None) to draw a different AGN.\n"
        "rng = np.random.default_rng(seed=42)\n"
        "row = specz[int(rng.integers(len(specz)))]",
        "# Hint: draw a random row from `specz`. Use np.random.default_rng() so your\n"
        "# draw is reproducible across kernel restarts.\n"
        "rng = ...        # <-- np.random.default_rng(seed=42) (or any seed)\n"
        "row = ...        # <-- pick one row using rng.integers(len(specz))",
    ),

    # ----------------------------------------------------------------------
    # Task 3a: project the candidate onto the SOM
    # ----------------------------------------------------------------------
    (
        "v_norm = (colors - mu) / sigma                      # z-score with training stats\n"
        "bmu = som.winner(v_norm)                            # (row, col) on the SOM grid",
        "# Hint: z-score the colors using the training `mu` and `sigma`, then pass\n"
        "# the result to `som.winner()`. The return is a (row, col) tuple.\n"
        "v_norm = ...    # <-- (colors - mu) / sigma\n"
        "bmu = ...       # <-- som.winner(v_norm)",
    ),

    # ----------------------------------------------------------------------
    # Task 3b: ADQL for the SIR spectrum-file association
    # ----------------------------------------------------------------------
    (
        "adql_assoc = (\n"
        "    'SELECT objectid, path, hdu '\n"
        "    f'FROM {SPECTRA_ASSOC_TABLE} '\n"
        "    f'WHERE objectid = {my_id}'\n"
        ")",
        "# Hints:\n"
        "#   - Table name is stored in the constant SPECTRA_ASSOC_TABLE (imported above);\n"
        "#     reference the constant rather than hard-coding the string.\n"
        "#   - This table's ID column is `objectid` (one word, no underscore).\n"
        "#     Naming conventions vary by IRSA table - watch for this in 3d.\n"
        "#   - Columns you need from each row: objectid, path, hdu.\n"
        "adql_assoc = ...   # <-- f-string SELECT ... FROM ... WHERE objectid = {my_id}",
    ),

    # ----------------------------------------------------------------------
    # Task 3c: SIA lookup + per-band IBE cutouts
    # ----------------------------------------------------------------------
    (
        "t0 = time.time()\n"
        "sia_rows = _mer_sia_lookup(ra, dec)\n"
        "print(f'SIA query: {time.time()-t0:.2f} s, {len(sia_rows)} mosaic rows at this position')",
        "# Hint: _mer_sia_lookup(ra_deg, dec_deg) issues an SIA query against IRSA\n"
        "# and returns the list of MER mosaic FITS rows that cover that position.\n"
        "t0 = time.time()\n"
        "sia_rows = ...   # <-- _mer_sia_lookup(ra, dec)\n"
        "print(f'SIA query: {time.time()-t0:.2f} s, {len(sia_rows)} mosaic rows at this position')",
    ),
    (
        "BANDS = ['VIS', 'Y', 'J', 'H']\n"
        "cutouts = {}\n"
        "t0 = time.time()\n"
        "for b in BANDS:\n"
        "    cutouts[b] = get_Q1_mer_cutout(ra, dec, band=b, size_arcsec=8.0, sia_rows=sia_rows, verbose=True)",
        "# Hint: get_Q1_mer_cutout(ra, dec, band=..., size_arcsec=..., sia_rows=...,\n"
        "# verbose=True) returns a dict with the cutout image data and metadata,\n"
        "# fetched live from IRSA IBE. Pass the same sia_rows each call to reuse the\n"
        "# SIA result (otherwise each call re-queries SIA).\n"
        "BANDS = ['VIS', 'Y', 'J', 'H']\n"
        "cutouts = {}\n"
        "t0 = time.time()\n"
        "for b in BANDS:\n"
        "    cutouts[b] = ...   # <-- get_Q1_mer_cutout(ra, dec, band=b, size_arcsec=8.0, sia_rows=sia_rows, verbose=True)",
    ),

    # ----------------------------------------------------------------------
    # Task 3d: ADQL for the SPE H-alpha line measurement
    # ----------------------------------------------------------------------
    (
        "adql_spe = (\n"
        "    'SELECT spe_rank, spe_line_central_wl_gf, spe_line_flux_gf, '\n"
        "    'spe_line_flux_err_gf, spe_line_snr_gf, spe_line_ew_gf '\n"
        "    f'FROM {SPE_LINES_TABLE} '\n"
        "    f'WHERE object_id = {my_id} '\n"
        "    \"AND spe_line_name = 'Halpha'\"\n"
        ")",
        "# Hints:\n"
        "#   - Table name is stored in the constant SPE_LINES_TABLE (imported above).\n"
        "#   - This table's ID column is `object_id` (with underscore) -- DIFFERENT from\n"
        "#     the spectrum-association table in 3b which used `objectid`. IRSA's\n"
        "#     naming conventions vary between catalogs and there's no shortcut.\n"
        "#   - You need TWO WHERE conditions: object_id = my_id AND spe_line_name = 'Halpha'.\n"
        "#   - Columns to SELECT: spe_rank, spe_line_central_wl_gf, spe_line_flux_gf,\n"
        "#     spe_line_flux_err_gf, spe_line_snr_gf, spe_line_ew_gf.\n"
        "adql_spe = ...   # <-- f-string with SELECT/FROM/WHERE/AND",
    ),
]


def main():
    with open(SRC) as fh:
        nb = json.load(fh)

    n_applied = 0
    for cell in nb['cells']:
        if cell['cell_type'] == 'code':
            cell['outputs'] = []
            cell['execution_count'] = None
        src = cell['source']
        if isinstance(src, list):
            src = ''.join(src)
        for find, replace in TODO_REPLACEMENTS:
            if find in src:
                src = src.replace(find, replace)
                n_applied += 1
        cell['source'] = src

    if n_applied != len(TODO_REPLACEMENTS):
        raise SystemExit(
            f'Expected {len(TODO_REPLACEMENTS)} TODO replacements, applied {n_applied}. '
            f'The solution notebook may have drifted - check the find_text strings.'
        )

    with open(DST, 'w') as fh:
        json.dump(nb, fh, indent=1)
    print(f'wrote {DST} ({len(nb["cells"])} cells, {n_applied} TODOs blanked)')


if __name__ == '__main__':
    main()
