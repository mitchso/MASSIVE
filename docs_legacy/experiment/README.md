## Experiment

**File:** `experiment.py`

Manages a full plate of `Sample` objects loaded from a folder of text files. Provides plate-level visualization and data export.

### Constructor

```python
Experiment(data_folder=None, noise_cutoff=0.95)
```

Pass `data_folder` as the path to a directory of `.txt` spectrum files. `noise_cutoff` sets the noise threshold applied to every `Sample`.

**Key attributes after construction:**

- `samples` — dict of `{well_id: Sample}` for chip=0 positions
- `calibrant_spots` — dict of `{well_id: Sample}` for chip=1 positions
- `exclude` — list of well IDs to skip in downstream analysis (set manually)

---

### Data Loading

#### `data_folder_to_sample_dict() → (dict, dict)`
Parses all `.txt` files in `data_folder`, constructs `Sample` objects from each, and returns two dicts: one for sample wells and one for calibrant spots. Filenames are parsed to extract chip, well, and replicate information.

#### `reinitiate_samples(samples) → dict`
Re-creates `Sample` objects from an existing list of samples, applying the current `Experiment.noise_cutoff`. Useful when loading pre-existing samples under a different noise setting.

---

### Visualization

#### `show_sample_positions(attrs, fontsize) → (fig, ax)`
Renders a 384-well plate diagram with text labels at each well position showing the values of specified `Sample` attributes. `attrs` is a list of attribute paths (e.g. `['enzyme.id', 'ntp']`). Nested attributes are supported using dot notation. Use this to verify that experimental conditions have been assigned correctly before analysis.

#### `heatmap(numerator, denominator, vmin, vmax, title, cmap, hide_excluded) → (fig, ax)`
Generates a colour-coded 384-well plate heatmap.

- `numerator` and `denominator` can be any `Sample` attribute or the name of an MOI (which maps to its MAI intensity).
- If only `numerator` is given, the raw value is plotted.
- If both are given, the ratio `numerator / denominator` is plotted.
- `cmap` accepts any matplotlib colourmap name.

#### `stacked_plot(wells, xlim, title, filtered, figsize, label, label_first_only, overlay, hide_excluded, label_peaks) → (fig, axs)`
Plots multiple spectra stacked vertically (or overlaid) for comparison.

- `wells` — list of well IDs to include
- `xlim` — `(start, end)` tuple to zoom into a mass range
- `overlay=True` — draws all spectra on a single axis instead of stacking
- `label_peaks=True` — calls peaks and annotates m/z values on each spectrum
- MOI regions are highlighted in red; MOI names are labelled on the first spectrum unless `label_first_only=False`

#### `total_ion_plot() → (fig, ax)`
Bar chart of total ion count (`sum(sample.i)`) for every sample in the experiment. Useful for a quick quality check across the plate.

#### `sorted_signal_plot(xlim, label_noisy_samples) → (fig, ax)`
Plots the sorted intensity distribution for all samples together on one axis, with a vertical line at `noise_cutoff`. Setting `label_noisy_samples=True` annotates wells whose signal is still elevated at the noise cutoff.

---

### Data Collection and Export

#### `collect_mois() → list`
Gathers all `Oligo` objects assigned as MOIs across all samples, removes duplicates, and returns them sorted alphabetically by name.

#### `collect_misc_conditions() → list`
Collects all keys and values from `Sample.misc_conditions` across all samples, returning a deduplicated sorted list.

#### `collect_sample_data(headers, i_type) → dict`
Returns a nested dict of `{well: {header: value}}` for the given list of `headers`. Headers can be `Sample` attributes, `Enzyme` attributes prefixed with `'enz_'`, or MOI names. `i_type` controls which intensity array to use for MOI values (`'raw'`, `'bg_sub'`, or `'filtered'`).

#### `write_to_excel(filename, overwrite) → None`
Exports all experimental data to an `.xlsx` file with the following worksheets:

- **readme** — description of background/noise calculations and data definitions
- **mois** — name, sequence, composition, and masses for each molecule of interest
- **raw** — unprocessed ion intensities per sample per MOI
- **bg_sub** — background-subtracted intensities
- **filtered** — noise-filtered intensities (values below noise threshold zeroed out)

Set `overwrite=True` to allow overwriting an existing file (default is to abort if the file exists).

---