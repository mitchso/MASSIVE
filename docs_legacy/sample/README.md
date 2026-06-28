## Sample

**File:** `sample.py`

Holds a single mass spectrum and performs background subtraction, noise filtering, and visualization. Typically created automatically by `Experiment`.

### Constructor

```python
Sample(file, well, mz, i, noise_cutoff, chip=0, mz_offset=0)
```

| Parameter | Description |
|---|---|
| `file` | Source filename |
| `well` | Well ID (e.g. `'E1'`) |
| `mz` | List of m/z values |
| `i` | List of intensity values |
| `noise_cutoff` | Percentile threshold for noise calculation (e.g. `0.95` = 95th percentile) |
| `chip` | Chip number (`0` = sample, `1` = calibrant) |
| `mz_offset` | Manual m/z correction offset |

**Key attributes after construction:**

- `i_bg_subtracted` — intensity after background subtraction
- `i_filtered` — intensity after noise filtering (values below noise set to 0)
- `background` — minimum intensity across the spectrum (used for background subtraction)
- `noise` — intensity threshold at `noise_cutoff` percentile
- `mois` — list of `Oligo` objects set as molecules of interest (assigned manually)
- `enzyme` — `Enzyme` object assigned to this sample
- `ntp` — nucleotide triphosphate condition label
- `misc_conditions` — dict for any additional experimental metadata

---

### Methods

#### `recalc_mz()`
Recomputes `self.mz` from the original raw values plus `self.mz_offset`. Call after changing `mz_offset` to correct a calibration shift.

#### `moi_intensities() → list[tuple]`
Returns a list of `(molecule_name, intensity)` tuples for each molecule of interest assigned to this sample. Intensity is determined by `Oligo.mai_intensity()`.

#### `calculate_noise(cutoff) → int`
Computes the noise level as the maximum intensity among the lowest `cutoff` fraction of background-subtracted values.

#### `filter_i()`
Populates `self.i_filtered` by zeroing out any background-subtracted intensity values that fall below the noise threshold.

#### `call_peaks(i_list) → list[int]`
Applies `scipy.signal.find_peaks` to find peak indices in a given intensity list. Returns a list of integer indices. Useful for exploratory inspection of a spectrum.

#### `slice_spectrum(start, end, mz, i) → tuple`
Extracts the portion of a spectrum between `start` and `end` m/z values. Returns `(mz_slice, i_slice)` as lists.

#### `plot(xlim, relative, label_peaks, theoretical_dist, label_mois, filtered, bg_subtracted, title) → (fig, ax)`
General-purpose spectrum plot.

- `xlim` — tuple `(start, end)` to zoom into a region; shows full spectrum if `None`
- `relative` — if `True`, normalizes intensities to 100%
- `label_peaks` — if `True`, marks and annotates called peaks
- `theoretical_dist` — if `True`, overlays theoretical isotope distributions for all MOIs
- `label_mois` — if `True`, annotates MOI positions by name
- `filtered` — if `True`, uses noise-filtered intensities

#### `plot_moi(moi, filtered) → (fig, ax)`
A convenience wrapper around `plot()` that zooms into the isotope distribution range of a specified `Oligo` MOI and adds a descriptive title with monoisotopic mass and charge state.

#### `sorted_signal_plot() → (fig, ax)`
Plots intensity values sorted from lowest to highest across the spectrum, with a vertical line at `noise_cutoff`. Useful for evaluating signal quality and choosing an appropriate noise threshold.

#### `plot_distance_between_points() → (fig, ax)`
Plots the m/z gap between consecutive data points. Helpful for identifying gaps or irregularities in instrument data collection.

#### `misc_conditions_to_attributes()`
Promotes all key-value pairs in `self.misc_conditions` to top-level attributes of the `Sample` object. Called automatically by `Experiment.write_to_excel()`.

---