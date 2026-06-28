## Oligo

**File:** `oligo.py`

Represents a DNA or RNA oligonucleotide. On instantiation, the chemical composition, isotopic distribution, monoisotopic mass, and average mass are all computed automatically.

### Constructor

```python
Oligo(seq, name=None, type='DNA', five_prime_end='OH', ps_bonds=0,
      error=0, charge=1, custom_mods=None)
```

| Parameter | Description |
|---|---|
| `seq` | Nucleotide sequence string (e.g. `'ATCG'`) |
| `name` | Human-readable label |
| `type` | `'DNA'` or `'RNA'` |
| `five_prime_end` | 5′ end modification: `'OH'`, `'P'`, `'PP'`, or `'PPP'` |
| `ps_bonds` | Number of phosphorothioate bonds |
| `error` | Manual m/z offset to apply to the isotopic distribution |
| `charge` | Charge state for m/z calculation |
| `custom_mods` | Dictionary of elemental composition adjustments |

**Key attributes after construction:**

- `composition` — elemental formula as a dict (e.g. `{'C': 40, 'H': 52, ...}`)
- `monoisotopic_mass` — mass of the most abundant isotope (3 decimal places)
- `average_mass` — intensity-weighted average mass (3 decimal places)
- `isotopic_distribution` — list of isotope peaks from `brainpy`
- `iso_dist_range` — mass range that corresponds to this Oligo. This is automatically calculated using `calc_iso_dist_range`, but the pre-computed range can be overwritten by assigning a new range to this attribute.

---

### Methods

#### `chem_composition(seq, five_prime_end, ps_bonds, custom_mods) → dict`
Calculates the elemental formula for the oligonucleotide. Called automatically on construction. Accounts for base composition, backbone phosphodiester bonds, phosphorothioate substitutions, 5′ end chemistry, and any custom modifications.

#### `calc_iso_dist(charge, error) → list`
Computes the isotopic distribution using `brainpy.isotopic_variants`. Returns a list of peaks, each with `.mz` and `.intensity`. Called automatically on construction.

#### `calc_monoisotopic_mass() → float`
Returns the m/z of the first (monoisotopic) peak, rounded to 3 decimal places.

#### `calc_avg_mass() → float`
Returns the intensity-weighted average m/z across all isotope peaks, rounded to 3 decimal places.

#### `iso_dist_plot(y_max, standalone, annotate, label, colour) → (fig, ax)`
Plots the theoretical isotopic distribution as a stem plot.

- If `standalone=True` (default), returns a new `(fig, ax)` tuple.
- If `standalone=False`, overlays onto an existing figure.
- `y_max` scales the distribution to match a specified intensity.
- `annotate=True` labels each peak with its m/z value.

#### `iso_dist_range(cumulative_threshold, left_pad, right_pad) → tuple`
Returns a `(start, end)` m/z range covering the cumulative isotope signal up to `cumulative_threshold` (default 95%), with optional padding on each side. Used internally to slice spectra around a molecule of interest.

#### `mai_intensity(sample, i_type) → float`
Returns the Most Abundant Isotope (MAI) intensity of this `Oligo` within a given `Sample`.

- `i_type` controls which intensity array to use: `'raw'`, `'bg_sub'`, or `'filtered'` (default).
- Returns `0` if no signal is found in the expected m/z range.

#### `composition_str() → str`
Returns the chemical composition as a human-readable string (e.g. `'C40 H52 N15 O24 P4 S0'`).

---