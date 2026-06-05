# Nucleic Acid MALDI Analysis — Wiki

A reference guide for the classes and functions in this codebase, which provides tools for visualizing and quantifying nucleic acid MALDI experiments.

---

## Table of Contents

1. [Overview](#overview)
2. [Data Format](#data-format)
3. [Oligo](#oligo)
4. [Enzyme](#enzyme)
5. [Sample](#sample)
6. [Experiment](#experiment)
7. [DataProcessor (quant.py)](#dataprocessor)
8. [Helper Functions](#helper-functions)

---

## Overview

Analysis is built around four primary objects that work together in a typical pipeline:

```
Oligo / Enzyme  →  Sample  →  Experiment  →  DataProcessor
(define species)   (load spectra) (manage plate)  (quantify & plot)
```

- **`Oligo`** — represents a nucleic acid molecule and its theoretical mass/isotope distribution
- **`Enzyme`** — represents an enzyme variant with mutation metadata and activity scoring
- **`Sample`** — holds a single mass spectrum (one well) and performs signal processing
- **`Experiment`** — manages a full plate of `Sample` objects and handles data export
- **`DataProcessor`** — processes Excel output from `Experiment` for quantitative analysis and plotting

---

## Data Format

Raw data files should be plain text (`.txt`) with two space-separated columns (m/z and intensity) and no header row. Filenames must follow this convention:

```
[year]_[month]_[day]_[run#]_[chip#]_[well]_[replicate].txt
# e.g. 2025_12_06_0001_0_E1_1.txt
```

`chip#` is `0` for regular sample spots and `1` for calibrant spots. `well` is the plate position (e.g. `A1`, `G12`).

---









## DataProcessor

**File:** `quant.py`

Processes the Excel output from `Experiment.write_to_excel()` for quantitative analysis. Accepts one or more pandas DataFrames and produces averaged, percentage-normalised, and activity-scored datasets.

### Constructor

```python
DataProcessor(dataframes, attributes, exclude_mois=None, collapse_mois=None,
              activity_scoring_function=None, activity_scoring_df_variables=None,
              activity_scoring_other_variables=None)
```

| Parameter | Description |
|---|---|
| `dataframes` | Dict of `{name: DataFrame}` — one entry per plate or dataset |
| `attributes` | List of column names that define a unique experimental condition (e.g. `['ntp', 'enz_name', 'enz_id']`). It is critical to include *all* grouping variables, or unrelated samples may be averaged together. |
| `exclude_mois` | List of MOI column names to drop before processing |
| `collapse_mois` | Dict mapping a new MOI name to a list of old MOI names that should be merged under it (e.g. `{'+N, PPP': ['+A, PPP', '+G, PPP']}`) |
| `activity_scoring_function` | A callable (e.g. `Enzyme.default_activity_score`) for scoring each row |
| `activity_scoring_df_variables` | Maps scoring function parameter names to DataFrame column names |
| `activity_scoring_other_variables` | Additional fixed keyword arguments passed to the scoring function |

**Key attributes after construction:**

| Attribute | Description |
|---|---|
| `data_ungrouped` | All input dataframes concatenated with a `source` column |
| `data_mois_collapsed` | After applying MOI renaming and exclusions |
| `data_replicates_averaged` | Technical replicates averaged by `attributes` |
| `data_as_percentage` | MOI intensities normalised to % of total per row |
| `dfs` | Dict containing all of the above DataFrames, plus `'activity_sums'` and `'evolvepro'` if scoring was enabled |

---

### Processing Methods

#### `collapse_input_dataframes() → DataFrame`
Concatenates all input DataFrames into one, adding a `source` column to track origin.

#### `collapse_mois(df, mapping_dict) → DataFrame`
Renames and merges MOI columns according to `mapping_dict`, then drops the original columns. Also drops any columns in `excluded_mois`.

#### `collapse_replicates(df, attributes) → DataFrame`
Groups by `attributes` and averages all numeric columns, collapsing technical replicates.

#### `convert_to_percentage(df) → DataFrame`
Converts absolute MOI intensities to fractional values (each row sums to 1.0) by dividing by the row total across all MOI columns.

#### `add_activity_scores() → DataFrame`
Applies `activity_scoring_function` to each row using the specified variable mappings, adding an `activity_score` column to `data_as_percentage`.

#### `sum_activities() → DataFrame`
Groups the scored data by `enz_id` and `enz_mutation`, summing activity scores across all nucleotide conditions. Returns a summary DataFrame.

---

### Visualization Methods

#### `stacked_bar(global_var, x_category, y_species, y_colours, x_labels, x_labels_overwrite, y_legend_overwrite, sort_by, sort_order, hide_points, figsize, title, show_legend, ax) → ax`
The primary visualization method. Generates a stacked bar chart of MOI fractions, with individual replicate data points overlaid as a swarm plot.

- `global_var` — `(column, value)` tuple to filter the dataset (e.g. `('ntp', 'G')`)
- `x_category` — column whose unique values become the x-axis categories
- `y_species` — list of MOI column names to stack
- `y_colours` — list of hex colour strings, one per MOI
- `sort_by` — MOI name(s) to sort x-axis categories by
- `sort_order` — `'ascending'` or `'descending'`
- `ax` — optionally pass an existing matplotlib `Axes` object to embed the plot in a larger figure

#### `score_bar(x_category, x_labels, y_species, y_colours, title, ax, show_legend) → ax`
Generates a stacked bar chart of activity scores broken down by nucleotide condition. Intended for comparing enzyme variants.

#### `sort_sample_labels(df_name, label_category, sort_by, global_var, ascending) → list`
Returns a sorted list of labels from a specified DataFrame column. Useful for controlling x-axis order when calling `stacked_bar`.

#### `write_to_excel(filename, overwrite) → None`
Writes all DataFrames in `self.dfs` to separate worksheets in a single `.xlsx` file.

---

## Helper Functions

**File:** `helper.py`

Utility functions for mapping plate well positions to enzyme and nucleotide assignments.

#### `key_to_enz_code(key, enz_rows, enz_cols) → int`
Converts a well ID (e.g. `'C5'`) to an enzyme code by cross-referencing which rows and columns that enzyme occupies, as defined in `enz_rows` and `enz_cols` dictionaries.

```python
# Example dictionaries:
enz_rows = {1: ['A', 'B', 'C'], 2: ['D', 'E', 'F']}
enz_cols = {1: [1, 2, 3], 2: [4, 5, 6]}

key_to_enz_code('B5', enz_rows, enz_cols)  # → ambiguous; row B → enz 1, col 5 → enz 2
```

#### `key_to_dntp(key, dntp_cols, exclude) → str`
Maps a well ID to a nucleotide triphosphate label by checking which column group it belongs to, as defined in `dntp_cols`. Wells in the `exclude` list are skipped and return `None`.

```python
dntp_cols = {'G': [1, 2], 'C': [3, 4], 'A': [5, 6], 'T': [7, 8]}
key_to_dntp('A3', dntp_cols)  # → 'C'
```

---

> **Note on deprecated functions:** Several methods across `helper.py`, `sample.py`, and `experiment.py` are marked `@deprecated`. These were used with an earlier XML-based data format from Bruker DataAnalysis software and are not compatible with the current plain-text file format. They are retained for reference but should not be used in new analyses.
