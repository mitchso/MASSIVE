import os
import random
import xlsxwriter
from . import helper
from .sample import Sample
from scipy.signal import find_peaks
from string import ascii_uppercase
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.axes import Axes
from matplotlib.figure import Figure


# TODO: implement peak matching algorithms
    # exploratory - sample to sample comparisons
    # directed - sample to analyte comparisons

class Collection:
    """
    General collection of [`Samples`][MASSIVE.sample.Sample]. All `Samples` in a `Collection` can be analyzed simultaneously, and `Collection` contains numerous methods for visualizing and quantifying the data.

    If all [`Samples`][MASSIVE.sample.Sample] were acquired in a single MALDI-ToF run, consider using the subclass [`Plate`][MASSIVE.collections.Plate]. This subclass contains methods for visualizing the data and metadata in heatmap format across the plate.

    """
    def __init__(self, data_folder:str|None=None, noise_cutoff:float=0.95):
        """

        Args:
            data_folder: Data folder path (see [Input data formatting assistance](../guides/data_formatting.md)). If `None` is given, instantiates an empty `Collection`.
            noise_cutoff: Sets a global cutoff for every `Sample` in the `Collection`. See [`Sample`][MASSIVE.sample.Sample].noise_cutoff for explanation of the effect.

        Attributes:
            data_folder (str|None): Data folder path.
            noise_cutoff (float): Global noise cutoff.
            samples (dict): Dictionary of `Sample` objects where `Sample.chip = 0`, keyed by sample name.
            calibrant_spots: Dictionary of `Sample` objects where `Sample.chip = 1`, keyed by sample name. See [`Sample`][MASSIVE.sample.Sample].chip for explanation of the difference between `chip = 0` and `chip = 1`.
        """
        self.data_folder = data_folder
        self.noise_cutoff = noise_cutoff
        self.samples, self.calibrant_spots = self._data_folder_to_sample_dict()
        self.exclude = []  # used to remove certain Samples from further analysis.


    def reinitiate_samples(self, samples):
        """
        Used to overwrite samples when loading them into an Experiment. Useful for when the
        noise_cutoff originally specified for the Sample is different than the global cutoff for the Experiment.
        """
        sample_dict = {}
        for s in samples:
            new_sample = Sample(file=s.file,
                                id=s.id,
                                mz=s.mz,
                                i=s.i,
                                noise_cutoff=self.noise_cutoff)
            sample_dict[s.id] = new_sample

        return sample_dict


    def _data_folder_to_sample_dict(self):
        """
        Takes a folder of txt files where each file corresponds to 1 spectra, encoded as m/z - intensity pairs.
        Returns a sample dictionary where each entry is 1 spectra, named according to filenames.

        example data file:
        /Users/mitchsyberg-olsen/Library/CloudStorage/Box-Box/lab/data/2025/12/2025-12-05_new_substrate_test/data/2025_12_06_0001_0_E1_1.txt
        """

        for root, dirs, files in os.walk(self.data_folder):
            root = root
            files = files
            break   # stop at first level of folder

        sample_dict = {}
        calibrant_dict = {}
        for file in files:
            if file == ".DS_Store":
                continue
            assert file.endswith(".txt"), f"Wrong file suffix: {file}"    # checks to make sure you are working with the correct files
            fields = file.split("_")
            with open(root + "/" + file, "r") as f:
                lines = f.readlines()

            # assumes name follows this structure:
            #   0      1      2      3      4      5        6
            # 2025_12_06_0001_0_E1_1.txt
            # [year]_[month]_[day]_[run#]_[chip]_[id]_[replicate].txt
            new_sample = Sample(chip=int(fields[4]),
                                file=file,
                                id=fields[5],
                                mz=[float(line.split(" ")[0]) for line in lines],
                                i=[int(line.split(" ")[1]) for line in lines],
                                noise_cutoff=self.noise_cutoff)

            if new_sample.chip == 0:
                sample_dict[new_sample.id] = new_sample

            elif new_sample.chip == 1:
                calibrant_dict[new_sample.id] = new_sample

        return sample_dict, calibrant_dict


    def collect_analytes(self) -> list:
        """
        Collects all analytes identified in experimental samples, removes duplicates,
        and returns them as an alphabetically sorted list.
        """
        analyte_list = []
        for s in self.samples.values():
            analyte_list.extend(s.analytes)

        unique_analytes = list(set(analyte_list))
        sorted_analytes = sorted(unique_analytes, key=lambda x: x.name)

        return sorted_analytes


    def scatter(self, x, y, regression=False) -> plt.axes:
        """
        Creates a scatter plot of any two variables across the entire collection.

        Warning:
            This section is under active development and may change without notice.

        """
        pass    # stub


    def plot(
            self,
            ids: list,
            xlim: tuple=None,
            title: str|None=None,
            filtered: bool=False,
            smoothed: bool=False,
            figsize:tuple|None=None,
            label_analytes:bool=True,
            highlight_analytes:bool=True,
            label_first_only:bool=True,
            overlay:bool=False,
            hide_excluded:bool=True,
            label_peaks:bool=False,
            custom_colours:dict|bool=None,
            base_colour:str='#1f77b4',
            analyte_colour:str='#d1495b',
    ) -> tuple[Figure, Axes]:
        """
        Takes a list of [`Sample.id`][MASSIVE.sample.Sample] and returns a plot of each spectrum stacked
        vertically (or overlaid on a single axis if `overlay=True`).


        Args:
            ids: list of IDs to plot, in order (top to bottom)
            figsize: (width, height) in inches; auto-sized if None
            title: figure suptitle
            overlay: draw all spectra on one axis instead of stacking
            label_first_only: only label MOIs on the first spectrum (reduces clutter)
            hide_excluded: skip IDs listed in `self.exclude`
            xlim: Passed to [`Sample.plot()`][MASSIVE.sample.Sample.plot]
            filtered: Passed to [`Sample.plot()`][MASSIVE.sample.Sample.plot]
            smoothed: Passed to [`Sample.plot()`][MASSIVE.sample.Sample.plot]
            label_analytes: Passed to [`Sample.plot()`][MASSIVE.sample.Sample.plot]
            highlight_analytes: Passed to [`Sample.plot()`][MASSIVE.sample.Sample.plot]
            label_peaks: Passed to [`Sample.plot()`][MASSIVE.sample.Sample.plot]
            custom_colours: Passed to [`Sample.plot()`][MASSIVE.sample.Sample.plot]
            base_colour: Passed to [`Sample.plot()`][MASSIVE.sample.Sample.plot]
            analyte_colour: Passed to [`Sample.plot()`][MASSIVE.sample.Sample.plot]

        """
        ### inital setup
        if hide_excluded:
            ids = [id for id in ids if id not in self.exclude]

        if figsize is None:
            figsize = (12, 3) if overlay or len(ids) == 1 else (12, len(ids))

        ### Make the figure
        plt.style.use('default')
        fig, axs = plt.subplots(
            nrows=1 if overlay else len(ids),
            ncols=1,
            figsize=figsize,
            sharex=True,
            sharey=True,
        )

        ### plot each sample
        for n, id in enumerate(ids):
            sample = self.samples[id]

            if overlay:
                ax = axs
                custom_colours = None
                base_colour = None
                analyte_colour = None
            elif len(ids) == 1:
                ax = axs
            else:
                ax = axs[n]

            if label_first_only and n != 0:
                label_analytes = False

            sample._base_plotter(ax=ax,
                                 xlim=xlim,
                                 filtered=filtered,
                                 smoothed=smoothed,
                                 normalized=True,
                                 overlay=overlay,
                                 highlight_analytes=highlight_analytes,
                                 label_analytes=label_analytes,
                                 label_peaks=label_peaks,
                                 custom_colours=custom_colours,
                                 base_colour=base_colour,
                                 analyte_colour=analyte_colour,
                                 linewidth=1.5)

            ### Formatting
            for spine in ['top', 'right', 'left']:
                ax.spines[spine].set_visible(False)
            ax.spines['bottom'].set_position(('data', -5))

            if not overlay:
                x_min, x_max = ax.get_xlim()
                ax.text(x=x_min, y=0, s=f"{sample.name} ", horizontalalignment='right')

        ### Final formatting
        ax.set_ylim(-10, 150)   # only have to set one time because share_y is True
        plt.xlabel('m/z')
        plt.yticks([])
        if title:
            plt.suptitle(title)

        return fig, axs


    def total_ion_plot(self) -> tuple[Figure, Axes]:
        """
        Generates a bar plot of total ion intensity for each sample (a sum of all ion counts across the m/z range).
        """

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(16, 3))

        categories = []
        values = []
        for id, sample in self.samples.items():
            categories.append(id)
            values.append(sum(sample.i))

        plt.bar(categories, values)
        ax.set_xlim(-1, len(categories) + 1)
        plt.setp(ax.get_xticklabels(), fontsize=5, rotation='vertical')
        plt.title(f"Total ion count for each sample.")

        return fig, ax

    def sorted_signal_plot(self, xlim=(50,100), label_noisy_samples=False) -> tuple[Figure, Axes]:
        """
        Visualizes the signal intensity of an entire sample, sorted from low to high and showing how the
        'noise_cutoff' parameter aligns with the data.
        """

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 3))

        # Collect all the data
        for id, sample in self.samples.items():
            x = np.linspace(0, 100, num=len(sample.i))
            y_max = max(sample.i)
            y = sorted([i * 100 / y_max for i in sample.i])
            plt.plot(x, y)

            if label_noisy_samples:
                y_at_cutoff = y[round(self.noise_cutoff * len(y))]
                if y_at_cutoff > 10:
                    plt.text(x=self.noise_cutoff * 100, y=y_at_cutoff + 1, s=f" {id}")

        plt.xlim(xlim)
        plt.yticks(np.linspace(0, 100, num=11))

        plt.axvline(x=self.noise_cutoff * 100)
        # plt.text(x=self.noise_cutoff * 100, y=0, s=" Noise cutoff")
        if label_noisy_samples:
            plt.text(x=self.noise_cutoff * 100, y=50, s=" Samples w/\n S/N<10")
        plt.ylim(ymin=0)
        plt.title(f"Per sample signal intensity distribution")
        plt.xlabel(f'Percent of data')
        plt.ylabel(f'Relative signal')

        return fig, ax

    def _collect_sample_data(self, headers: list, i_type='filtered') -> dict:
        """
        Collects general sample information and the most abundant ion intensity for each sample.
        """

        data = {}
        for k, v in self.samples.items():
            if k in self.exclude:
                continue

            data[k] = {}
            for item in headers:
                # first, try to see if it's a normal attribute
                try:
                    attr = getattr(v, item)
                    if isinstance(attr, (int, float)):
                        data[k][item] = attr
                    else:
                        data[k][item] = str(attr)
                except AttributeError:
                    pass

                # if that doesn't work, it could be an analyte
                for a in v.analytes:
                    if a.name == item:
                        data[k][item] = a._peak_intensity(sample=v, i_type=i_type)

        return data

    def _write_excel_readme(self, workbook):
        """Adds a README worksheet to the data output workbook. Use with write_to_excel()."""

        # write README
        readme = f"""
                analytes: Contains all analyte and their associated information.

                raw: Unprocessed data. Each entry corresponds to the ion intensity of the most abundant isotope of a given analyte.

                bg_sub: Raw data with background subtraction.

                filtered: Background subtracted data with noise removed. Any value that appears here is in the 90th percentile 
                of signal intensity for a given sample.

                Background calculation:
                    Background = lowest value across the entire spectrum for a given sample.

                Noise calculation: 
                    Noise = For a given spectrum, anything less than the 90th percentile of all data, in terms of absolute ion intensity.
                    90th percentile is chosen by default, but can be changed when an Experiment object is instantiated using the 
                    'noise_cutoff' parameter.
                """

        worksheet = workbook.add_worksheet('readme')
        worksheet.insert_textbox('B2', readme, {'width': 1000, 'height': 1000})

    def _write_excel_analytes(self, workbook):
        """Adds a worksheet to the data output workbook, containing information about.
         All the analytes in the experiment.
         Use with write_to_excel()."""

        worksheet = workbook.add_worksheet('analytes')

        # write analyte data
        headers = ['name', 'sequence', 'composition', 'monoisotopic mass', 'average mass']
        for i, header in enumerate(headers):  # write headers
            worksheet.write(0, i, header)

        analytes = self.collect_analytes()
        for i, a in enumerate(analytes):  # i + 1 for row
            data = [a.name, a.seq, a.composition_str(), a.monoisotopic_mass, a.average_mass]
            for j, ele in enumerate(data):  # j for column
                worksheet.write(i + 1, j, ele)

    def _collect_attributes_for_excel(self):
        """Collects all attributes from samples that should be written to excel."""

        # collect regular attributes
        attribute_keys = []
        for s in self.samples.values():
            for k, v in s._collect_attributes().items():
                if k not in attribute_keys:
                    attribute_keys.append(k)

        # collect moi names
        analyte_attributes = [a.name for a in self.collect_analytes()]

        return attribute_keys + analyte_attributes

    def _write_excel_data(self, workbook):
        """Adds worksheets to the data output workbook, containing all experimental data."""
        # write data sheets
        headers = self._collect_attributes_for_excel()

        i_types = ['raw', 'bg_sub', 'filtered']
        for i_type in i_types:
            worksheet = workbook.add_worksheet(i_type)  # create a worksheet
            for i, header in enumerate(headers):  # write headers
                worksheet.write(0, i, header)

            # write data
            i = 1
            data = self._collect_sample_data(headers=headers, i_type=i_type)  # collect data according to i_type param

            # lambda function sorts keys by letter (row) and then number (column)
            sorted_keys = sorted(data.keys(), key=lambda x: (x[0], int(x[1:])))

            for key in sorted_keys:
                value = data[key]
                for j, header in enumerate(headers):  # iterate through headers
                    try:
                        worksheet.write(i, j, value[header])  # if data matches header, write
                    except KeyError:
                        worksheet.write(i, j, "n/a")  # if no data exists for this header, write n/a
                i += 1

    def write_to_excel(self, filename: str, overwrite:bool=False) -> None:
        """
        Compiles data about the entire collection and writes to an excel file. Call this method after assigning [`Analytes`][MASSIVE.analytes.Analyte] to the [`Samples`][MASSIVE.sample.Sample] in the [`Collection`][MASSIVE.collections.Collection], so that the peak intensity values for each [`Analyte`][MASSIVE.analytes.Analyte] will be calculated and included.

        Args:
            filename: Absolute path to write the file.
            overwrite: If True, overwrite the file if it already exists.

        """

        # create file
        outfile = filename if filename.endswith(".xlsx") else filename + ".xlsx"

        if not overwrite:
            if os.path.exists(outfile):
                print(f"\'write_to_excel()\' did not execute because the file already exists.\n"
                      f"To proceed anyway, use parameter \'overwrite=True\'.")
                return

        workbook = xlsxwriter.Workbook(outfile)
        self._write_excel_readme(workbook)
        self._write_excel_analytes(workbook)
        self._write_excel_data(workbook)
        workbook.close()


class Plate(Collection):
    """
    Extends the [`Collection`][MASSIVE.collections.Collection] class to include methods for visualizing and analyzing on a 96-, 384-, or 1536-target plate.
    """
    def __init__(self, data_folder:str=None, num_targets:int=384, noise_cutoff:float=0.95):
        """
        Args:
            data_folder: Data folder path (see [Input data formatting assistance](../guides/data_formatting.md)). If `None` is given, instantiates an empty `Collection`.
            noise_cutoff: Sets a global cutoff for every `Sample` in the `Collection`. See [`Sample`][MASSIVE.sample.Sample].noise_cutoff for explanation of the effect.
            num_targets: Number of targets on the plate. Must be 96, 384, or 1536.

        """
        super().__init__(data_folder=data_folder, noise_cutoff=noise_cutoff)
        self.num_targets = num_targets
        self.rows, self.columns = self._assign_rows_columns()
        self._verify_plate_positions()



    def _verify_plate_positions(self):
        accepted_positions = self._plate_positions()

        bad_positions = []
        for id in self.samples.keys():
            if id not in accepted_positions:
                bad_positions.append(id)

        if bad_positions:
            raise ValueError(f"The following sample IDs do not correspond to a location on a {self.num_targets}-target plate: {str(bad_positions)}.\n"
                             f"Please correct the positions, or change Plate.num_targets.\n"
                             f"If your samples are not organized by plate map, use the parent class Collection instead of Plate.\n")

    def _assign_rows_columns(self) -> tuple:
        configs = {
            96:   (8, 12),  # row, col
            384:  (16, 24),
            1536: (32, 48),
        }

        if self.num_targets not in configs:
            raise ValueError(f"Plate.num_targets must be 96, 384, or 1536, got {self.num_targets}.")

        num_rows, num_cols = configs[self.num_targets]
        letters = list(ascii_uppercase[:num_rows]) + ['AA', 'AB', 'AC', 'AD', 'AE', 'AF']
        rows = letters[:num_rows]
        cols = [n+1 for n in range(0, num_cols)]

        return rows, cols

    def _plate_positions(self) -> list:
        positions = []
        for row in self.rows:
            for col in self.columns:
                positions.append(row + str(col))
        return positions

    def show_sample_positions(self, attrs: list, fontsize=8) -> plt.axes:
        """
        Used to visualize sample definitions on the target plate, so you know you've done it right.

        Warning:
            This section is under active development and may change without notice.

        """
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(12, 18))
        im = ax.imshow(np.zeros((len(self.rows), len(self.columns))), vmax=0, vmin=-1, cmap='gray')

        ntp_colours = {'G': 'blue',
                       'C': 'red',
                       'A': 'green',
                       'T': '#FFBF00'} # yellow

        for i, row in enumerate(self.rows):
            for j, col in enumerate(self.columns):
                id = row + str(col)

                if id in self.exclude:
                    ax.text(x=j, y=i, s=f"excl.", horizontalalignment='center', verticalalignment='center', fontsize=5)
                    continue

                try:
                    sample = self.samples[id]
                    strings = []
                    for attr in attrs:
                        attr_value = helper._get_nested_attr(sample, attr)

                        # if the attribute is a list, convert it to a string.
                        if isinstance(attr_value, list):
                            attr_value = ','.join(map(str, attr_value))

                        strings.append(str(attr_value))

                    final_string = "\n".join(strings)
                    ax.text(x=j, y=i, s=final_string, horizontalalignment='center', verticalalignment='center', fontsize=fontsize)

                except KeyError:
                    pass

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(range(len(self.columns)), labels=self.columns)
        ax.set_yticks(range(len(self.rows)), labels=self.rows)
        ax.tick_params(top=True, bottom=True,
                       labeltop=True, labelbottom=True, right=True, labelright=True)
        ax.tick_params(length=0, which='minor')

        ax.xaxis.set_minor_locator(tck.AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator(2))
        ax.grid(which="minor", color="black", linestyle='-', linewidth=0.25, snap=False)

        fig.tight_layout()
        return fig, ax


    def heatmap(self, numerator, denominator=None, vmin=0, vmax=None, title=None, cmap='viridis',
                hide_excluded=False) -> tuple[Figure, Axes]:
        """
        Warning:
            This section is under active development and may change without notice.

        Creates a heatmap of a variable across the 384 target sample plate.

        numerator and denominator
            Possible variables:
            Sample total_ion, noise, background, or any MOI.
            If an MOI is indicated, the value for moi.mai_intensity() for that MOI will be plotted.
            If only numerator is specified, the value of that variable will be plotted.
            If numerator and denominator are specified, the ratio of the two will be plotted.

        vmin: lower bound of data range covered by colour map
        vmax: upper bound of data range covered by colour map
        cmap: colormap to use (see docs online for pyplot.imshow()
        title: title of the plot
        """
        n_rows = len(self.rows)
        n_cols = len(self.columns)

        data = []
        for row in self.rows:
            row_data = []
            for col in self.columns:
                key = row + str(col) # ie 'C12'

                if hide_excluded is True and key in self.exclude:
                    row_data.append(np.nan)
                    continue

                try:
                    sample = self.samples[key]  # If key exists, collect data
                    for variable in [numerator, denominator]:
                        d = np.nan  # by default, return NaN
                        try:    # check if variable is an attribute of the sample
                            d = getattr(sample, variable)
                        except AttributeError:  # if it isn't, check MOIs
                            for moi in sample.analytes:
                                if moi.name == variable:
                                    d = moi._peak_intensity(sample)
                                    break
                        except TypeError:   # occurs if denominator=None
                            pass

                        if variable == numerator:
                            num = d
                        elif variable == denominator:
                            denom = d

                    if denominator:  # return the ratio of num and denom
                        if numerator == denominator:
                            d = 1
                        else:
                            try:
                                d = num / denom
                            except ZeroDivisionError:
                                d = np.nan
                    else:   # return numerator only
                        d = num

                except KeyError:    # If key does not exist, leave blank
                    d = np.nan

                row_data.append(d)
            data.append(row_data)

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(n_cols * 0.4, n_rows * 0.4))
        im = ax.imshow(data, vmin=vmin, vmax=vmax, cmap=cmap)
        cbar = ax.figure.colorbar(im, ax=ax, shrink=0.7)

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(range(n_cols), labels=self.columns)
        ax.set_yticks(range(n_rows), labels=self.rows)
        ax.tick_params(top=True, labeltop=True,
                       bottom=True, labelbottom=True,
                       right=True, labelright=True)
        ax.tick_params(length=0, which='minor')

        # Manually draw grid lines
        for x in np.arange(-0.5, n_cols, 1):
            ax.axvline(x, color='black', linewidth=0.25)
        for y in np.arange(-0.5, n_rows, 1):
            ax.axhline(y, color='black', linewidth=0.25)

        if title is None:
            if denominator:
                title = f"{numerator} / {denominator}"
            else:
                title = numerator

        ax.set_title(title)
        fig.tight_layout()
        return fig, ax

