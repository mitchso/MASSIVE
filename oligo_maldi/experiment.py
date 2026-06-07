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


class Experiment:
    def __init__(self, data_folder=None, noise_cutoff=0.95):
        self.data_folder = data_folder
        self.noise_cutoff = noise_cutoff
        self.samples, self.calibrant_spots = self._data_folder_to_sample_dict()
        self.exclude = [] # used to remove certain wells from further analysis.

    def reinitiate_samples(self, samples):
        """
        Used to overwrite samples when loading them into an Experiment. Useful for when the
        noise_cutoff originally specified for the Sample is different than the global cutoff for the Experiment.
        """
        sample_dict = {}
        for s in samples:
            new_sample = Sample(file=s.file,
                                well=s.well,
                                mz=s.mz,
                                i=s.i,
                                noise_cutoff=self.noise_cutoff)
            sample_dict[s.well] = new_sample

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
            # [year]_[month]_[day]_[run#]_[chip]_[well]_[replicate].txt
            new_sample = Sample(chip=int(fields[4]),
                                file=file,
                                well=fields[5],
                                mz=[float(line.split(" ")[0]) for line in lines],
                                i=[int(line.split(" ")[1]) for line in lines],
                                noise_cutoff=self.noise_cutoff)

            if new_sample.chip == 0:
                sample_dict[new_sample.well] = new_sample

            elif new_sample.chip == 1:
                calibrant_dict[new_sample.well] = new_sample

        return sample_dict, calibrant_dict

    def collect_analytes(self) -> list:
        """
        Collects all analytes identified in experimental samples, removes duplicates,
        and returns them as an alphabetically sorted list.
        """
        analyte_list = []
        for s in self.samples.values():
            analyte_list.extend(s.analytes)

        unique_mois = list(set(analyte_list))
        sorted_mois = sorted(unique_mois, key=lambda x: x.name)

        return sorted_mois

    def scatter(self, x, y, regression=False) -> plt.axes:
        """
        Creates a scatter plot of any two variables across the entire experiment.

        """
        pass    # stub

    def show_sample_positions(self, attrs: list, fontsize=8) -> plt.axes:
        """
        Used to visualize sample definitions on a 384-well plate, so you know you've done it right.
        :return:
        """
        rows = list(ascii_uppercase[:16])  # ABCDEFGHIJKLMOP
        columns = [str(i) for i in range(1, 25)]  # 1-24

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(12, 18))
        im = ax.imshow(np.zeros((16, 24)), vmax=0, vmin=-1, cmap='gray')

        ntp_colours = {'G': 'blue',
                       'C': 'red',
                       'A': 'green',
                       'T': '#FFBF00'} # yellow

        for i, row in enumerate(rows):
            for j, col in enumerate(columns):
                well = row + str(col)

                if well in self.exclude:
                    ax.text(x=j, y=i, s=f"excl.", horizontalalignment='center', verticalalignment='center', fontsize=5)
                    continue

                try:
                    sample = self.samples[well]
                    strings = []
                    for attr in attrs:
                        attr_value = helper.get_nested_attr(sample, attr)

                        # if the attribute is a list, convert it to a string.
                        if isinstance(attr_value, list):
                            attr_value = ','.join(map(str, attr_value))

                        strings.append(str(attr_value))

                    final_string = "\n".join(strings)
                    ax.text(x=j, y=i, s=final_string, horizontalalignment='center', verticalalignment='center', fontsize=fontsize)

                except KeyError:
                    pass

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(range(len(columns)), labels=columns)
        ax.set_yticks(range(len(rows)), labels=rows)
        ax.tick_params(top=True, bottom=True,
                       labeltop=True, labelbottom=True, right=True, labelright=True)
        ax.tick_params(length=0, which='minor')

        ax.xaxis.set_minor_locator(tck.AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator(2))
        ax.grid(which="minor", color="black", linestyle='-', linewidth=0.25, snap=False)

        fig.tight_layout()
        return fig, ax


    def heatmap(self, numerator, denominator=None, vmin=0, vmax=None, title=None, cmap='viridis',
                hide_excluded=False) -> plt.axes:
        """
        Creates a heatmap of a variable across the 384 well sample plate.

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
        rows = list(ascii_uppercase[:16])   # ABCDEFGHIJKLMOP
        columns = [str(i) for i in range(1, 25)]    # 1-24

        data = []
        for row in rows:
            row_data = []
            for col in columns:
                key = row + col # ie 'C12'

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
                                    d = moi._mai_intensity(sample)
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
        fig, ax = plt.subplots()
        im = ax.imshow(data, vmin=vmin, vmax=vmax, cmap=cmap)
        cbar = ax.figure.colorbar(im, ax=ax, shrink=0.7)

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(range(len(columns)), labels=columns)
        ax.set_yticks(range(len(rows)), labels=rows)
        ax.tick_params(top=True, labeltop=True,
                       bottom=True, labelbottom=True,
                       right=True, labelright=True)
        ax.tick_params(length=0, which='minor')
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator(2))
        ax.grid(which="minor", color="black", linestyle='-', linewidth=0.25, snap=False)

        if title is None:
            if denominator:
                title = f"{numerator} / {denominator}"
            else:
                title = numerator

        ax.set_title(title)
        fig.tight_layout()
        return fig, ax


    def stacked_plot(
            self,
            wells: list,
            xlim=None,
            title=None,
            filtered=False,
            smoothed=False,
            figsize=None,
            label_analytes=True,
            highlight_analytes=True,
            label_first_only=True,
            overlay=False,
            hide_excluded=True,
            label_peaks=False,
            custom_colours=None,
            base_colour='#1f77b4',
            analyte_colour='#d1495b',
    ) -> plt.axes:
        """
        Takes a list of well IDs and returns a plot of each spectrum stacked
        vertically (or overlaid on a single axis if overlay=True).

        Each point is assigned exactly one colour by priority:
            base_colour (lowest) < MOI colour ('#d1495b') < custom_colour (highest)
        Contiguous same-colour runs are drawn as a single ax.plot() call, so
        no pixel is painted twice and there are no overlapping-line artifacts.

        Parameters
        ----------
        wells            : list of well IDs to plot
        xlim             : (start_mz, end_mz) window; uses full spectrum if None
        title            : figure suptitle
        filtered         : use noise-filtered intensities
        smoothed         : apply Savitzky-Golay smoothing
        figsize          : (width, height) in inches; auto-sized if None
        label            : annotate MOI names above their peaks
        label_first_only : only label MOIs on the first spectrum (reduces clutter)
        overlay          : draw all spectra on one axis instead of stacking
        hide_excluded    : skip wells listed in self.exclude
        label_peaks      : annotate detected peaks with their m/z value
        custom_colours   : {(start_mz, end_mz): colour_str, ...}
                           e.g. {(4830, 4843): '#0072B2', (4847, 4856): '#CC6677'}
        base_colour      : colour for regions with no special assignment
        """
        ### inital setup
        if hide_excluded:
            wells = [well for well in wells if well not in self.exclude]

        if figsize is None:
            figsize = (12, 3) if overlay or len(wells) == 1 else (12, len(wells))

        ### Make the figure
        plt.style.use('default')
        fig, axs = plt.subplots(
            nrows=1 if overlay else len(wells),
            ncols=1,
            figsize=figsize,
            sharex=True,
            sharey=True,
        )

        ### plot each sample
        for n, well in enumerate(wells):
            sample = self.samples[well]

            if overlay:
                ax = axs[1]
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


    def total_ion_plot(self) -> plt.axes:
        """
        Generates a bar plot of total ion intensity for each sample (a sum of all ion counts across the m/z range).
        """

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(16, 3))

        categories = []
        values = []
        for well, sample in self.samples.items():
            categories.append(well)
            values.append(sum(sample.i))

        plt.bar(categories, values)
        ax.set_xlim(-1, len(categories) + 1)
        plt.setp(ax.get_xticklabels(), fontsize=5, rotation='vertical')
        plt.title(f"Total ion count for each sample.")

        return fig, ax

    def sorted_signal_plot(self, xlim=(50,100), label_noisy_samples=False) -> plt.axes:
        """
        Visualizes the signal intensity of an entire sample, sorted from low to high and showing how the
        'noise_cutoff' parameter aligns with the data.
        """

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 3))

        # Collect all the data
        for well, sample in self.samples.items():
            x = np.linspace(0, 100, num=len(sample.i))
            y_max = max(sample.i)
            y = sorted([i * 100 / y_max for i in sample.i])
            plt.plot(x, y)

            if label_noisy_samples:
                y_at_cutoff = y[round(self.noise_cutoff * len(y))]
                if y_at_cutoff > 10:
                    plt.text(x=self.noise_cutoff * 100, y=y_at_cutoff + 1, s=f" {well}")

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
                        data[k][item] = a._mai_intensity(sample=v, i_type=i_type)

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

    def write_to_excel(self, filename: str, overwrite=False) -> None:
        """
        Writes data to an excel file.
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

