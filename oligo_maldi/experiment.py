import os
import random
import xlsxwriter
import pyopenms as oms
from .helper import *
from warnings import deprecated
from scipy.signal import find_peaks
from string import ascii_uppercase
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck


class Experiment:
    def __init__(self, data_folder=None, noise_cutoff=0.95):
        self.data_folder = data_folder
        self.noise_cutoff = noise_cutoff
        self.samples, self.calibrant_spots = self.data_folder_to_sample_dict()
        self.exclude = [] # used to remove certain wells from further analysis.

    def reinitiate_samples(self, samples):
        """
        Used to overwrite samples when loading them into an Experiment. Useful for when the
        noise_cutoff originally specified for the Sample is different than the global cutoff for the Experiment.
        """
        sample_dict = {}
        for s in samples:
            new_sample = Sample(well=s.well,
                                mz=s.mz,
                                i=s.i,
                                noise_cutoff=self.noise_cutoff)
            sample_dict[s.well] = new_sample

        return sample_dict

    def data_folder_to_sample_dict(self):
        """
        Takes a folder of txt files where each file corresponds to 1 spectra, encoded as m/z - intensity pairs.
        Returns a sample dictionary where each entry is 1 spectra, named according to filenames.

        example data file:
        /Users/mitchsyberg-olsen/Library/CloudStorage/Box-Box/lab/data/2025/12/2025-12-05_new_substrate_test/data/2025_12_06_0001_0_E1_1.txt
        """

        for root, dirs, files in os.walk(self.data_folder):
            root = root
            files = files

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

    def collect_mois(self) -> list:
        """
        Collects all MOIs identified in experimental samples, removes duplicates,
        and returns them as an alphabetically sorted list.
        """
        moi_list = []
        for s in self.samples.values():
            moi_list.extend(s.mois)

        unique_mois = list(set(moi_list))
        sorted_mois = sorted(unique_mois, key=lambda x: x.name)

        return sorted_mois

    def scatter(self, x, y, regression=False) -> plt.axes:
        """
        Creates a scatter plot of any two variables across the entire experiment.

        """
        pass    # stub

    @staticmethod
    def get_nested_attr(object, attr):
        """Helper function for show_sample_positions. Takes a nested attribute string and finds the value.
        i.e. 'enzyme.id' will return the ID of a particular enzyme."""

        attr_list = attr.split(sep='.')

        try:
            # recursively travel through nested objects
            for i, attr in enumerate(attr_list):
                attr_value = getattr(object, attr)
                object = attr_value

            return attr_value

        except AttributeError:
            print(f"Unable to access attribute {attr}. Skipping.")
            return None

    def show_sample_positions(self, attrs: list = ['enzyme.id', 'ntp'], fontsize=8) -> plt.axes:
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
                        attr_value = self.get_nested_attr(sample, attr)

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
                            for moi in sample.mois:
                                if moi.name == variable:
                                    d = moi.mai_intensity(sample)
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
        else:
            pass

        ax.set_title(title)
        fig.tight_layout()
        return fig, ax

    def stacked_plot(self, wells: list, xlim=None, title=None, filtered=False, smoothed=False, figsize=None,
                     label=True, label_first_only=True, overlay=False, hide_excluded=True, label_peaks=False, custom_colours=None, base_colour='#1f77b4') -> plt.axes:
        """
        Takes a list of well IDs and returns a plot of each spectrum stacked vertically.

        X-values are determined by xlim parameter. If none provided, use the whole spectrum.
        Y-values are background corrected and always relative.
        """
        plt.style.use('default')
        if figsize is None:
            if overlay:
                figsize = (12, 3)
            if not overlay:
                if len(wells) > 1:
                    figsize = (12, len(wells) * 1)  # (x size, y size)
                else:
                    figsize = (12, 3)

        fig, axs = plt.subplots(nrows=1 if overlay else len(wells),
                                ncols=1,
                                figsize=figsize,  # width, height in inches
                                sharex=True, sharey=True)

        for n, well in enumerate(wells):
            if hide_excluded is True and well in self.exclude:
                continue

            sample = self.samples[well]
            mz_plot = sample.mz  # initial mz
            i_plot = sample.i # initial i

            if filtered:
                i_plot = sample.i_filtered

            if smoothed:
                i_plot = sample.savitzky_golay(i=i_plot)

            if xlim:  # slice by xlim
                mz_plot, i_plot = sample.slice_spectrum(start=xlim[0], end=xlim[1], mz=mz_plot, i=i_plot)


            try:
                i_max = max(i_plot)  # scale everything relative to what is in range
                i_plot = [i * 100 / i_max for i in i_plot]
                if overlay:
                    # many overlapping peaks that are too uniform becomes hard to read.
                    # this section will randomly scale down each spectrum to create some variation
                    small_offset = random.uniform(0.90, 1.00)
                    i_plot = [i * small_offset for i in i_plot]
            except ZeroDivisionError:
                pass

            try:
                axis = axs[n]
            except TypeError:
                axis = axs

            if label_peaks:
                label = False   # turn of MOI species labels
                plt.subplots_adjust(hspace=0.5) # add some extra space between plots
                peaks, _ = find_peaks(i_plot, height=20)    # simple peak finding algorithm, only detects peaks >20% max signal
                for i, peak in enumerate(peaks):
                    try:    # catches an error that occurs when checking the final peak in the list
                        if peaks[i + 1] - peaks[i] < 5: # skip peaks that are less than 5 daltons separated
                            plot_peak = False
                        else:
                            plot_peak = True

                    except IndexError:
                        plot_peak = True

                    if plot_peak:
                        axis.plot(mz_plot[peak], i_plot[peak], 'x', color='black')
                        axis.vlines(mz_plot[peak], ymin=-9, ymax=i_plot[peak], colors='grey')
                        axis.text(x=mz_plot[peak], y=i_plot[peak] + 10, s=str(int(round(mz_plot[peak], 0))),
                                  horizontalalignment='center', rotation=75, verticalalignment='bottom',
                                  fontsize=10)

            if overlay:
                axis.plot(mz_plot, i_plot, linewidth=1.5)
            else:
                axis.plot(mz_plot, i_plot, linewidth=1.5, c=base_colour)

            axis.set_ylim(-10, 150)
            axis.set_xlim(min(mz_plot), max(mz_plot))
            axis.spines['top'].set_visible(False)
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_visible(False)
            axis.spines['bottom'].set_position(('data', -5))

            if not overlay:
                axis.text(x=min(mz_plot), y=0, s=f"{sample.name} ", horizontalalignment='right')

            for moi in sample.mois:  # colour mois:
                moi_range = moi.iso_dist_range
                mz_moi, i_moi = sample.slice_spectrum(start=moi_range[0], end=moi_range[1], mz=mz_plot, i=i_plot)
                if not overlay:
                    axis.plot(mz_moi, i_moi, c='#d1495b', linewidth=1.5)

                if label:
                    if overlay:
                        moi_label_height = 105
                    else:
                        try:
                            moi_label_height = max(i_moi) + 20
                        except ValueError:
                            moi_label_height = 20

                    if label_first_only:
                        if n == 0:
                            if min(mz_plot) < moi.monoisotopic_mass < max(mz_plot):  # Add species labels
                                axis.text(x=moi.monoisotopic_mass, y=moi_label_height, s=moi.name, horizontalalignment='center',
                                          fontsize=10)
                    else:
                        if min(mz_plot) < moi.monoisotopic_mass < max(mz_plot):  # Add species labels
                            axis.text(x=moi.monoisotopic_mass, y=moi_label_height, s=moi.name, horizontalalignment='center',
                                      fontsize=10)

            try:
                for colour, c_range in custom_colours.items():
                    mz_colour, i_colour = sample.slice_spectrum(start=c_range[0], end=c_range[1], mz=mz_plot, i=i_plot)
                    axis.plot(mz_colour, i_colour, c=colour, linewidth=1.5)
            except:
                pass

        plt.xlabel(f'm/z')
        plt.yticks([])

        if title:
            plt.suptitle(title)
        # plt.tight_layout()
        # plt.subplots_adjust(wspace=0, hspace=0)

        return fig, axs

    # claude /start
    # TODO: verify this code and integrate it into stacked_plot
    def stacked_plot_segmented(self, wells: list, xlim=None, title=None, filtered=False, smoothed=False,
                               figsize=None, label=True, label_first_only=True, overlay=False,
                               hide_excluded=True, label_peaks=False, custom_colours=None,
                               base_colour='#1f77b4') -> plt.axes:
        """
        Identical to stacked_plot() but uses a segmented drawing strategy to eliminate
        overlapping line artifacts. Each point on the spectrum is assigned exactly one colour
        according to priority: base_colour < moi colour < custom_colour. Contiguous runs of
        the same colour are then drawn as a single ax.plot() call, so no pixel is painted twice.

        custom_colours format: {(start, end): colour_str, ...}
            Keys are (start_mz, end_mz) tuples; values are colour strings. This allows the
            same colour to be reused across multiple ranges, and guarantees each range has
            exactly one colour. Example:
                custom_colours = {(4830, 4843): '#0072B2',
                                  (4847, 4856): '#CC6677',
                                  (4856, 4865): '#CC6677'}  # same colour, two ranges

        All other parameters are identical to stacked_plot().
        """

        def build_colour_array(mz_arr, moi_ranges, custom_colour_ranges, base_colour):
            """
            Build a list of colours, one per m/z point, respecting priority:
                1. base_colour       (lowest priority)
                2. moi colour        ('#d1495b')
                3. custom colours    (highest priority)

            moi_ranges:           list of (start, end) tuples
            custom_colour_ranges: dict of {(start, end): colour_str}
            """
            colours = [base_colour] * len(mz_arr)

            # Layer 2: MOI regions
            for start, end in moi_ranges:
                for idx, mz_val in enumerate(mz_arr):
                    if start <= mz_val <= end:
                        colours[idx] = '#d1495b'

            # Layer 3: custom colours (highest priority)
            if custom_colour_ranges:
                for (start, end), colour in custom_colour_ranges.items():
                    for idx, mz_val in enumerate(mz_arr):
                        if start <= mz_val <= end:
                            colours[idx] = colour

            return colours

        def plot_segmented(axis, mz_arr, i_arr, colours):
            """
            Walk the colour array and emit one ax.plot() call per contiguous
            same-colour run. Adjacent segments share their boundary point so
            there are no gaps and no overlapping draws.
            """
            mz_arr = list(mz_arr)
            i_arr = list(i_arr)
            n = len(mz_arr)
            if n == 0:
                return

            seg_start = 0
            for idx in range(1, n):
                if colours[idx] != colours[seg_start]:
                    # Plot the segment, including the current point as a shared
                    # boundary so adjacent segments connect seamlessly
                    axis.plot(mz_arr[seg_start:idx + 1],
                              i_arr[seg_start:idx + 1],
                              c=colours[seg_start], linewidth=1.5)
                    seg_start = idx  # next segment starts at the boundary point

            # Plot the final segment
            axis.plot(mz_arr[seg_start:],
                      i_arr[seg_start:],
                      c=colours[seg_start], linewidth=1.5)

        # ------------------------------------------------------------------ #
        #  Figure setup — identical to stacked_plot()                         #
        # ------------------------------------------------------------------ #
        plt.style.use('default')
        if figsize is None:
            if overlay:
                figsize = (12, 3)
            else:
                figsize = (12, len(wells) * 1) if len(wells) > 1 else (12, 3)

        fig, axs = plt.subplots(nrows=1 if overlay else len(wells),
                                ncols=1,
                                figsize=figsize,
                                sharex=True, sharey=True)

        for n, well in enumerate(wells):
            if hide_excluded and well in self.exclude:
                continue

            sample = self.samples[well]
            mz_plot = sample.mz
            i_plot = sample.i

            if filtered:
                i_plot = sample.i_filtered
            if smoothed:
                i_plot = sample.savitzky_golay(i=i_plot)
            if xlim:
                mz_plot, i_plot = sample.slice_spectrum(start=xlim[0], end=xlim[1],
                                                        mz=mz_plot, i=i_plot)

            try:
                i_max = max(i_plot)
                i_plot = [i * 100 / i_max for i in i_plot]
                if overlay:
                    small_offset = random.uniform(0.90, 1.00)
                    i_plot = [i * small_offset for i in i_plot]
            except ZeroDivisionError:
                pass

            try:
                axis = axs[n]
            except TypeError:
                axis = axs

            # ---------------------------------------------------------- #
            #  Collect colour regions for this sample                      #
            # ---------------------------------------------------------- #
            moi_ranges = []
            for moi in sample.mois:
                if not overlay:  # MOI colouring is suppressed in overlay mode (matches original)
                    moi_ranges.append(moi.calc_iso_dist_range())

            custom_colour_ranges = custom_colours if custom_colours else {}
            # Normalise keys to plain tuples (support both list and tuple from the caller)
            custom_colour_ranges = {tuple(r): c for r, c in
                                    custom_colour_ranges.items()} if custom_colour_ranges else {}

            # ---------------------------------------------------------- #
            #  Build per-point colour array and draw segmented spectrum    #
            # ---------------------------------------------------------- #
            colours = build_colour_array(mz_plot, moi_ranges, custom_colour_ranges, base_colour)
            plot_segmented(axis, mz_plot, i_plot, colours)

            # ---------------------------------------------------------- #
            #  Axis cosmetics — identical to stacked_plot()               #
            # ---------------------------------------------------------- #
            axis.set_ylim(-10, 150)
            axis.set_xlim(min(mz_plot), max(mz_plot))
            axis.spines['top'].set_visible(False)
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_visible(False)
            axis.spines['bottom'].set_position(('data', -5))

            if not overlay:
                axis.text(x=min(mz_plot), y=0, s=f"{sample.name} ",
                          horizontalalignment='right')

            # ---------------------------------------------------------- #
            #  Peak labels                                                 #
            # ---------------------------------------------------------- #
            if label_peaks:
                plt.subplots_adjust(hspace=0.5)
                peaks, _ = find_peaks(i_plot, height=20)
                for i, peak in enumerate(peaks):
                    try:
                        plot_peak = (peaks[i + 1] - peaks[i]) >= 5
                    except IndexError:
                        plot_peak = True
                    if plot_peak:
                        axis.plot(mz_plot[peak], i_plot[peak], 'x', color='black')
                        axis.vlines(mz_plot[peak], ymin=-9, ymax=i_plot[peak], colors='grey')
                        axis.text(x=mz_plot[peak], y=i_plot[peak] + 10,
                                  s=str(int(round(mz_plot[peak], 0))),
                                  horizontalalignment='center', rotation=75,
                                  verticalalignment='bottom', fontsize=10)

            # ---------------------------------------------------------- #
            #  MOI and custom-colour region labels                         #
            # ---------------------------------------------------------- #
            if label and not label_peaks:
                for moi in sample.mois:
                    moi_range = moi.calc_iso_dist_range()
                    _, i_moi = sample.slice_spectrum(start=moi_range[0], end=moi_range[1],
                                                     mz=mz_plot, i=i_plot)
                    if overlay:
                        moi_label_height = 105
                    else:
                        try:
                            moi_label_height = max(i_moi) + 20
                        except ValueError:
                            moi_label_height = 20

                    if label_first_only:
                        if n == 0 and min(mz_plot) < moi.monoisotopic_mass < max(mz_plot):
                            axis.text(x=moi.monoisotopic_mass, y=moi_label_height, s=moi.name,
                                      horizontalalignment='center', fontsize=10)
                    else:
                        if min(mz_plot) < moi.monoisotopic_mass < max(mz_plot):
                            axis.text(x=moi.monoisotopic_mass, y=moi_label_height, s=moi.name,
                                      horizontalalignment='center', fontsize=10)

        plt.xlabel('m/z')
        plt.yticks([])
        if title:
            plt.suptitle(title)

        return fig, axs

    # claude /end

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

    def sorted_signal_plot(self, xlim=[50,100], label_noisy_samples=False) -> plt.axes:
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

    def collect_sample_data(self, headers: list, i_type='filtered') -> dict:
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

                # if that doesn't work, it could be a MOI species
                for moi in v.mois:
                    if moi.name == item:
                        data[k][item] = moi.mai_intensity(sample=v, i_type=i_type)

        return data

    def write_excel_readme(self, workbook):
        """Adds a README worksheet to the data output workbook. Use with write_to_excel()."""

        # write README
        readme = f"""
                mois: Contains all molecules of interest and their associated information.

                raw: Unprocessed data. Each entry corresponds to the ion intensity of the most abundant isotope of a given molecule.

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

    def write_excel_mois(self, workbook):
        """Adds a worksheet to the data output workbook, containing information about.
         All the molecules of interest in the experiment.
         Use with write_to_excel()."""

        worksheet = workbook.add_worksheet('mois')

        # write MOI data
        headers = ['name', 'sequence', 'composition', 'monoisotopic mass', 'average mass']
        for i, header in enumerate(headers):  # write headers
            worksheet.write(0, i, header)

        mois = self.collect_mois()
        for i, moi in enumerate(mois):  # i + 1 for row
            data = [moi.name, moi.seq, moi.composition_str(), moi.monoisotopic_mass, moi.average_mass]
            for j, ele in enumerate(data):  # j for column
                worksheet.write(i + 1, j, ele)

    def collect_attributes_for_excel(self):
        """Collects all attributes from samples that should be written to excel."""

        # collect regular attributes
        attribute_keys = []
        for s in self.samples.values():
            for k, v in s.collect_attributes().items():
                if k not in attribute_keys:
                    attribute_keys.append(k)

        # collect moi names
        moi_attributes = [m.name for m in self.collect_mois()]

        return attribute_keys + moi_attributes

    def write_excel_data(self, workbook):
        """Adds worksheets to the data output workbook, containing all experimental data."""
        # write data sheets
        headers = self.collect_attributes_for_excel()

        i_types = ['raw', 'bg_sub', 'filtered']
        for i_type in i_types:
            worksheet = workbook.add_worksheet(i_type)  # create a worksheet
            for i, header in enumerate(headers):  # write headers
                worksheet.write(0, i, header)

            # write data
            i = 1
            data = self.collect_sample_data(headers=headers, i_type=i_type)  # collect data according to i_type param

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
        self.write_excel_readme(workbook)
        self.write_excel_mois(workbook)
        self.write_excel_data(workbook)
        workbook.close()

