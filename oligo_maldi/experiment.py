import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import numpy as np
import xlsxwriter
from string import ascii_uppercase
import os
import threading
from .sample import Sample
from .helper import *

"""
TODO
plots should automatically include processing in title (bg subtracted?)
Salt to be included in analysis?


    - plotting function to do scatter of any two variables against each other
    - interactive plot to choose noise_cutoff
    - More intelligent signal processing
        - adjust signal based on expected distribution
        - differentiate between noise and low signal based on distribution
    - increase the speed of the report() function
        - graphing faster via multithreading?
"""

class Experiment:
    def __init__(self, samples, noise_cutoff=0.90):
        # self.filename = xml_file
        # self.xml_tree = self.xml_to_tree(xml_file)
        self.noise_cutoff = noise_cutoff
        # self.samples = self.trees_to_sample_dict()
        self.samples = self.reinitiate_samples(samples)

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

    def xml_to_tree(self, xml_file: str) -> list:
        """
        takes one xml file and converts it to a data tree.
        """
        tree = ET.parse(xml_file)
        root = tree.getroot()
        return root

    def trees_to_sample_dict(self) -> dict:
        """
        converts data tree derived from a MALDI xml file to a dictionary
        key = sample well ID
        value = Sample object
        """
        sample_dict = {}
        # tree = root -> [2] = analysis -> [1] = DataAnalysis -> [1:] = list of ms_spectrum
        ms_spectra = self.xml_tree[2][1][1:]
        for spectrum in ms_spectra:
            sample = Sample(spectrum,
                            noise_cutoff=self.noise_cutoff)  # instantiate a Sample object using the data from the ms_spectrum
            sample_dict[sample.well] = sample  # add to dictionary, key=well, value=Sample object

        return sample_dict

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
        pass

    def heatmap(self, numerator, denominator=None) -> plt.axes:
        """
        Creates a heatmap of a variable across the 386 well sample plate.
        Possible variables:
        Sample total_ion, noise, background, or any MOI.
        If an MOI is indicated, the value for moi.mai_intensity() for that MOI will be plotted.

        If only numerator is specified, the value of that variable will be plotted.
        If numerator and denominator are specified, the ratio of the two will be plotted.
        """
        rows = list(ascii_uppercase[:16])   # ABCDEFGHIJKLMOP
        columns = [str(i) for i in range(1, 25)]    # 1-24

        data = []
        for row in rows:
            row_data = []
            for col in columns:
                key = row + col # ie 'C12'
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
        im = ax.imshow(data)
        cbar = ax.figure.colorbar(im, ax=ax)

        # Show all ticks and label them with the respective list entries
        ax.set_xticks(range(len(columns)), labels=columns)
        ax.set_yticks(range(len(rows)), labels=rows)

        if denominator:
            title = f"{numerator} / {denominator}"
        else:
            title = numerator

        ax.set_title(title)
        fig.tight_layout()
        return fig, ax

    def stacked_plot(self, wells: list, xlim=None, title=None, filtered=True) -> plt.axes:
        """
        Takes a list of well IDs and returns a plot of each spectrum stacked vertically.

        X-values are determined by xlim parameter. If none provided, use the whole spectrum.
        Y-values are background corrected and always relative.
        """

        plt.style.use('default')
        if len(wells) > 1:
            figsize = (8, len(wells) * 1.5)
        else:
            figsize = (8, 3)

        fig, axs = plt.subplots(nrows=len(wells), ncols=1,
                                figsize=figsize,  # width, height in inches
                                sharex=True, sharey=True)

        for n, well in enumerate(wells):
            sample = self.samples[well]
            mz_plot = sample.mz  # initial mz
            if filtered:
                i_plot = sample.i_filtered  # initial i
            else:
                i_plot = sample.i

            if xlim:  # slice by xlim
                mz_plot, i_plot = sample.slice_spectrum(start=xlim[0], end=xlim[1], mz=mz_plot, i=i_plot)

            try:
                i_max = max(i_plot)  # scale everything relative to what is in range
                i_plot = [i * 100 / i_max for i in i_plot]
            except ZeroDivisionError:
                pass

            try:
                axis = axs[n]
            except TypeError:
                axis = axs

            axis.plot(mz_plot, i_plot, linewidth=1)
            axis.set_ylim(0, 125)
            axis.set_xlim(min(mz_plot), max(mz_plot))
            axis.text(x=max(mz_plot), y=110, s=f"{sample.name} ", horizontalalignment='right')

            for moi in sample.mois:  # colour mois:
                moi_range = moi.iso_dist_range()
                mz_moi, i_moi = sample.slice_spectrum(start=moi_range[0], end=moi_range[1], mz=mz_plot, i=i_plot)
                axis.plot(mz_moi, i_moi, c='#d1495b', linewidth=1)
                try:
                    moi_label_height = max(i_moi) + 10
                except ValueError:
                    moi_label_height = 10

                if min(mz_plot) < moi.monoisotopic_mass < max(mz_plot):  # Add species labels
                    axis.text(x=moi.monoisotopic_mass, y=moi_label_height, s=moi.name, horizontalalignment='center',
                              fontsize=10)

        fig.text(0, 0.5, f'Relative intensity (au)', ha='left', va='center', rotation='vertical')
        plt.xlabel(f'm/z')
        if title:
            plt.suptitle(title)
        plt.tight_layout()
        plt.subplots_adjust(wspace=0, hspace=0)

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

    def sorted_signal_plot(self, xlim=[50,100]) -> plt.axes:
        """
        Generates .
        """

        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 3))

        # Collect all the data
        for well, sample in self.samples.items():
            x = np.linspace(0, 100, num=len(sample.i))
            y_max = max(sample.i)
            y = sorted([i * 100 / y_max for i in sample.i])
            plt.plot(x, y)

            y_at_cutoff = y[round(self.noise_cutoff * len(y))]
            if y_at_cutoff > 10:
                plt.text(x=self.noise_cutoff * 100, y=y_at_cutoff + 1, s=f" {well}")

        plt.xlim(xlim)
        plt.yticks(np.linspace(0, 100, num=11))

        plt.axvline(x=self.noise_cutoff * 100)
        plt.text(x=self.noise_cutoff * 100, y=90, s=" Noise cutoff")
        plt.text(x=self.noise_cutoff * 100, y=50, s=" Samples w/\n S/N<10")

        plt.title(f"Per sample signal intensity distribution")
        plt.xlabel(f'Percent of data')
        plt.ylabel(f'Relative signal')

        return fig, ax

    def plot_wells(self, wells: dict, folder_path: str):
        """Helper function for report()."""

        # Individual well plots
        for k, v in wells.items():
            subfolder = f"{folder_path}/wells/{k}"
            os.makedirs(subfolder)

            # overall
            fig, ax = v.plot(filtered=False)
            figsave_fast(fig, f"{subfolder}/raw.png")
            # plt.savefig(f"{subfolder}/raw.png")
            plt.close()

            fig, ax = v.plot(filtered=True)
            figsave_fast(fig, f"{subfolder}/filtered.png")
            # plt.savefig(f"{subfolder}/filtered.png")
            plt.close()

            fig, ax = v.plot_distance_between_points()
            figsave_fast(fig, f"{subfolder}/distance_between_points.png")
            # plt.savefig(f"{subfolder}/distance_between_points.png")
            plt.close()

            fig, ax = v.sorted_signal_plot()
            figsave_fast(fig, f"{subfolder}/sorted_signal_plot.png")
            # plt.savefig(f"{subfolder}/sorted_signal_plot.png")
            plt.close()

            subsubfolder = f"{subfolder}/mois"
            os.makedirs(subsubfolder)
            for moi in v.mois:
                try:
                    fig, ax = v.plot_moi(moi)
                    figsave_fast(fig, f"{subsubfolder}/{moi.name}.png")
                    # plt.savefig(f"{subsubfolder}/{moi.name}.png")
                    plt.close()
                except ValueError:  # occurs when there is no data at a particular MOI region
                    pass

    def report(self, threads=4):
        """
        General:
            - TODO: implement multithreading to speed up plotting.
            - heatmaps:
                total ion, noise, background
                TODO: signal to noise ratio for each MOI
            - sorted signal plot
        MOIs:
            - theoretical isotope distribution
        Wells:
            - overall spectra:
                - raw, bg_sub, filtered
            - moi overlays
            - signal intensity distribution
            - distance between data
        """
        root_dir = os.path.dirname(self.filename)
        folder_path = root_dir + "/report"

        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        else:
            raise FileExistsError(f"Folder '{folder_path}' already exists.")

        # Turn interactive plotting off
        plt.ioff()
        plt.close('all')

        # Make heatmaps
        fig, ax = self.heatmap('total_ion')
        figsave_fast(fig, folder_path + "/ion_totals.png")
        # plt.savefig(folder_path + "/ion_totals.png")
        plt.close()
        fig, ax = self.heatmap('noise')
        figsave_fast(fig, folder_path + "/noise.png")
        # plt.savefig(folder_path + "/noise.png")
        plt.close()
        fig, ax = self.heatmap('background')
        figsave_fast(fig, folder_path + "/background.png")
        # plt.savefig(folder_path + "/background.png")
        plt.close()

        # Make signal plot
        fig, ax = self.sorted_signal_plot()
        figsave_fast(fig, folder_path + "/sorted_signal.png")
        # plt.savefig(folder_path + "/sorted_signal.png")
        plt.close()

        # Make MOI plots
        subfolder = f"{folder_path}/mois"
        os.makedirs(subfolder)
        for moi in self.collect_mois():
            fig, ax = moi.iso_dist_plot()
            figsave_fast(fig, f"{subfolder}/{moi.name}.png")
            # plt.savefig(f"{subfolder}/{moi.name}.png")
            plt.close()

        # Individual well plots
        self.plot_wells(wells=self.samples, folder_path=folder_path)

        # Turn interactive plotting back on
        plt.ion()

    def collect_sample_data(self, i_type='filtered') -> dict:
        """
        Collects general sample information and the most abundant ion intensity for each sample.
        """

        data = {}
        for k, v in self.samples.items():
            data[k] = {'well': v.well,
                       'name': v.name,
                       'background': v.background,
                       'noise': v.noise,
                       'noise_cutoff': v.noise_cutoff,
                       'total_ion_count': sum(v.i)}

            for moi in v.mois:
                data[k][moi.name] = moi.mai_intensity(sample=v, i_type=i_type)

        return data

    def write_to_excel(self) -> None:
        """
        Writes data to an excel file.
        """

        # create file
        outfile = self.filename.replace(".xml", ".xlsx")
        workbook = xlsxwriter.Workbook(outfile)
        worksheet0 = workbook.add_worksheet('readme')
        worksheet1 = workbook.add_worksheet('mois')

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
        worksheet0.insert_textbox('B2', readme,
                                  {'width': 1000, 'height': 1000})
        # worksheet0.write(0, 0, readme)

        # write MOI data
        headers = ['name', 'sequence', 'composition', 'monoisotopic mass', 'average mass']
        for i, header in enumerate(headers):  # write headers
            worksheet1.write(0, i, header)

        mois = self.collect_mois()
        for i, moi in enumerate(mois):  # i + 1 for row
            data = [moi.name, moi.seq, moi.composition_str(), moi.monoisotopic_mass, moi.average_mass]
            for j, ele in enumerate(data):  # j for column
                worksheet1.write(i + 1, j, ele)

        # write data sheets
        headers = ['well', 'name', 'background', 'noise', 'noise_cutoff', 'total_ion_count'] + [m.name for m in
                                                                                                self.collect_mois()]

        i_types = ['raw', 'bg_sub', 'filtered']
        for i_type in i_types:
            worksheet = workbook.add_worksheet(i_type)  # create a worksheet
            for i, header in enumerate(headers):  # write headers
                worksheet.write(0, i, header)

            # write data
            i = 1
            data = self.collect_sample_data(i_type=i_type)  # collect data according to i_type param
            for k, v in data.items():  # iterate through samples
                for j, header in enumerate(headers):  # iterate through headers
                    try:
                        worksheet.write(i, j, v[header])  # if data matches header, write
                    except KeyError:
                        worksheet.write(i, j, "n/a")  # if no data exists for this header, write n/a
                i += 1

        workbook.close()

