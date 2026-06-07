import matplotlib.pyplot as plt
import numpy as np
import random
from . import helper
from . import analytes


class Sample:
    def __init__(self, file, well, mz, i, noise_cutoff: float, chip=0, mz_offset=0):
        self.file = file  # source data file
        self.chip = chip    # for AnchorChip plates, chip=0 refers to regular sample spots and chip=1 refers to calibrant spots
        self.well = well
        self.name = self.well   # default naming but this can be changed
        self.noise_cutoff = noise_cutoff
        self.mz_offset = mz_offset  # used to manually adjust the spectra if calibration is off.
        self.analytes = []

        # process spectrum
        self.mz_raw = mz    # Keeps a record of the original mz values
        self.mz = mz    # can be updated using the recalc_mz() method
        self.i = i
        self.i_bg_subtracted = []
        self.i_filtered = []

        self.background = min(self.i)
        self._i_background_subtracted()  # populates self.i_bg_subtracted

        self.noise = self._calculate_noise(self.noise_cutoff)
        self._filter_i()   # populates self.i_filtered

        # params for plotting
        self.total_ion = sum(self.i)
        self.max_i = max(self.i)
        self.num_points = len(self.i)

    def recalc_mz(self):
        """Updates self.mz, useful when mz_offset has been changed."""
        self.mz = [point + self.mz_offset for point in self.mz_raw]

    def _collect_attributes(self):
        """
        Collects all attributes of the Sample object and returns them as a dictionary.
        Excludes attributes that are only meant to be seen internally.
        """
        include = ['file', 'chip', 'well', 'name', 'noise_cutoff', 'mz_offset', 'background', 'noise']
        exclude = ['misc_conditions', 'analytes', 'mz_raw', 'mz', 'i', 'i_bg_subtracted', 'i_filtered', 'total_ion',
                   'max_i', 'num_points']

        attrs = {}

        for key in include: # first, add essential attributes, important for ordering
            attrs[key] = getattr(self, key)

        for k, v in self.__dict__.items():  # then iterate through all attributes
            if k in include:
                pass    # already added
            elif k not in exclude:
                attrs[k] = v

        return attrs

    def analyte_intensities(self) -> list[tuple]:
        """
        Returns a list of tuples: (molecule name, intensity value), one for each molecule of interest.
        The intensity value is determined by Oligo.mai_intensity()
        """
        intensities = []
        for a in self.analytes:
            intensities.append(
                (a.name, a._mai_intensity(sample=self))
            )
        return intensities

    def _calculate_noise(self, cutoff: float) -> int:
        """
        Calculates the noise value of the spectra.
        For our purposes, this is a constant value rather than a function.
        """
        i_vals = self.i_bg_subtracted
        lowest_vals = sorted(i_vals)[:int(len(i_vals)*cutoff)]
        return max(lowest_vals)

    def _i_background_subtracted(self) -> None:
        """
        Subtracts the background signal from the spectrum.
        """
        background = self.background
        self.i_bg_subtracted = [i - background for i in self.i]

    def _filter_i(self) -> None:
        """
        Keeps only intensity values above the noise of the instrument.
        """

        i_filtered = []
        for i in self.i_bg_subtracted:
            if i > self.noise:
                i_filtered.append(i)
            else:
                i_filtered.append(0)

        self.i_filtered = i_filtered

    def plot_distance_between_points(self) -> plt.axes:
        """
        Generates a plot of the distance (daltons) between each data point collected in the spectrum.
        """
        delta_mz = []
        i = 1
        while i < len(self.mz):
            delta_mz.append(self.mz[i] - self.mz[i-1])
            i += 1

        fig, ax = plt.subplots(figsize=(16, 3))
        plt.plot(self.mz[:-1], delta_mz)
        plt.title('Distance between adjacent mz points in spectrum data')
        plt.ylabel(f'Distance (daltons)')
        plt.xlabel('m/z')

        return fig, ax


    # TODO: make this the main plot that gets called by Experiment
    def sorted_signal_plot(self) -> plt.axes:
        """
        Generates # TODO finish describing.
        """
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(8, 3))

        # Collect all the data
        x = np.linspace(0, 100, num=len(self.i))
        y = sorted(self.i)
        plt.plot(x,y)

        y_at_cutoff = y[round(self.noise_cutoff*len(y))]
        # plt.text(x=self.noise_cutoff*100, y=y_at_cutoff+1, s=f" {self.noise}")

        plt.xlim(50, 100)
        # plt.xticks(np.linspace(50, 100, num=6))
        # plt.yticks(np.linspace(0, 100, num=11))

        plt.axvline(x=self.noise_cutoff*100)

        plt.title(f"Signal intensity distribution")
        plt.xlabel(f'Percent of data')
        plt.ylabel(f'Relative signal')

        return fig, ax

    def plot_analyte(self, analyte: analytes.Analyte, filtered=True) -> plt.axes:
        """
        Generates a plot for an analyte, zoomed in to visualize the isotope distribution and mass accuracy.
        """
        return self.plot(xlim=analyte.iso_dist_range,
                         relative=True,
                         theoretical_dist=True,
                         filtered=filtered,
                         label_analytes=False,
                         title=f"{analyte.name}, monoisotopic mass = {analyte.monoisotopic_mass}, charge = +{analyte.charge}")



    def _base_plotter(
            self,
            ax=plt.axes,
            xlim=None,
            filtered=False,
            smoothed=False,
            normalized=False,
            overlay=False,
            highlight_analytes=True,
            label_analytes=True,
            label_peaks=False,
            custom_colours=None,
            base_colour='#1f77b4',
            analyte_colour='#d1495b',
            linewidth=1.5
    ) -> None:
        """
        base function for generating and mz vs i plot. Wrapped by other functions.
        """

        ### Collect correct mz & i values
        mz_plot = self.mz
        i_plot = self.i

        if filtered:    # i is now filtered based on noise cutoff
            i_plot = self.i_filtered

        if xlim:    # adjust x axis limits
            mz_plot, i_plot = helper._slice_spectrum(start=xlim[0],
                                                   end=xlim[1],
                                                   mz=mz_plot,
                                                   i=i_plot)

        if smoothed:    # apply a savitzky-golay filter
            i_plot = helper.savitzky_golay(i=i_plot)

        if normalized:  # normalize to 100 as maximum value
            try:
                i_plot = [i*100 / max(i_plot) for i in i_plot]
            except ZeroDivisionError:
                print("Normalization failed, all intensity values are 0. \n"
                      "Try plotting your data without normalization.")

        if overlay: # apply a random scaling to prevent overlapping lines from becoming unreadable
            i_plot = [i * random.uniform(0.90, 1.00) for i in i_plot]



        ### Assign a colour to every data point
        if custom_colours:
            # make sure the custom colours are in the correct format
            custom_colours = {tuple(x_range): colour for x_range, colour in custom_colours.items()}
        else:
            custom_colours = {}

        if highlight_analytes:
            analyte_ranges = [a.iso_dist_range for a in self.analytes]
        else:
            analyte_ranges = []

        # build the colour array
        colours = helper._build_colour_array(x_vals=mz_plot,
                                             analyte_ranges=analyte_ranges,
                                             custom_colour_ranges=custom_colours,
                                             base_colour=base_colour,
                                             analyte_colour=analyte_colour)

        ### Plot the line
        helper._plot_segmented(axis=ax,
                               mz=mz_plot,
                               i=i_plot,
                               colour_array=colours,
                               linewidth=linewidth)

        ### Add labels and annotations
        if label_peaks:
            peaks, _ = helper._call_peaks(i_plot)
            for peak in peaks:
                ax.plot(mz_plot[peak], i_plot[peak], 'x', color='black')
                ax.vlines(mz_plot[peak], ymin=-9, ymax=i_plot[peak], colors='grey')
                ax.text(
                    x=mz_plot[peak], y=i_plot[peak] + 10,
                    s=str(int(round(mz_plot[peak], 0))),
                    horizontalalignment='center', rotation=75,
                    verticalalignment='bottom', fontsize=10,
                )

        if label_analytes:
            for a in self.analytes:
                if min(mz_plot) < a.monoisotopic_mass < max(mz_plot):
                    a_range = a.iso_dist_range
                    _, i_slice = helper._slice_spectrum(start=a_range[0],
                                                        end=a_range[1],
                                                        mz=mz_plot,
                                                        i=i_plot)

                    if overlay:
                        analyte_label_height = 105
                    elif normalized:
                        analyte_label_height = max(i_slice) + 10
                    else:
                        y_min, y_max = ax.get_ylim()
                        analyte_label_height = max(i_slice) + (y_max - y_min) * 0.10

                    ax.text(x=a.monoisotopic_mass,
                            y=analyte_label_height,
                            s=a.name,
                            horizontalalignment='center',
                            fontsize=10)


        ### Formatting
        ymax = round(max(i_plot) * 1.2 + 1)  # +1 accounts for possibly all 0 values
        ax.set_ylim(ymin=min(i_plot), ymax=ymax)
        ax.set_xlim(min(mz_plot), max(mz_plot))


    def plot(
            self,
            xlim=None,
            normalized=False,
            label_peaks=False,
            theoretical_dist=False,
            highlight_analytes=True,
            label_analytes=True,
            filtered=False,
            smoothed=False,
            title=None,
            custom_colours=None,
            base_colour='#1f77b4',
            analyte_colour='#d1495b',
            theoretical_colour='#CC6677',
            linewidth=1.5,
            ax=None,
            figsize=(8, 3)
    ) -> plt.axes:

        """
        General purpose plotting function for a Sample's spectrum.
        """
        ### Settings
        if theoretical_dist:
            normalized = True

        if not title:
            title = f'Sample: {self.name}'

        if normalized:
            ylabel = 'Relative intensity (au)'
        else:
            ylabel = 'Ion count'

        if ax is None:  # create a standalone figure
            plt.style.use('default')
            fig, ax = plt.subplots(figsize=figsize)

        # Base plot
        self._base_plotter(ax=ax,
                           xlim=xlim,
                           filtered=filtered,
                           smoothed=smoothed,
                           normalized=normalized,
                           overlay=False,   # Always false for standalone plots
                           highlight_analytes=highlight_analytes,
                           label_analytes=label_analytes,
                           label_peaks=label_peaks,
                           custom_colours=custom_colours,
                           base_colour=base_colour,
                           analyte_colour=analyte_colour,
                           linewidth=linewidth)

        if theoretical_dist:   # Overlay theoretical isotope distributions
            for a in self.analytes:
                a.iso_dist_plot(ax=ax, y_max=100, colour=theoretical_colour, annotate=False)

            # Create a legend
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys())

        # Formatting
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel('m/z')

        return ax

