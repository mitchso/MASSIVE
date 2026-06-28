import matplotlib.pyplot as plt
import matplotlib.axes
import numpy as np
import random
from . import helper
from . import analytes
from typing import Literal

# TODO: implement actual baseline subtraction algorithms
# TODO: implement actual peak identification algorithms
# TODO: improve peak intensity calling
# TODO: implement peak area calculation and visualization

class Sample:
    """
    Represents a single MALDI-ToF spectrum.

    `Sample` can be assigned a list of [`Analytes`][MASSIVE.analytes.Analyte], which then enables each Analyte to be automatically detected and quantified.

    Additionally, `Sample` objects can be assigned to a [`Collection`][MASSIVE.collections.Collection], which then enables the groups of samples to be visualized and analyzed together.

    When creating a [`Collection`][MASSIVE.collections.Collection], `Sample` objects are generated automatically from the data files.
    This is the typical workflow, but it is also possible to manually create `Sample` objects if you wish.
    """
    def __init__(self, file:str, id:str, mz:list[float], i:list[float], noise_cutoff: float=0.95, chip:Literal[0,1]=0, mz_offset:float=0):
        """

        Args:
            file: Data file path (see [Input data formatting assistance](../guides/data_formatting.md))
            id: Unique sample identifier. For samples on a [`Plate`][MASSIVE.collections.Plate], this should be the position of the sample on the plate (e.g. 'A1')
            mz: m/z values of the spectrum, encoded as a list of floats.
            i: Intensity values of the spectrum, encoded as a list of floats.
            noise_cutoff: The noise cutoff value used to filter out low-signal data points. At 0.95, every intensity value along the spectrum except for the highest 5% is set to 0. This can be tuned to eliminate noise at low signal intensities.
            chip: Bruker MALDI-ToF plates are divided into two 'chips', where 0 refers to regular sample spots and 1 refers to calibrant spots.
            mz_offset: Used to manually adjust the spectra if calibration is off. This adjustment simply slides the spectra by the given value along the x-axis.

        Attributes:
            file (str): Data file path (see [Input data formatting assistance](../guides/data_formatting.md))
            chip (int): Bruker MALDI-ToF plates are divided into two 'chips', where 0 refers to regular sample spots and 1 refers to calibrant spots.
            id (str): Unique sample identifier. For samples on a [`Plate`][MASSIVE.collections.Plate], this should be the position of the sample on the plate (e.g. 'A1')
            name (str): Human-readable name of the sample (defaults to `self.id`, can be changed by the user)
            noise_cutoff (float): The noise cutoff value used to filter out low-signal data points.
            background (float): Background signal level of the spectrum.
            noise (float): Noise of the spectrum (anything below this value is considered noise).
            analytes (list): List of [`Analytes`][MASSIVE.analytes.Analyte] to search for in the sample.
            mz_raw (list): Original m/z values of the spectrum, encoded as a list of floats. Unaffected by `mz_offset`.
            mz (list): Adjusted m/z values of the spectrum, encoded as a list of floats. Adjusted by `mz_offset`. If `mz_offset` is changed after initializing the object, you must call `self.recalc_mz()` to update this attribute.
            i (list): Intensity values of the spectrum, encoded as a list of floats.
            i_bg_subtracted (list): Background subtracted intensity values of the spectrum, encoded as a list of floats.
            i_filtered (list): Intensity values of the spectrum, anything below `noise` set to 0, encoded as a list of floats.
        """

        self.file = file  # source data file
        self.chip = chip    # for AnchorChip plates, chip=0 refers to regular sample spots and chip=1 refers to calibrant spots
        self.id = id
        self.name = self.id   # default naming but this can be changed
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

    def recalc_mz(self):
        """Updates self.mz to align with the current `self.mz_offset` value."""
        self.mz = [point + self.mz_offset for point in self.mz_raw]
        return True

    def _collect_attributes(self):
        """
        Collects all attributes of the Sample object and returns them as a dictionary.
        Excludes attributes that are only meant to be seen internally.
        """
        include = ['file', 'chip', 'id', 'name', 'noise_cutoff', 'mz_offset', 'background', 'noise']
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
        Returns a list of tuples: (analyte name, intensity value), one for each molecule of interest.
        The intensity value is determined by [`Analyte`][MASSIVE.analytes.Analyte]`._peak_intensity()`.
        """
        intensities = []
        for a in self.analytes:
            intensities.append(
                (a.name, a._peak_intensity(sample=self))
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

    def plot_distance_between_points(self) -> matplotlib.axes.Axes:
        """
        Generates a plot of the distance (daltons) between each data point collected in the spectrum.
        This can be revealing for understanding how your instrument creates m/z bins.
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
    def sorted_signal_plot(self) -> matplotlib.axes.Axes:
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

    def plot_analyte(self, analyte: analytes.Analyte, filtered:bool=True) -> matplotlib.axes.Axes:
        """
        Generates a plot for an analyte, zoomed in to visualize the isotope distribution and mass accuracy.
        """
        return self.plot(xlim=analyte.iso_dist_range,
                         theoretical_dist=True,
                         filtered=filtered,
                         label_analytes=False,
                         title=f"{analyte.name}, monoisotopic mass = {analyte.monoisotopic_mass}, charge = {analyte.charge}")


    def _base_plotter(
            self,
            ax=matplotlib.axes.Axes,
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
            i_plot = helper._savitzky_golay(i=i_plot)

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
            xlim:tuple=None,
            normalized:bool=False,
            label_peaks:bool=False,
            theoretical_dist:bool=False,
            highlight_analytes:bool=True,
            label_analytes:bool=True,
            filtered:bool=False,
            smoothed:bool=False,
            title:str|None=None,
            custom_colours:dict|None=None,
            base_colour:str='#1f77b4',
            analyte_colour:str='#d1495b',
            theoretical_colour:str='#CC6677',
            linewidth:float=1.5,
            ax:matplotlib.axes.Axes|None=None,
            figsize:tuple=(8, 3)
    ) -> matplotlib.axes.Axes:

        """
        Main plotting function for [`Sample`][MASSIVE.sample.Sample]. Highly customizable, includes automatic highlighting and labeling of `Analytes` contained in `Sample.analytes`.

        Each position is assigned exactly one colour by priority:
            base_colour < analyte_colour < custom_colour

        Args:
            xlim: X-axis limits.
            normalized: Converts y-axis values to relative intensities (au) between 0 and 100.
            label_peaks: Auto-detect and label peaks in the spectrum.
            theoretical_dist: Overlays the spectrum with theoretical isotope distributions for each `Analyte` in `Sample.analytes`.
            highlight_analytes: Highlights the area in the spectrum where each `Analyte` is located.
            label_analytes: Labels each `Analyte` with `Analyte.name`.
            filtered: If `True`, only plot data points above the noise cutoff. If `False`, plot all data points.
            smoothed: Apply a savitzky-golay filter to the spectrum.
            title: Optional title for the plot.
            custom_colours: Dictionary of custom colour ranges, e.g. {(100, 200): 'red', (200, 300): 'blue'}. Enables the user to bring attention to specified locations on the spectrum.
            base_colour: Base colour of the plot.
            analyte_colour: Colour for areas covered by `highlight_analytes`.
            theoretical_colour: Colour for the theoretical isotope distributions.
            linewidth: Line width of the plot.
            ax: Optional matplotlib axes object. If not provided, a new figure and axes object will be created.
            figsize: Figure size.

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
                a.plot(ax=ax, y_max=100, colour=theoretical_colour, annotate=False)

            # Create a legend
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys())

        # Formatting
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel('m/z')

        return ax

