import matplotlib.pyplot as plt
import numpy as np
from warnings import deprecated
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from . import oligo
from . import enzyme


class Sample:
    def __init__(self, file, well, mz, i, noise_cutoff: float, chip=0, mz_offset=0):
        self.file = file  # source data file
        self.chip = chip    # for AnchorChip plates, chip=0 refers to regular sample spots and chip=1 refers to calibrant spots
        self.well = well
        self.name = self.well   # default naming but this can be changed
        self.noise_cutoff = noise_cutoff
        self.mz_offset = mz_offset  # used to manually adjust the spectra if calibration is off.

        # Optional, experiment dependent attributes
        self.enzyme = enzyme.Enzyme(id=None, name="None")
        self.ntp = None
        # allows flexibility for other experimental designs. Can include anything here, and then it will show up in the excel output
        self.misc_conditions = {
            # 'incubation_time': None,  # examples
            # 'cation': None
        }

        self.mois = []  # MOI = molecule of interest

        # process spectrum
        self.mz_raw = mz    # Keeps a record of the original mz values
        self.mz = mz    # can be updated using the recalc_mz() method
        self.i = i
        self.i_bg_subtracted = []
        self.i_filtered = []

        self.background = min(self.i)
        self.i_background_subtracted()  # populates self.i_bg_subtracted

        self.noise = self.calculate_noise(self.noise_cutoff)
        self.filter_i()   # populates self.i_filtered

        # params for plotting
        self.total_ion = sum(self.i)
        self.max_i = max(self.i)
        self.num_points = len(self.i)

    def recalc_mz(self):
        """Updates self.mz, useful when mz_offset has been changed."""
        self.mz = [point + self.mz_offset for point in self.mz_raw]

    def collect_attributes(self):
        """
        Collects all attributes of the Sample object and returns them as a dictionary.
        Excludes attributes that are only meant to be seen internally.
        """
        include = ['file', 'chip', 'well', 'name', 'noise_cutoff', 'mz_offset', 'background', 'noise']
        exclude = ['misc_conditions', 'mois', 'mz_raw', 'mz', 'i', 'i_bg_subtracted', 'i_filtered', 'total_ion',
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



    def misc_conditions_to_attributes(self):
        """Converts the misc_conditions dictionary to attributes."""
        for k, v in self.misc_conditions.items():
            setattr(self, k, v)

    @deprecated("This was used when processing XML files. The new plain text format does not work with this function.")
    def unpack_spectrum(self):
        """
        Populates self.mz and self.i with values from a spectrum.
        This occurs when the Sample object is instantiated.
        """
        peaks = self.spectrum[0]
        for peak in peaks:
            self.mz.append(float(peak.attrib['mz']))
            self.i.append(float(peak.attrib['i']))

    def moi_intensities(self) -> list[tuple]:
        """
        Returns a list of tuples: (molecule name, intensity value), one for each molecule of interest.
        The intensity value is determined by Oligo.mai_intensity()
        """
        intensities = []
        for moi in self.mois:
            intensities.append(
                (moi.name, moi.mai_intensity(sample=self))
            )
        return intensities

    def calculate_noise(self, cutoff: float) -> int:
        """
        Calculates the noise value of the spectra.
        For our purposes, this is a constant value rather than a function.
        """
        i_vals = self.i_bg_subtracted
        lowest_vals = sorted(i_vals)[:int(len(i_vals)*cutoff)]
        return max(lowest_vals)

    def i_background_subtracted(self) -> None:
        """
        Subtracts the background signal from the spectrum.
        """
        background = self.background
        self.i_bg_subtracted = [i - background for i in self.i]

    def filter_i(self) -> None:
        """
        Keeps only intensity values above the noise of the instrument.
        """
        noise = self.noise

        i_filtered = []
        for i in self.i_bg_subtracted:
            if i - noise > 0:
                i_filtered.append(i)
            else:
                i_filtered.append(0)

        self.i_filtered = i_filtered

    def savitzky_golay(self, i=None, window=15, polyorder=3):
        if i is None:   # allows for smoothing on spectra that have been processed in other ways as well
            i = self.i

        return savgol_filter(i, window_length=window, polyorder=polyorder)

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

    def call_peaks(self, i_list: list) -> list[int]:
        """
        Uses scipy.signal 'find_peaks' to attempt to call peaks on the graph.
        Returns index positions of the peaks based on the list of intensity values provided.
        """
        peak_indices, properties = find_peaks(i_list,
                                              prominence=200,
                                              height=max(i_list) / 10,
                                              threshold=max(i_list) / 10,
                                              distance=100)

        return list(peak_indices)

    def slice_spectrum(self, start: int, end: int, mz: list, i: list) -> tuple:
        """
        Slices the spectrum at indicated start and end values, returns a tuple of lists
        """
        mz_slice = []
        i_slice = []
        for j in range(len(mz)):
            if start <= mz[j] <= end:
                mz_slice.append(mz[j])
                i_slice.append(i[j])
        return mz_slice, i_slice

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

    def plot_moi(self, moi: oligo.Oligo, filtered=True) -> plt.axes:
        """
        Generates a plot for an MOI, zoomed in to visualize the isotope distribution and mass accuracy.
        """
        return self.plot(xlim=moi.iso_dist_range,
                         relative=True,
                         theoretical_dist=True,
                         filtered=filtered,
                         label_mois=False,
                         title=f"{moi.name}, monoisotopic mass = {moi.monoisotopic_mass}, charge = +{moi.charge}")

    def plot(
            self,
            xlim=None,
            relative=False,
            label_peaks=False,
            theoretical_dist=False,
            label_mois=False,
            filtered=False,
            bg_subtracted=False,
            smoothed=False,
            range_of_interest=False,
            title=None,
            colour='#00798c',
            ax=None,
    ) -> plt.axes:

        """
        General purpose plotting function for a Sample's spectrum.
        """

        mz_plot = self.mz
        i_plot = self.i


        if bg_subtracted:   # only works if filtered=False
            i_plot = self.i_bg_subtracted

        if filtered:  # Choose filtered or raw values
            i_plot = self.i_filtered

        if xlim:    # Adjust x axis limits
            mz_plot, i_plot = self.slice_spectrum(start=xlim[0], end=xlim[1], mz=mz_plot, i=i_plot)
        else:
            xlim = (min(mz_plot), max(mz_plot))  # Needed to set plt.xlim later

        if relative:    # Choose relative or absolute intensity
            try:
                i_max = max(i_plot)
                i_plot = [i*100/i_max for i in i_plot]
            except ZeroDivisionError:  # occurs when all i values are 0
                pass


        if ax is None:
            plt.style.use('default')
            fig, ax = plt.subplots(figsize=(8, 3))

        # Base plot
        ax.plot(mz_plot, i_plot, label=f'Exp data', c=colour, linewidth=1, zorder=10)

        if label_peaks:  # Label peaks
            peak_indices = self.call_peaks(i_plot)
            peak_x = [mz_plot[x] for x in peak_indices]
            peak_y = [i_plot[x] for x in peak_indices]
            ax.scatter(peak_x, peak_y, color='red', marker='x', s=2)

            for x, y in zip(peak_x, peak_y):
                ax.annotate(f'{round(x)}', (x, y*1.02), ha='center')

        if theoretical_dist:   # Overlay theoretical isotope distributions
            for moi in self.mois:
                moi.iso_dist_plot(y_max=max(i_plot), standalone=False, colour='#CC6677')

                if range_of_interest:
                    ax.plot(moi.iso_dist_range,  # x1, x2
                            [(max(i_plot)+1)/10] * 2,  # y1, y2 (both values will be max(i_plot)+1/10
                            label='Range of interest', marker='.', markersize=10)

        if label_mois:  # Add indicators for each MOI
            for moi in self.mois:
                xy = moi.monoisotopic_mass, moi.mai_intensity(sample=self) + 0.1 * max(i_plot)
                ax.annotate(text=f"{moi.name}", xy=xy, rotation=0, ha='center')

        # Formatting
        ymax = round(max(i_plot) * 1.2 + 1)  # +1 accounts for possibly all 0 values
        default_title = f'Sample: {self.well}'

        ax.set_ylim(0, ymax)
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_title(default_title) if title is None else plt.title(title)
        ax.set_ylabel(f'Relative intensity (au)') if relative else plt.ylabel(f'Ion count')
        ax.set_xlabel('m/z')

        # Create a legend
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())

        return ax

