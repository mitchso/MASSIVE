import matplotlib.pyplot as plt
from brainpy import isotopic_variants


class Analyte:
    def __init__(self, name, composition, charge=1, custom_mods=None):
        self.name = name
        self.custom_mods = custom_mods
        self.charge = charge    # Assumes this is a +1 charged ion
        self.composition = self._chem_composition(composition, custom_mods)
        self.isotopic_distribution = self._calc_iso_dist()
        self.iso_dist_range = self._calc_iso_dist_range()
        self.monoisotopic_mass = self._calc_monoisotopic_mass()
        self.average_mass = self._calc_avg_mass()

    def __str__(self):
        return self.name

    def composition_str(self):
        return " ".join([k + str(v) for k,v in self.composition.items()])

    def _chem_composition(self, composition:dict, custom_mods:dict=None) -> dict:
        """
        Takes a base chemical composition and combines it with any custom modifications.
        """
        if custom_mods is None:
            custom_mods = {}    # make it a dict so that we can treat it like one

        all_elements = set(composition.keys()) | set(custom_mods.keys())

        final_comp = {e: 0 for e in all_elements}

        for dictionary in [composition, custom_mods]:
            for element, value in dictionary.items():
                final_comp[element] += value

        return final_comp

    def _calc_iso_dist(self, charge=0, error=0) -> tuple:
        """
        Calculates the isotopic distribution of a molecule
        """
        dist = isotopic_variants(self.composition, charge=charge)
        for peak in dist:
            peak.mz += error
        return dist

    def _calc_monoisotopic_mass(self) -> float:
        """
        Calculates the monoisotopic mass to 3 decimal places.
        """
        return round(self.isotopic_distribution[0].mz, 3)

    def _calc_avg_mass(self):
        """
        Calculates the average mass based on isotopic distribution
        """
        avg_mass = 0
        for peak in self.isotopic_distribution:
            avg_mass += peak.mz * peak.intensity
        return round(avg_mass, 3)

    def iso_dist_plot(self, ax=None, y_max=None, annotate=True, label='Theoretical', colour='#d1495b', cumulative_threshold=0.99999) -> plt.axes:
        x = []
        y = []

        mass_start, mass_end = self._calc_iso_dist_range(cumulative_threshold=cumulative_threshold, left_pad=0, right_pad=0)

        for peak in self.isotopic_distribution:
            if mass_start <= peak.mz <= mass_end:
                x.append(peak.mz)
                y.append(peak.intensity)

        if y_max:   # if y_max is specified, scale everything to match that.
            scale_factor = y_max / max(y)
            y = [n * scale_factor for n in y]

        if ax is None:
            plt.style.use('default')
            fig, ax = plt.subplots(figsize=(8, 3))
            plt.ylim(0, round(max(y) * 1.5, 2))
            plt.title(f'{self.name}', loc='left')
            plt.ylabel(f'Intensity (au)')
            plt.xlabel('m/z')

        ax.stem(x, y, markerfmt='.', linefmt=colour, label=label)

        if annotate:
            for x,y in zip(x, y):
                ax.annotate(f'{round(x, 2)}', (x, y + 0.025), rotation=90, ha='center')

        return ax


    def _calc_iso_dist_range(self, cumulative_threshold=0.95, left_pad=10, right_pad=10) -> tuple:
        """
        Calculates the mass range where isotopes of the same molecule may be observed.

        Since isotope distributions can have long tails, cumulative_threshold cuts the range to where
        95% of the signal is expected.

        Padding allows some extra room (daltons) on either side of the calculated distribution.

        Returns a tuple of (start, end) range.
        """

        cumulative_signal = 0
        cumulative_mz_vals = []
        i = 0
        while i < len(self.isotopic_distribution):
            mz = self.isotopic_distribution[i].mz
            intensity = self.isotopic_distribution[i].intensity
            if cumulative_signal + intensity <= cumulative_threshold:
                cumulative_mz_vals.append(mz)
                cumulative_signal += intensity
                i += 1
            else:
                break

        cumulative_range = (min(cumulative_mz_vals)-left_pad, max(cumulative_mz_vals)+right_pad)
        return cumulative_range

    def _mai_intensity(self, sample, i_type='filtered') -> float:
        """
        Returns the Most Abundant Isotope (MAI) intensity of this molecule in a given sample.
        """
        if i_type == 'raw':
            i_vals = sample.i
        elif i_type == 'bg_sub':
            i_vals = sample.i_bg_subtracted
        elif i_type == 'filtered':
            i_vals = sample.i_filtered
        else:
            raise ValueError("i_type must be raw, bg_sub, or filtered.")

        start, end = self.iso_dist_range
        mz = []
        i = []
        for j in range(len(sample.mz)):
            if start <= sample.mz[j] <= end:
                mz.append(sample.mz[j])
                i.append(i_vals[j])

        try:
            mai = max(i)
        except ValueError:
            mai = 0

        return mai


class Oligo(Analyte):
    def __init__(self, name, seq: str, type='DNA', five_prime_end='OH', charge=1, custom_mods=None):
        # set up oligo specific attributes
        self.seq = seq
        self.type = type
        self.five_prime_end = five_prime_end

        # use these to calculate elemental composition
        composition = self._oligo_composition()

        # pass this elemental composition to the base Analyte class
        super().__init__(name, composition, charge, custom_mods)

    def _oligo_composition(self) -> dict:
        """
        Takes an oligo sequence and returns an elemental composition map.
        """

        bases = {'A': {'C': 10, 'H': 13, 'N': 5, 'O': 3, 'P': 0, 'S': 0},
                 'T': {'C': 10, 'H': 14, 'N': 2, 'O': 5, 'P': 0, 'S': 0},
                 'C': {'C': 9, 'H': 13, 'N': 3, 'O': 4, 'P': 0, 'S': 0},
                 'G': {'C': 10, 'H': 13, 'N': 5, 'O': 4, 'P': 0, 'S': 0}
                 }

        ends = {'OH': {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0, 'S': 0},
                'P': {'C': 0, 'H': 1, 'N': 0, 'O': 3, 'P': 1, 'S': 0},
                'PP': {'C': 0, 'H': 2, 'N': 0, 'O': 6, 'P': 2, 'S': 0},
                'PPP': {'C': 0, 'H': 3, 'N': 0, 'O': 9, 'P': 3, 'S': 0}
                }

        bonds = {'PS': {'C': 0, 'H': 0, 'N': 0, 'O': -1, 'P': 0, 'S': 1},  # Difference relative to PO bond
                 'PO': {'C': 0, 'H': -1, 'N': 0, 'O': 2, 'P': 1, 'S': 0}
                 }

        final_comp = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0, 'S': 0}

        # Add up base compositions
        for n in self.seq:
            for k, v in final_comp.items():
                final_comp[k] += bases[n][k]

        # type adjustment
        if self.type == 'DNA':
            pass
        elif self.type == 'RNA':
            final_comp['O'] += len(self.seq)
        else:
            raise NotImplementedError("Unexpected Oligo.type")

        # Add phosphodiester bonds
        po_bonds = len(self.seq)-1
        for k, v in final_comp.items():
            final_comp[k] += bonds['PO'][k] * po_bonds

        # Make end adjustments
        for k, v in final_comp.items():
            final_comp[k] += ends[self.five_prime_end][k]

        return final_comp