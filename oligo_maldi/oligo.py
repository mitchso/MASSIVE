import matplotlib.pyplot as plt
from brainpy import isotopic_variants

class Oligo:
    def __init__(self, seq: str, name=None, five_prime_end='OH', ps_bonds=0, error=0, charge=1):
        self.name = name
        self.seq = seq
        self.fiveprime = five_prime_end
        self.composition = self.chem_composition(seq, five_prime_end, ps_bonds)
        self.isotopic_distribution = self.calc_iso_dist(charge=charge, error=error)
        self.monoisotopic_mass = self.calc_monoisotopic_mass()
        self.average_mass = self.calc_avg_mass()
        self.charge = charge

    def __str__(self):
        return f"5' ({self.fiveprime[::-1]})-{self.seq} 3'"

    def composition_str(self):
        return " ".join([k + str(v) for k,v in self.composition.items()])

    def chem_composition(self, seq=None, five_prime_end='OH', ps_bonds=0) -> dict:
        """
        Takes a DNA sequence and returns an elemental composition map.
        """
        if seq is None:
            seq = self.seq

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
        for n in seq:
            for k, v in final_comp.items():
                final_comp[k] += bases[n][k]

        # Add phosphodiester bonds
        po_bonds = len(seq)-1
        for k, v in final_comp.items():
            final_comp[k] += bonds['PO'][k] * po_bonds

        # Substitute phosphorothioate bonds
        for k, v in final_comp.items():
            final_comp[k] += bonds['PS'][k] * ps_bonds

        # Make end adjustments
        for k, v in final_comp.items():
            final_comp[k] += ends[five_prime_end][k]

        return final_comp

    def calc_iso_dist(self, charge=0, error=0) -> tuple:
        """
        Calculates the isotopic distribution of a molecule
        """
        dist = isotopic_variants(self.composition, charge=charge)
        for peak in dist:
            peak.mz += error
        return dist

    def calc_monoisotopic_mass(self) -> float:
        """
        Calculates the monoisotopic mass to 3 decimal places.
        """
        return round(self.isotopic_distribution[0].mz, 3)

    def calc_avg_mass(self):
        """
        Calculates the average mass based on isotopic distribution
        """
        avg_mass = 0
        for peak in self.isotopic_distribution:
            avg_mass += peak.mz * peak.intensity
        return round(avg_mass, 3)

    def iso_dist_plot(self, y_max=None, standalone=True) -> plt.axes:
        x = []
        y = []
        for peak in self.isotopic_distribution:
            x.append(peak.mz)
            y.append(peak.intensity)

        if y_max:   # if y_max is specified, scale everything to match that.
            scale_factor = y_max / max(y)
            y = [n * scale_factor for n in y]

        plt.style.use('default')
        if standalone:  # return a figure
            fig, ax = plt.subplots(figsize=(8,3))
            plt.ylim(0, round(max(y) * 1.5, 2))
            plt.title(f'{self.name}', loc='left')
            plt.ylabel(f'Intensity (au)')
            plt.xlabel('m/z')
            plt.stem(x, y, markerfmt='.', linefmt='#d1495b', label='Theoretical')

            for x,y in zip(x, y):
                plt.annotate(f'{round(x, 2)}', (x, y + 0.025), rotation=90, ha='center')

            return fig, ax

        else:   # add a stem plot to existing figure, return nothing
            plt.stem(x, y, markerfmt='.', linefmt='#d1495b', label='Theoretical')

    def iso_dist_range(self, cumulative_threshold=0.95, left_pad=15, right_pad=15) -> tuple:
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
        while True:
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

    def mai_intensity(self, sample, i_type='filtered') -> float:
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

        start, end = self.iso_dist_range()
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

