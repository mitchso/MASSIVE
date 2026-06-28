import matplotlib.pyplot as plt
from brainpy import isotopic_variants


class Analyte:
    """
    Base class representing any molecule with a specific elemental composition.

    Attributes such as monoisotopic mass, average mass, and isotopic distribution are automatically calculated.

    The base Analyte class is typically useful for small molecules, while the subclasses Oligo and Peptide are helpful for defining larger molecules in terms of their sequences (i.e. 'ACTGTA') and their modifications (i.e. methylation, phosphorylation)
    """
    def __init__(self, name:str, composition:dict, charge:int=1, mods:None|list|dict|str=None):
        """
        Args:
            name: Human-readable name for the analyte.
            composition: Elemental composition as a dict e.g. {'C': 10, 'H': 13, 'N': 5}.
            charge: Ion charge state. Defaults to 1.
            mods: Optional modifications.

        Note:
        `mods` accepts several input formats:

        - **None** — no modification applied.
        - **str** — a single named modification e.g. `'methyl'`.
        - **dict** — elemental changes e.g. `{'C': 1, 'O':1, 'H': -2}`.
        - **list** — multiple named modifications, which can be either strings or dicts e.g. `['methyl', `{'C': 1, 'O':1, 'H': -2}`]`.

        Modifications given as strings are resolved via `KNOWN_MODIFICATIONS` on the subclass.
        For `Oligo`, valid names include `'methyl'`, `'PS'`, `'3P'`, `'5PPP'`, etc.
        See `Oligo.KNOWN_MODIFICATIONS` for the full list.

        Attributes:
            name: Human-readable name for the analyte.
            mods: A record of any modifications applied to the molecule, in the format they were given as input.
            composition: Final elemental composition as a dict, including any modifications.
            charge: Ion charge state.
            monoisotopic_mass: Mass of the most abundant isotopologue in Daltons,
                rounded to 3 decimal places.
            average_mass: Intensity-weighted average mass across the isotopic
                distribution, rounded to 3 decimal places.
            isotopic_distribution: Full isotopic distribution as a list of peaks,
                each with `.mz` and `.intensity` attributes.
            iso_dist_range: Defaults to a tuple of (start, end) mass range covering 95% of the
                isotopic signal, with 10 Da padding on each side.

        """

        self.name = name
        self.mods = mods
        self.charge = charge    # Assumes this is a +1 charged ion
        self.composition = self._chem_composition(composition, mods)
        self.isotopic_distribution = self._calc_iso_dist()
        self.iso_dist_range = self._calc_iso_dist_range()
        self.monoisotopic_mass = self._calc_monoisotopic_mass()
        self.average_mass = self._calc_avg_mass()

    def __str__(self):
        return self.name

    def composition_str(self):
        """Returns a string representation of the elemental composition, i.e. C146 H182 N67 O85 P15"""
        return " ".join([k + str(v) for k,v in self.composition.items()])

    def _resolve_modifications(self, mods, known_modifications: dict) -> dict | None:
        """
        Resolves modification names into elemental compositions.

        mods can be:
          - None
          - a string: 'methyl'
          - a list of strings: ['methyl', '3phos']
          - a dict (already resolved): {'C': 1, 'H': 2, ...}
        """
        # Already a resolved elemental composition — pass through
        if isinstance(mods, dict) or mods is None:
            return mods

        # make sure mods is a list
        if not isinstance(mods, list):
            mods = [mods]

        resolved = {}
        for mod in mods:
            if isinstance(mod, dict):
                for element, count in mod.items():
                    resolved[element] = resolved.get(element, 0) + count

            elif isinstance(mod, str):
                if mod in known_modifications:
                    for element, count in known_modifications[mod].items():
                        resolved[element] = resolved.get(element, 0) + count
                else:
                    raise ValueError(
                        f"Unknown modification '{mod}'. "
                        f"Known modifications are listed in the KNOWN_MODIFICATIONS class attribute."
                    )

        return resolved


    def _chem_composition(self, composition:dict, mods: dict | None=None) -> dict:
        """
        Takes a base chemical composition and combines it with any custom modifications.
        """
        if mods is None:
            mods = {}    # make it a dict so that we can treat it like one

        # get all elements used between the base oligo and the modifications
        all_elements = set(composition.keys()) | set(mods.keys())

        final_comp = {e: 0 for e in all_elements}

        for dictionary in [composition, mods]:
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
        """
        Test

        """
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

        ax.stem(x, y, markerfmt='.', basefmt=colour, linefmt=colour, label=label)

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
    BASES = {
        'A': {'C': 10, 'H': 13, 'N': 5, 'O': 4, 'P': 0},
        'T': {'C': 10, 'H': 14, 'N': 2, 'O': 6, 'P': 0},
        'C': {'C': 9,  'H': 13, 'N': 3, 'O': 5, 'P': 0},
        'G': {'C': 10, 'H': 13, 'N': 5, 'O': 5, 'P': 0},
        'U': {'C': 9,  'H': 12, 'N': 2, 'O': 6, 'P': 0},
        'I': {'C': 10, 'H': 12, 'N': 4, 'O': 5, 'P': 0},
    }

    BONDS = {
        'PO': {'C': 0, 'H': -1, 'N': 0, 'O': 2, 'P': 1, 'S': 0}
    }

    KNOWN_MODIFICATIONS = {
        'methyl':           {'C': 1, 'H': 2, 'N': 0, 'O': 0, 'P': 0, 'S': 0},
        'hydroxymethyl':    {'C': 1, 'H': 3, 'N': 0, 'O': 1, 'P': 0, 'S': 0},
        'formyl':           {'C': 1, 'H': 1, 'N': 0, 'O': 1, 'P': 0, 'S': 0},
        'carboxy':          {'C': 1, 'H': 0, 'N': 0, 'O': 2, 'P': 0, 'S': 0},
        'PS':               {'C': 0, 'H': 0, 'N': 0, 'O': -1, 'P': 0, 'S': 1},
        'adenylation':      {'C': 10, 'H': 14, 'N': 5, 'O': 5, 'P': 1, 'S': 0},
        '3P':               {'C': 0, 'H': 1, 'N': 0, 'O': 3, 'P': 1, 'S': 0},
        '5P':               {'C': 0, 'H': 1, 'N': 0, 'O': 3, 'P': 1, 'S': 0},
        '5PP':              {'C': 0, 'H': 2, 'N': 0, 'O': 6, 'P': 2, 'S': 0},
        '5PPP':             {'C': 0, 'H': 3, 'N': 0, 'O': 9, 'P': 3, 'S': 0},
        '5AmMC12':          {'C': 12, 'H': 26, 'N': 1, 'O': 3, 'P': 1, 'S': 0}   # IDT
    }

    def __init__(self, name, seq: str, type='DNA', charge=1, mods=None):
        # set up oligo specific attributes
        self.seq = seq
        self.type = type

        # use these to calculate elemental composition
        composition = self._oligo_composition()
        mods = self._resolve_modifications(mods, self.KNOWN_MODIFICATIONS)


        # pass this elemental composition to the base Analyte class
        super().__init__(name, composition, charge, mods)

    def _oligo_composition(self) -> dict:
        """
        Takes an oligo sequence and returns an elemental composition map.
        """
        if len(self.seq) < 1:
            raise ValueError("Oligo sequence must be at least one base long.")

        composition = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0}

        # Add up base compositions
        for N in self.seq:
            for ele in composition.keys():
                composition[ele] += self.BASES[N][ele]

        # type adjustment
        if self.type == 'DNA':
            composition['O'] -= len(self.seq)   # remove 1 oxygen per base
        elif self.type == 'RNA':
            pass
        else:
            raise NotImplementedError("Unexpected Oligo.type")

        # Add phosphodiester bonds
        po_bonds = len(self.seq)-1
        for ele in composition.keys():
            composition[ele] += self.BONDS['PO'][ele] * po_bonds

        return composition


# TODO: build this class
class Peptide(Analyte):
    AMINO_ACIDS = {}

    BONDS = {}

    KNOWN_MODIFICATIONS = {}

    def __init__(self, name, seq: str, charge=1, mods=None):
        self.seq = seq

        # use these to calculate elemental composition
        composition = self._peptide_composition()
        mods = self._resolve_modifications(mods, self.KNOWN_MODIFICATIONS)

        # pass this elemental composition to the base Analyte class
        super().__init__(name, composition, charge, mods)

    def _peptide_composition(self) -> dict:
        """
        Takes an oligo sequence and returns an elemental composition map.
        """
        composition = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0, 'S': 0}

        return composition  # stub