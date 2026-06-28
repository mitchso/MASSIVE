import matplotlib.pyplot as plt
import matplotlib.axes
from brainpy import isotopic_variants


class Analyte:
    """
    Base class representing a molecule with a defined elemental composition.

    The molecule's monoisotopic mass, average mass, and isotopic distribution are automatically calculated (see `Attributes`).

    Analytes can be assigned to [`Sample`][MASSIVE.sample.Sample] objects ([`Sample`][MASSIVE.sample.Sample].analytes), which then enables the Analyte to be automatically detected and quantified.

    The base Analyte class is typically used for small molecules, while the subclasses [`Oligo`][MASSIVE.analytes.Oligo] and [`Peptide`][MASSIVE.analytes.Peptide] are helpful abstractions for defining larger molecules in terms of their sequences (i.e. 'ACTGTA') and their modifications (i.e. methylation, phosphorylation)
    """
    def __init__(self, name:str, composition:dict, charge:int = +1, mods:None|list|dict|str=None):
        """
        Args:
            name: Human-readable name for the analyte.
            composition: Elemental composition as a dict e.g. {'C': 10, 'H': 13, 'N': 5}.
            charge: Ion charge state. Defaults to a singly charged positive ion.
            mods: Optional modifications.

        Note:
        `mods` accepts several input formats:

        - **None** — no modification applied.
        - **str** — a single named modification e.g. `'methyl'`.
        - **dict** — elemental changes e.g. `{'C': 1, 'O':1, 'H': -2}`.
        - **list** — multiple modifications, which can be either strings or dicts e.g. `['methyl', `{'C': 1, 'O':1, 'H': -2}`]`.

        Modifications given as strings are resolved depending on the subclass.
        `Oligo` and `Peptide` have their own sets of modifications, while `Analyte` is too broad to reasonably define a set of modifications, so none are accepted.
        See subclasses such as [`Oligo`][MASSIVE.analytes.Oligo] for more details.

        Each item in a list of modifications is resolved additively, according to the logic above.

        Attributes:
            name (str): Human-readable name for the analyte.
            mods (None|list|dict|str): A record of any modifications applied to the molecule, in the format they were given as input.
            composition (dict): Final elemental composition as a dict, including any modifications.
            charge (int): Ion charge state.
            monoisotopic_mass (float): Mass of the most abundant isotopologue in Daltons,
                rounded to 3 decimal places.
            average_mass (float): Intensity-weighted average mass across the isotopic
                distribution, rounded to 3 decimal places.
            isotopic_distribution (list): Full isotopic distribution as a list of peaks,
                each with `.mz` and `.intensity` attributes.
            iso_dist_range(tuple): Defaults to a tuple of (start, end) mass range covering 95% of the
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
        """Returns a string representation of the elemental composition, i.e. `C146 H182 N67 O85 P15`"""
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

    def plot(self, ax:matplotlib.axes.Axes=None, y_max: int | float=None, annotate:bool=True, label:str= 'Theoretical', colour:str= '#d1495b', cumulative_threshold:float=0.99999) -> matplotlib.axes.Axes:
        """
        Generates a stem plot of the isotopic distribution of the molecule. Useful for visualizing the isotopic distribution.

        Args:
            ax: Optionally, provide an existing axes object to plot on. If no axes object is provided, a new figure and axes are created.
            y_max: Optionally, provide a maximum y-axis value to scale the plot to.
            annotate: If True, annotate each peak with its exact m/z value (2 decimal places).
            label: Label for the plot legend.
            colour: Colour for the stem plot.
            cumulative_threshold: Cutoff for the isotopic distribution. Since distributions can have long tails, this keeps the plot more manageable.

        Returns:
            An axes object containing the stem plot.

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

    def _peak_intensity(self, sample, i_type='filtered') -> float:
        """
        Returns the peak intensity of this molecule in a given sample.
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
    """
    Represents an oligonucleotide (DNA or RNA).

    `Oligo` is useful for defining an analyte in terms of its nucleotide sequence, as well as any modifications (e.g. methylation, phosphorylation).

    Common bases and modifications are defined in `Attributes`. These dictionaries can be expanded as you wish to include additional bases and modifications.
    For `Oligo`, valid names include `'methyl'`, `'PS'`, `'3P'`, etc.

    """
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

    def __init__(self, name:str, seq: str, type:str='DNA', charge:int = +1, mods:None|list|dict|str=None):
        """
        Args:
            name: Human-readable name for the Oligo.
            seq: Sequence of nucleotides, e.g. 'ACTGTA'.
            type: `DNA` or `RNA`.
            charge: Ion charge state. Defaults to singly charged positive ion.
            mods: See [`Analyte`][MASSIVE.analytes.Analyte] for information about how to define modifications.


        Attributes:
            seq (str): Sequence of nucleotides, e.g. 'ACTGTA'.
            type (str): `DNA` or `RNA`.
            composition (dict): Elemental composition of the oligo, including any modifications.
            mods (dict): Elemental composition of modifications applied to the oligo.
            BASES (dict): Dictionary of nucleotides and their corresponding elemental compositions. Includes `'A', 'T', 'C', 'G', 'U', 'I'`. For an exact list and molecular compositions, see `Oligo.BASES` in the source code.
            KNOWN_MODIFICATIONS (dict): Dictionary of modifications and their corresponding elemental composition changes. For an exact list and molecular compositions, see `Oligo.KNOWN_MODIFICATIONS` in the source code.

        | Example | Description | Net elemental change |
        |------|-------------|------------------|
        | `methyl` | Methylation | C: +1, H: +2 (gain a carbon and 3 hydrogens but lose 1 hydrogen) |
        | `PS` | Phosphorothioation | O: -1, S: +1 (replace oxygen with sulphur) |
        """

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
    """
    <b>STILL UNDER CONSTRUCTION.</b>

    Represents a peptide.

    `Peptide` is useful for defining an analyte in terms of its amino acid sequence, as well as any modifications (e.g. methylation, phosphorylation).

    Common amino acids and modifications are defined in `Attributes`. These dictionaries can be expanded as you wish to include additional bases and modifications.

    """

    AMINO_ACIDS = {}

    BONDS = {}

    KNOWN_MODIFICATIONS = {}

    def __init__(self, name: str, seq: str, charge:int = +1, mods:None|list|dict|str=None):
        """
        Args:
            name: Human-readable name for the Peptide.
            seq: Sequence of amino acids, e.g. 'METKAV'.
            charge: Ion charge state. Defaults to a singly charged positive ion.
            mods: See [`Analyte`][MASSIVE.analytes.Analyte] for information about how to define modifications.

        Attributes:
            seq (str): Sequence of amino acids, e.g. 'METKAV'.
            composition (dict): Elemental composition of the peptide, including any modifications.
            mods (dict): Elemental composition of modifications applied to the peptide.
            AMINO_ACIDS (dict): Dictionary of amino acids and their corresponding elemental compositions. For an exact list and molecular compositions, see `Oligo.AMINO_ACIDS` in the source code.
            KNOWN_MODIFICATIONS (dict): Dictionary of modifications and their corresponding elemental composition changes. For an exact list and molecular compositions, see `Peptide.KNOWN_MODIFICATIONS` in the source code.

        """

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