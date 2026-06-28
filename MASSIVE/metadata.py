"""
Meant to contain metadata classes that can be added to [`Samples`][MASSIVE.sample.Sample] if desired. `Enzyme` is an example. Metadata can enrich visualizations (see [`Plate`][MASSIVE.collections.Plate].show_sample_positions() as an example).
"""

class Enzyme(object):
    """
    Used to standardize naming and retrieval of enzyme variant attributes.
    """
    def __init__(self, id:int, name:str, mutation:str="WT", organism:str=None):
        """

        Args:
            id: Unique, incremental identifier for each variant.
            name: Enzyme name, e.g. 'T7_RNAP'
            mutation: Amino acid substitution relative to WT, like "T47A" or "WT"
            organism: Organism name, e.g. 'T7_bacteriophage'
        """
        self.id = id
        self.name = name
        self.mutation = mutation
        self.organism = organism

    def __str__(self):
        return f"{self.id}, {self.mutation}"

    def get_aa_position(self) -> int|None:
        """
        returns integer value of amino acid position from `Enzyme.mutation`.

        Examples:<br>
            "D83A" return int(83)<br>
            "S340P" returns int(340)<br>
            "WT" returns None<br>
            "blah" returns None<br>
        """
        try:
            pos = int(self.mutation[1:-1])
            return pos
        except ValueError:  # happens if there is text entered in self.mutation
            return None
        except TypeError:   # happens if self.mutation is None
            return None

