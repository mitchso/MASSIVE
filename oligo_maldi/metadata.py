"""
Meant to contain metadata classes that can be added to Samples if desired.
"""

class Enzyme(object):
    """
    Used to standardize naming and retrieval of enzyme variant attributes
    """
    def __init__(self, id, name, mutation="WT", organism=None):
        self.id = id    # unique, incremental identifier for each variant
        self.name = name    # ie T7_RNAP
        self.mutation = mutation    # ie "T47A" or "WT"
        self.organism = organism # derived from?

    def __str__(self):
        return f"{self.id}, {self.mutation}"

    def get_aa_position(self):
        """
        returns integer value of amino acid position from single character labels.
        EXAMPLE:    "D83A" return int(83)
                    "S340P" returns int(340)
                    "blah" returns 0
        """
        try:
            pos = int(self.mutation[1:-1])
            return pos
        except ValueError:  # happens if there is text entered in self.mutation
            return 0
        except TypeError:   # happens if self.mutation is None
            return None

