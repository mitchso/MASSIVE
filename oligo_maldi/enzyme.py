import pandas as pd

class Enzyme(object):
    """
    Used to standardize naming and retrieval of enzyme variant attributes
    """
    def __init__(self, id, name, mutation="WT", round=None, organism=None):
        self.id = id    # unique, incremental identifier for each variant
        self.name = name    # ie DdiTLP4
        self.mutation = mutation    # ie "T47A" or "WT"
        self.round = round  # for directed evolution experiments
        self.organism = organism # derived from?

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
