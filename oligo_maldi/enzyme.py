import math

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

    @staticmethod
    def default_activity_score(unreacted: float, correct: float, a1=1, a2=1, b1=2, b2=2):
        """Converts measurements of unreacted, correct and incorrect products into an aggregate score."""
        if unreacted == 1:
            score = 0

        else:
            reacted = 1 - unreacted
            correctness = correct / reacted

            score = (a1 * reacted**b1 + a2 * correctness**b2) / (a1 + a2)

        return score

    @staticmethod
    def activity_score_exponential(unreacted: float, correct: float, a1=1, a2=1, b1=2, b2=2):
        """Converts measurements of unreacted, correct and incorrect products into an aggregate score."""
        # TODO: verify implementation
        reacted = 1 - unreacted

        if unreacted == 1:
            correctness = 1
        else:
            correctness = correct / reacted

        score = math.exp(reacted) + correctness**b2

        return score

    @staticmethod
    def activity_score_multiplied(unreacted: float, correct: float, a1=1, a2=1, b1=2, b2=2):
        """Converts measurements of unreacted, correct and incorrect products into an aggregate score."""
        # TODO: verify implementation
        if unreacted == 1:
            score = 0

        else:
            reacted = 1 - unreacted
            correctness = correct / reacted

            score = (reacted**b1 * correctness**b2)

        return score

    @staticmethod
    def activity_score_r_plus_rc(unreacted: float, correct: float, a1=1, a2=1, b1=1, b2=2):
        """Converts measurements of unreacted, correct and incorrect products into an aggregate score."""
        if unreacted == 1:
            score = 0

        else:
            reacted = 1 - unreacted
            correctness = correct / reacted

            score = (a1 * reacted ** b1 + a2 * reacted * correctness ** b2) / (a1 + a2)

        return score

    @staticmethod
    def activity_score_c_only(unreacted: float, correct: float, a1=1, a2=1, b1=1, b2=2):
        """Converts measurements of unreacted, correct and incorrect products into an aggregate score."""
        if unreacted == 1:
            score = 0

        else:
            reacted = 1 - unreacted
            correctness = correct / reacted

            score = reacted * correctness

        return score

