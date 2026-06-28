## Enzyme

**File:** `enzyme.py`

Represents an enzyme variant. Used to attach metadata to `Sample` objects and to compute activity scores.

### Constructor

```python
Enzyme(id, name, mutation='WT', round=None, organism=None)
```

| Parameter | Description |
|---|---|
| `id` | Unique integer identifier |
| `name` | Enzyme name (e.g. `'DdiTLP4'`) |
| `mutation` | Mutation label (e.g. `'T47A'`) or `'WT'` |
| `round` | Directed evolution round number |
| `organism` | Source organism |

---

### Methods

#### `get_aa_position() → int`
Parses `self.mutation` and returns the integer amino acid position. For example, `'D83A'` returns `83`. Returns `0` if the mutation string is not in the expected format, or `None` if `mutation` is `None`.

#### `default_activity_score(unreacted, correct, a1, a2, b1, b2) → float` *(static)*
The primary activity scoring function. Computes a composite score from the fraction of unreacted substrate and the fraction of correctly processed product.

- Returns `0` if `unreacted == 1` (no reaction occurred).
- Otherwise: `score = (a1 * reacted^b1 + a2 * correctness^b2) / (a1 + a2)`

Additional scoring variants are available if needed:

| Method | Notes |
|---|---|
| `activity_score_exponential` | Uses `exp(reacted) + correctness^b2` |
| `activity_score_multiplied` | Uses `reacted^b1 * correctness^b2` |
| `activity_score_r_plus_rc` | Uses a linear + quadratic composite |

All scoring methods are `@staticmethod` and share the same signature: `(unreacted, correct, a1, a2, b1, b2)`.

---