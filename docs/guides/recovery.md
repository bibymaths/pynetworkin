# False-Negative Recovery

## Background

The original NetworKIN used network context (STRING) only as a **false-positive filter** applied on top of motif predictions. Any kinase–substrate pair the motif step missed was permanently excluded.

pynetworkin introduces a **false-negative recovery** step to rescue such pairs.

## Algorithm

A pair `(kinase, substrate)` is **recovered** if:

1. It was **not** scored by the motif step (absent from `motif_scores`).
2. Its context score ≥ `CONTEXT_RECOVERY_THRESHOLD` (default: 0.6).

The context score is computed as:

```
context_score = 1 / (1 + shortest_path_distance)
```

where `shortest_path_distance` is from the Floyd–Warshall all-pairs matrix on the STRING network.

A threshold of 0.6 corresponds to an effective path length ≤ 0.67.

## Output Columns

Recovered predictions include:

| Column | Value |
|--------|-------|
| `recovered` | `True` |
| `recovery_method` | `"context_proximity"` |
| `motif_score` | `-1.0` (sentinel, not scored) |

## Tuning

Adjust the threshold via:

```python
from recovery import CONTEXT_RECOVERY_THRESHOLD
# or set programmatically before calling recover_false_negatives()
```

Higher threshold → fewer but more confident recoveries.
Lower threshold → more recoveries, higher false-positive rate.
