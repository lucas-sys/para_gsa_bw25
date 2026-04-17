"""
lca_matrices.py -- Sparse matrix update utilities

This module safely updates technosphere/biosphere sparse matrices under parameter perturbation:
- aggregate_updates_by_cell(): Sum exchange updates mapping to the same matrix cell
- apply_updates_to_matrices(): Apply parameterized update deltas to base matrices

Core formula (preserving non-parameterized contributions):
    new_cell = base_cell - baseline_param_sum + new_param_sum
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, Mapping, Sequence, Tuple


def aggregate_updates_by_cell(
    updates: Sequence[dict],
) -> Dict[Tuple[str, int, int], float]:
    """
    Sum exchange updates that map to the same matrix cell.

    Returns dict keyed by (matrix_name, row_index, col_index).
    """
    grouped = defaultdict(float)
    for item in updates:
        key = (item["matrix"], int(item["row_index"]), int(item["col_index"]))
        grouped[key] += float(item["matrix_value"])
    return dict(grouped)


def apply_updates_to_matrices(
    base_tech,
    base_bio,
    baseline_grouped_updates: Mapping[Tuple[str, int, int], float],
    new_grouped_updates: Mapping[Tuple[str, int, int], float],
):
    """
    Update matrix cells by replacing only the parameterized portion:

        new_cell = base_cell - baseline_param_sum + new_param_sum

    This preserves any non-parameterized contribution already present in the
    stored matrix and correctly handles multiple parameterized exchanges per cell.

    Parameters
    ----------
    base_tech : sparse matrix
        Baseline technosphere matrix (will NOT be mutated).
    base_bio : sparse matrix
        Baseline biosphere matrix (will NOT be mutated).
    baseline_grouped_updates : dict
        {('technosphere_matrix'|'biosphere_matrix', row, col): sum_of_baseline_param_values}
    new_grouped_updates : dict
        Same shape as above but with new parameter values.

    Returns
    -------
    (tech, bio) : tuple of sparse CSR matrices
    """
    tech = base_tech.tolil(copy=True)
    bio = base_bio.tolil(copy=True)

    touched_keys = set(baseline_grouped_updates) | set(new_grouped_updates)

    for matrix, r, c in touched_keys:
        old_param = baseline_grouped_updates.get((matrix, r, c), 0.0)
        new_param = new_grouped_updates.get((matrix, r, c), 0.0)

        if matrix == "technosphere_matrix":
            base_val = float(base_tech[r, c])
            tech[r, c] = base_val - old_param + new_param
        elif matrix == "biosphere_matrix":
            base_val = float(base_bio[r, c])
            bio[r, c] = base_val - old_param + new_param

    return tech.tocsr(), bio.tocsr()


# Backward-compatible alias (used by Monte Carlo module)
apply_grouped_updates_to_matrices = apply_updates_to_matrices
