"""
lca_contribution.py -- Contribution analysis, OAT sensitivity, and analytical variance

This module provides:
- contribution_analysis(): Recursive contribution analysis using bw2analyzer
- oat_sensitivity(): One-at-a-time finite-difference sensitivity analysis
- analytical_variance_from_parameter(): Compute analytical variance from uncertainty metadata
- parameter_metadata(): Export parameter metadata to a DataFrame
- combine_oat_and_analytical_variance(): Merge OAT coefficients with analytical variance

Internal utilities:
- _reset_and_recalculate_mlca(): Clear cached attributes and re-run LCI/LCIA
- _recompute_scores_with_param_context(): Recompute scores with a perturbed parameter context
"""

from __future__ import annotations

from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

import bw2analyzer as bwa
from bw2calc.multi_lca import MultiLCA

from utils.lca_config import LCASetup
from utils.lca_project_setup import (
    initialize_multilca,
    flatten_multilca_scores,
)
from utils.lca_parameters import (
    get_activity_parameters,
    build_parameter_formula_map,
    build_exchange_formula_map,
    resolve_parameter_context,
    evaluate_exchange_formulas,
)
from utils.lca_matrices import aggregate_updates_by_cell, apply_updates_to_matrices


# -----------------------------------------------------------------------------
# Contribution analysis
# -----------------------------------------------------------------------------

def contribution_analysis(
    functional_units: Sequence,
    methods: Sequence[tuple],
    max_level: int = 1,
    cutoff: float = 0.001,
) -> pd.DataFrame:
    """Recursive contribution analysis for each (FU, method) pair."""
    rows = []
    for method in methods:
        for fu in functional_units:
            obj = bwa.utils.recursive_calculation_to_object(
                activity=fu,
                lcia_method=method,
                amount=1,
                max_level=max_level,
                cutoff=cutoff,
            )
            df = pd.DataFrame(obj)
            root_score = df.loc[df["label"] == "root", "score"].iloc[0]
            sub = df[df["parent"] == "root"].copy()
            row = sub.set_index("name")["score"].to_dict()
            row.update(
                {
                    "Scenario": fu.get("name", fu["code"]),
                    "Method": method[2] if len(method) > 2 else str(method),
                    "Total score": root_score,
                }
            )
            rows.append(row)
    return pd.DataFrame(rows).fillna(0)


# -----------------------------------------------------------------------------
# Matrix recalculation helpers
# -----------------------------------------------------------------------------

def _reset_and_recalculate_mlca(mlca: MultiLCA, tech_matrix, bio_matrix) -> None:
    """Replace matrices and force a clean recalculation."""
    mlca.technosphere_matrix = tech_matrix.copy().tocsr()
    mlca.biosphere_matrix = bio_matrix.copy().tocsr()

    for attr in (
        "solver",
        "supply_array",
        "inventory",
        "characterized_inventory",
        "normalized_inventory",
        "weighted_inventory",
    ):
        if hasattr(mlca, attr):
            try:
                delattr(mlca, attr)
            except Exception:
                pass

    mlca.lci()
    mlca.lcia()


def _recompute_scores_with_param_context(
    setup: LCASetup,
    group_name: str,
    param_context: Mapping[str, float],
    exchange_formula_map: Optional[Sequence[dict]] = None,
    baseline_param_context: Optional[Mapping[str, float]] = None,
) -> Dict[Tuple[str, tuple], float]:
    """Recompute LCA scores using a perturbed parameter context."""
    exchange_formula_map = list(exchange_formula_map or build_exchange_formula_map(group_name, setup.mlca))

    if baseline_param_context is None:
        baseline_param_context = resolve_parameter_context(group_name, {})

    baseline_updates = evaluate_exchange_formulas(exchange_formula_map, baseline_param_context)
    new_updates = evaluate_exchange_formulas(exchange_formula_map, param_context)

    baseline_grouped = aggregate_updates_by_cell(baseline_updates)
    new_grouped = aggregate_updates_by_cell(new_updates)

    tech, bio = apply_updates_to_matrices(
        setup.base_tech,
        setup.base_bio,
        baseline_grouped,
        new_grouped,
    )

    _reset_and_recalculate_mlca(setup.mlca, tech, bio)

    return flatten_multilca_scores(setup.mlca.scores)


# -----------------------------------------------------------------------------
# OAT sensitivity
# -----------------------------------------------------------------------------

def oat_sensitivity(
    functional_units: Sequence,
    methods: Sequence[tuple],
    group_name: str,
    rel_step: float = 0.1,
    absolute_step: Optional[float] = None,
) -> pd.DataFrame:
    """One-at-a-time finite-difference sensitivities for group parameters."""
    setup = initialize_multilca(functional_units, methods)
    params = get_activity_parameters(group_name)
    param_names = [p.name for p in params]
    parameter_formula_map = build_parameter_formula_map(group_name)
    exchange_formula_map = build_exchange_formula_map(group_name, setup.mlca)

    base_context = resolve_parameter_context(group_name, {}, parameter_formula_map)
    base_scores = flatten_multilca_scores(setup.mlca.scores)

    rows = []
    for name in param_names:
        x0 = base_context[name]
        dx = absolute_step if absolute_step is not None else max(abs(x0) * rel_step, 1e-12)
        up_context = dict(base_context)
        up_context[name] = x0 + dx
        up_context = resolve_parameter_context(group_name, up_context, parameter_formula_map)
        up_scores = _recompute_scores_with_param_context(
            setup, group_name, up_context, exchange_formula_map
        )

        for key, y0 in base_scores.items():
            y1 = up_scores[key]
            rows.append(
                {
                    "Parameter": name,
                    "Functional unit": key[0],
                    "Method": key[1][2] if len(key[1]) > 2 else str(key[1]),
                    "Base parameter value": x0,
                    "Delta x": dx,
                    "Base score": y0,
                    "Perturbed score": y1,
                    "Delta y": y1 - y0,
                    "Sensitivity ratio": ((y1 - y0) / dx) * abs(x0 / y0) if x0 != 0 and y0 != 0 else np.nan,
                    "Sensitivity coefficient": (y1 - y0) / dx,
                    "Sensitivity coefficient squared": ((y1 - y0) / dx) ** 2,
                }
            )
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Analytical variance from uncertainty metadata
# -----------------------------------------------------------------------------

def analytical_variance_from_parameter(par) -> float:
    """
    Compute the analytical variance for a single Brightway ActivityParameter.

    Maps Brightway uncertainty type codes to scipy.stats distributions.
    """
    p = par.dict
    code = p.get("uncertainty type", 0)
    loc = p.get("loc", 0)
    scale = p.get("scale", 1)
    shape = p.get("shape", 1)
    p_min = p.get("minimum", 0)
    p_max = p.get("maximum", 1)
    shape2 = p.get("shape2", 1)

    if code in (0, 1):
        return 0.0
    if code == 2:
        return float(stats.lognorm.var(s=scale, loc=0, scale=np.exp(loc)))
    if code == 3:
        return float(stats.norm.var(loc=loc, scale=scale))
    if code == 4:
        return float(stats.uniform.var(loc=p_min, scale=p_max - p_min))
    if code == 5:
        width = p_max - p_min
        if width == 0:
            return 0.0
        c = (loc - p_min) / width
        return float(stats.triang.var(c=c, loc=p_min, scale=width))
    if code == 6:
        return float(stats.bernoulli.var(p=loc))
    if code == 7:
        return float(stats.randint.var(low=p_min, high=p_max + 1))
    if code == 8:
        return float(stats.weibull_min.var(c=shape, loc=loc, scale=scale))
    if code == 9:
        return float(stats.gamma.var(a=shape, loc=loc, scale=scale))
    if code == 10:
        return float(stats.beta.var(a=shape, b=shape2, loc=loc, scale=scale))
    if code == 11:
        return float(stats.genextreme.var(c=shape, loc=loc, scale=scale))
    if code == 12:
        return float(stats.t.var(df=shape, loc=loc, scale=scale))
    return np.nan


def parameter_metadata(group_name: Optional[str] = None) -> pd.DataFrame:
    """Export parameter metadata including analytical variance to a DataFrame."""
    rows = []
    for par in get_activity_parameters(group_name):
        p = par.dict
        rows.append(
            {
                "Parameter": par.name,
                "Parameter amount": par.amount,
                "Group": getattr(par, "group", None),
                "uncertainty type": p.get("uncertainty type", 0),
                "loc": p.get("loc", 0),
                "scale": p.get("scale", 1),
                "shape": p.get("shape", 1),
                "shape2": p.get("shape2", 1),
                "minimum": p.get("minimum", 0),
                "maximum": p.get("maximum", 1),
                "Analytical variance": analytical_variance_from_parameter(par),
            }
        )
    return pd.DataFrame(rows)


def combine_oat_and_analytical_variance(
    sensitivity_df: pd.DataFrame,
    param_meta_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merge OAT sensitivity results with parameter metadata and compute variance contributions."""
    out = sensitivity_df.merge(param_meta_df, on="Parameter", how="left")
    if "Sensitivity coefficient squared" not in out.columns:
        out["Sensitivity coefficient squared"] = out["Sensitivity coefficient"] ** 2
    out["Analytical variance contribution"] = (
        out["Sensitivity coefficient squared"] * out["Analytical variance"]
    )
    out["Analytical variance contribution"] = out[
        "Analytical variance contribution"
    ].replace([np.inf, -np.inf], np.nan)
    return out
