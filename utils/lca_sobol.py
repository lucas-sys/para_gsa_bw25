"""
lca_sobol.py -- Sobol global sensitivity analysis

This module provides:
- build_sobol_problem(): Build a SALib Sobol problem from Brightway ActivityParameters
- sanitize_sobol_problem(): Validate and repair Sobol problem parameters
- generate_sobol_samples(): Generate Sobol samples
- expected_sobol_sample_rows(): Calculate expected number of sample rows
- init_worker_sobol(): Sobol multiprocessing worker initialization
- sobol_worker(): Sobol multiprocessing worker processing function
- recalculate_scores_for_sobol_sample(): Re-evaluate LCA scores for a single sample (Sobol version)
- run_parallel_sobol_from_samples(): Parallel evaluation of Sobol samples
- sobol_indices_from_results(): Compute Sobol indices (S1, ST) from results
- sobol_indices_to_dataframe(): Convert Sobol indices to a tidy DataFrame
- run_full_sobol_workflow(): One-click Sobol workflow (build -> sample -> evaluate)
"""

from __future__ import annotations

import math
import multiprocessing as mp
import warnings
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from SALib.analyze import sobol as sobol_analyze
from SALib.sample import sobol as sobol_sample_mod

import bw2data as bd

from utils.lca_config import (
    VERBOSE,
    G_GROUP,
    G_PARAMETER_FORMULA_MAP,
    G_EXCHANGE_FORMULA_MAP,
    G_MULTILCA,
    G_BASE_TECH,
    G_BASE_BIO,
    G_METHODS,
    G_SOBOL_PARAMETER_NAMES,
    G_BASELINE_PARAM_CONTEXT,
    G_BASELINE_UPDATES,
    G_BASELINE_GROUPED_UPDATES,
    Timer,
    _log,
)
from utils.lca_project_setup import set_project, initialize_multilca
from utils.lca_parameters import (
    get_activity_parameters,
    build_parameter_formula_map,
    build_exchange_formula_map,
    resolve_parameter_context,
    evaluate_exchange_formulas,
)
from utils.lca_matrices import aggregate_updates_by_cell, apply_grouped_updates_to_matrices
from utils.lca_contribution import _reset_and_recalculate_mlca
from utils.lca_monte_carlo import _parse_score_key


# -----------------------------------------------------------------------------
# Sobol problem construction
# -----------------------------------------------------------------------------

def build_sobol_problem(group_name: str) -> Tuple[dict, list]:
    """
    Build a SALib Sobol problem from Brightway ActivityParameters in one group.

    Supported uncertainty types: 2 (lognormal), 3 (normal), 4 (uniform), 5 (triangular).

    Returns
    -------
    (problem, parameters_used)
        problem: SALib problem dictionary
        parameters_used: ordered list of Brightway parameters included
    """
    params = get_activity_parameters(group_name)

    names: List[str] = []
    bounds: List[list] = []
    dists: List[str] = []
    used_params: List = []

    EPS = 1e-12  # avoid c == 1 and zero-width issues

    for p in params:
        meta = p.dict
        ut = meta.get("uncertainty type", 0)

        if ut in (0, 1):
            _log(f"Sobol skipped parameter {p.name}: uncertainty type {ut}")
            continue

        if ut == 2:  # lognormal
            loc = meta.get("loc")
            scale = meta.get("scale")

            if loc is None or scale is None or not np.isfinite(loc) or not np.isfinite(scale) or scale <= 0:
                warnings.warn(
                    f"[Sobol] skipping parameter '{p.name}' "
                    f"(invalid lognormal loc/scale: loc={loc}, scale={scale})"
                )
                continue

            names.append(p.name)
            bounds.append([float(loc), float(scale)])
            dists.append("lognorm")
            used_params.append(p)

        elif ut == 3:  # normal
            loc = meta.get("loc")
            scale = meta.get("scale")

            if loc is None or scale is None or not np.isfinite(loc) or not np.isfinite(scale) or scale <= 0:
                warnings.warn(
                    f"[Sobol] skipping parameter '{p.name}' "
                    f"(invalid normal loc/scale: loc={loc}, scale={scale})"
                )
                continue

            names.append(p.name)
            bounds.append([float(loc), float(scale)])
            dists.append("norm")
            used_params.append(p)

        elif ut == 4:  # uniform
            minimum = meta.get("minimum")
            maximum = meta.get("maximum")

            if (
                minimum is None or maximum is None
                or not np.isfinite(minimum) or not np.isfinite(maximum)
                or maximum <= minimum
            ):
                warnings.warn(
                    f"[Sobol] skipping parameter '{p.name}' "
                    f"(invalid uniform bounds: minimum={minimum}, maximum={maximum})"
                )
                continue

            names.append(p.name)
            bounds.append([float(minimum), float(maximum)])
            dists.append("unif")
            used_params.append(p)

        elif ut == 5:  # triangular
            minimum = meta.get("minimum")
            maximum = meta.get("maximum")
            mode = meta.get("loc")  # Brightway often stores triangular mode in 'loc'

            if (
                minimum is None or maximum is None
                or not np.isfinite(minimum) or not np.isfinite(maximum)
                or maximum <= minimum
            ):
                warnings.warn(
                    f"[Sobol] skipping parameter '{p.name}' "
                    f"(invalid triangular bounds: minimum={minimum}, maximum={maximum})"
                )
                continue

            # SALib requires maximum >= 0 for triangular distributions
            # (its internal check: b1 = maximum, rejects b1 < 0).
            # If the entire range is negative, fall back to uniform.
            if maximum < 0:
                warnings.warn(
                    f"[Sobol] parameter '{p.name}' has negative triangular "
                    f"bounds [{minimum}, {maximum}]; falling back to uniform."
                )
                names.append(p.name)
                bounds.append([minimum, maximum])
                dists.append("unif")
                used_params.append(p)
                continue

            minimum = float(minimum)
            maximum = float(maximum)
            width = maximum - minimum

            # Missing mode -> use midpoint
            if mode is None or not np.isfinite(mode):
                mode = 0.5 * (minimum + maximum)

            mode = float(mode)

            # Clamp mode into [minimum, maximum]
            if mode < minimum:
                warnings.warn(
                    f"[Sobol] parameter '{p.name}' triangular mode {mode} < minimum {minimum}; clamped."
                )
                mode = minimum
            elif mode > maximum:
                warnings.warn(
                    f"[Sobol] parameter '{p.name}' triangular mode {mode} > maximum {maximum}; clamped."
                )
                mode = maximum

            # SALib triang expects [lower, upper, c], where c is relative peak position in [0, 1)
            c = (mode - minimum) / width

            # SALib rejects c >= 1, so slightly back off from the right endpoint
            if c >= 1.0:
                c = 1.0 - EPS
            elif c < 0.0:
                c = 0.0

            names.append(p.name)
            bounds.append([minimum, maximum, c])
            dists.append("triang")
            used_params.append(p)

        else:
            warnings.warn(
                f"[Sobol] skipping parameter '{p.name}' "
                f"(unsupported uncertainty type {ut} for Sobol)"
            )

    problem = {
        "num_vars": len(names),
        "names": names,
        "bounds": bounds,
        "dists": dists,
    }

    print(
        f"Sobol problem setup complete: "
        f"{problem['num_vars']} parameters, distributions={sorted(set(dists))}"
    )
    return problem, used_params


# -----------------------------------------------------------------------------
# Sobol problem sanitization
# -----------------------------------------------------------------------------

def sanitize_sobol_problem(problem: Mapping[str, object]) -> dict:
    """
    Validate and repair a SALib problem before Sobol sampling.

    Notes
    -----
    - SALib uses:
        * unif    -> [low, high]
        * norm    -> [mean, std]
        * lognorm -> [ln_mean, ln_std]
        * triang  -> [low, high, c], c in [0, 1)
    - For invalid parameters we replace with a very narrow legal uniform interval.
    """
    names = list(problem.get("names", []))
    bounds = list(problem.get("bounds", []))
    dists = list(problem.get("dists", ["unif"] * len(bounds)))

    if len(names) != len(bounds):
        raise ValueError(
            f"Invalid Sobol problem: len(names)={len(names)} != len(bounds)={len(bounds)}"
        )

    if len(dists) != len(bounds):
        raise ValueError(
            f"Invalid Sobol problem: len(dists)={len(dists)} != len(bounds)={len(bounds)}"
        )

    EPS = 1e-12

    def _tiny_legal_uniform(center: float) -> Tuple[list, str]:
        if not np.isfinite(center):
            center = 0.0
        delta = max(abs(center) * 1e-12, 1e-12)
        return [center - delta, center + delta], "unif"

    new_names: List[str] = []
    new_bounds: List[list] = []
    new_dists: List[str] = []

    for i, (name, b, dist) in enumerate(zip(names, bounds, dists)):
        pname = name if name is not None else f"param_{i}"

        try:
            if dist == "triang":
                if not isinstance(b, (list, tuple)) or len(b) != 3:
                    raise ValueError("triang bounds must be [low, high, c]")

                low, high, c = b
                low = float(low)
                high = float(high)
                c = float(c)

                if not np.isfinite(low) or not np.isfinite(high) or not np.isfinite(c):
                    raise ValueError("triang bounds contain non-finite values")
                if high <= low:
                    raise ValueError("triang requires high > low")
                if high < 0:
                    raise ValueError("triang requires high >= 0 in present SALib")
                if c < 0.0:
                    raise ValueError("triang requires c >= 0")
                if c >= 1.0:
                    c = 1.0 - EPS

                # SALib rejects triangular distributions where high (b1) < 0.
                # Fall back to uniform if the entire range is negative.
                if high < 0:
                    warnings.warn(
                        f"[Sobol sanitize] Parameter '{pname}' has negative triangular "
                        f"bounds [low={low}, high={high}]; falling back to uniform."
                    )
                    new_names.append(pname)
                    new_bounds.append([low, high])
                    new_dists.append("unif")
                else:
                    new_names.append(pname)
                    new_bounds.append([low, high, c])
                    new_dists.append("triang")

            elif dist == "unif":
                if not isinstance(b, (list, tuple)) or len(b) != 2:
                    raise ValueError("unif bounds must be [low, high]")

                low, high = float(b[0]), float(b[1])

                if not np.isfinite(low) or not np.isfinite(high):
                    raise ValueError("unif bounds contain non-finite values")
                if high <= low:
                    raise ValueError("unif requires high > low")

                new_names.append(pname)
                new_bounds.append([low, high])
                new_dists.append("unif")

            elif dist in ("norm", "lognorm"):
                if not isinstance(b, (list, tuple)) or len(b) != 2:
                    raise ValueError(f"{dist} bounds must be [loc, scale]")

                loc, scale = float(b[0]), float(b[1])

                if not np.isfinite(loc) or not np.isfinite(scale):
                    raise ValueError(f"{dist} bounds contain non-finite values")
                if scale <= 0:
                    raise ValueError(f"{dist} requires scale > 0")

                new_names.append(pname)
                new_bounds.append([loc, scale])
                new_dists.append(dist)

            else:
                raise ValueError(f"unsupported SALib dist '{dist}'")

        except Exception as e:
            if isinstance(b, (list, tuple)) and len(b) >= 1:
                try:
                    center = float(np.mean([x for x in b[:2] if np.isfinite(x)]))
                except Exception:
                    center = 0.0
            else:
                center = 0.0

            fixed_bounds, fixed_dist = _tiny_legal_uniform(center)
            warnings.warn(
                f"[Sobol sanitize] Parameter '{pname}' invalid ({e}); "
                f"replaced with narrow legal uniform interval {fixed_bounds}"
            )
            new_names.append(pname)
            new_bounds.append(fixed_bounds)
            new_dists.append(fixed_dist)

    problem_clean = dict(problem)
    problem_clean["names"] = new_names
    problem_clean["bounds"] = new_bounds
    problem_clean["dists"] = new_dists
    problem_clean["num_vars"] = len(new_names)

    return problem_clean


# -----------------------------------------------------------------------------
# Sobol sampling
# -----------------------------------------------------------------------------

def generate_sobol_samples(
    sobol_problem: Mapping[str, object],
    N: int = 128,
    calc_second_order: bool = False,
    scramble: bool = True,
    seed: Optional[int] = 42,
) -> np.ndarray:
    """
    Generate Sobol samples safely.

    Returns shape:
    - (N * (D + 2), D) if calc_second_order=False
    - (N * (2D + 2), D) if calc_second_order=True
    """
    sobol_problem_clean = sanitize_sobol_problem(sobol_problem)

    if sobol_problem_clean["num_vars"] == 0:
        raise ValueError("No valid parameters available for Sobol sampling.")

    samples = sobol_sample_mod.sample(
        sobol_problem_clean,
        N,
        calc_second_order=calc_second_order,
        scramble=scramble,
        seed=seed,
    )
    return samples


def expected_sobol_sample_rows(
    problem: Mapping[str, object],
    N: int,
    calc_second_order: bool = False,
) -> int:
    """Expected number of Sobol sample rows for a SALib problem."""
    D = int(problem["num_vars"])
    return N * (2 * D + 2) if calc_second_order else N * (D + 2)


# -----------------------------------------------------------------------------
# Sobol multiprocessing workers
# -----------------------------------------------------------------------------

def init_worker_sobol(
    project_name: str,
    group_name: str,
    fus_keys: Sequence[Tuple[str, str]],
    methods: Sequence[tuple],
    sobol_parameter_names: Sequence[str],
    verbose: bool = False,
):
    """Initialize a Sobol multiprocessing worker."""
    global VERBOSE
    global G_GROUP
    global G_PARAMETER_FORMULA_MAP
    global G_EXCHANGE_FORMULA_MAP
    global G_MULTILCA
    global G_BASE_TECH
    global G_BASE_BIO
    global G_METHODS
    global G_SOBOL_PARAMETER_NAMES
    global G_BASELINE_PARAM_CONTEXT
    global G_BASELINE_UPDATES
    global G_BASELINE_GROUPED_UPDATES

    VERBOSE = verbose

    with Timer("sobol_worker_init"):
        set_project(project_name)
        G_GROUP = group_name
        G_METHODS = list(methods)
        G_SOBOL_PARAMETER_NAMES = list(sobol_parameter_names)

        fus = [bd.get_node(database=db, code=code) for db, code in fus_keys]

        G_PARAMETER_FORMULA_MAP = build_parameter_formula_map(group_name)

        setup = initialize_multilca(fus, G_METHODS)
        G_MULTILCA = setup.mlca
        G_BASE_TECH = setup.base_tech
        G_BASE_BIO = setup.base_bio

        G_EXCHANGE_FORMULA_MAP = build_exchange_formula_map(group_name, G_MULTILCA)

        G_BASELINE_PARAM_CONTEXT = resolve_parameter_context(
            G_GROUP,
            sampled_values_by_name={},
            parameter_formula_map=G_PARAMETER_FORMULA_MAP,
        )
        G_BASELINE_UPDATES = evaluate_exchange_formulas(
            G_EXCHANGE_FORMULA_MAP,
            G_BASELINE_PARAM_CONTEXT,
        )
        G_BASELINE_GROUPED_UPDATES = aggregate_updates_by_cell(G_BASELINE_UPDATES)

        _log(
            f"initialized Sobol worker group={group_name} "
            f"fus={len(fus)} methods={len(methods)} "
            f"param_formulas={len(G_PARAMETER_FORMULA_MAP)} "
            f"exchange_formulas={len(G_EXCHANGE_FORMULA_MAP)} "
            f"sobol_params={len(G_SOBOL_PARAMETER_NAMES)}"
        )


def recalculate_scores_for_sobol_sample(
    sample_index: int,
    sampled_values: Mapping[str, float],
) -> List[dict]:
    """Evaluate one Sobol sample row and return tidy score rows."""
    param_context = resolve_parameter_context(
        G_GROUP,
        sampled_values_by_name=sampled_values,
        parameter_formula_map=G_PARAMETER_FORMULA_MAP,
    )

    new_updates = evaluate_exchange_formulas(G_EXCHANGE_FORMULA_MAP, param_context)
    new_grouped_updates = aggregate_updates_by_cell(new_updates)

    tech_matrix, bio_matrix = apply_grouped_updates_to_matrices(
        base_tech=G_BASE_TECH,
        base_bio=G_BASE_BIO,
        baseline_grouped_updates=G_BASELINE_GROUPED_UPDATES,
        new_grouped_updates=new_grouped_updates,
    )

    _reset_and_recalculate_mlca(G_MULTILCA, tech_matrix, bio_matrix)

    rows = []
    for key, score in G_MULTILCA.scores.items():
        fu_label, method = _parse_score_key(key)
        rows.append(
            {
                "Sample index": sample_index,
                "Method": method[2] if isinstance(method, tuple) and len(method) > 2 else str(method),
                "Method full": method,
                "Functional unit": fu_label,
                "LCA Score": float(score),
            }
        )
    return rows


def sobol_worker(sample_chunk: Sequence[Tuple[int, np.ndarray]]) -> List[dict]:
    """Multiprocessing worker for Sobol sample chunks."""
    rows = []
    with Timer(f"sobol chunk size={len(sample_chunk)}"):
        for sample_index, sample_row in sample_chunk:
            sampled_values = {
                G_SOBOL_PARAMETER_NAMES[j]: float(sample_row[j])
                for j in range(len(G_SOBOL_PARAMETER_NAMES))
            }
            rows.extend(
                recalculate_scores_for_sobol_sample(
                    sample_index=sample_index,
                    sampled_values=sampled_values,
                )
            )
    return rows


# -----------------------------------------------------------------------------
# Parallel Sobol evaluation
# -----------------------------------------------------------------------------

def run_parallel_sobol_from_samples(
    project_name: str,
    group_name: str,
    functional_units: Sequence,
    methods: Sequence[tuple],
    sobol_problem: Mapping[str, object],
    sobol_samples: np.ndarray,
    n_workers: Optional[int] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Evaluate an existing Sobol sample matrix in parallel.

    Returns a tidy DataFrame with one row per sample/FU/method combination.
    """
    set_project(project_name)

    if sobol_samples.ndim != 2:
        raise ValueError("sobol_samples must be a 2D array")

    if sobol_samples.shape[1] != len(sobol_problem["names"]):
        raise ValueError(
            f"sobol_samples has {sobol_samples.shape[1]} columns but "
            f"sobol_problem has {len(sobol_problem['names'])} names"
        )

    indexed_samples = [(i, sobol_samples[i, :]) for i in range(sobol_samples.shape[0])]

    if n_workers is None:
        n_workers = min(max(1, mp.cpu_count() - 1), len(indexed_samples))

    chunk_size = max(1, math.ceil(len(indexed_samples) / n_workers))
    chunks = [
        indexed_samples[i : i + chunk_size]
        for i in range(0, len(indexed_samples), chunk_size)
    ]

    fus_keys = [(fu["database"], fu["code"]) for fu in functional_units]

    ctx = mp.get_context("spawn")
    with ctx.Pool(
        processes=len(chunks),
        initializer=init_worker_sobol,
        initargs=(
            project_name,
            group_name,
            fus_keys,
            methods,
            sobol_problem["names"],
            verbose,
        ),
    ) as pool:
        nested_results = pool.map(sobol_worker, chunks)

    return pd.DataFrame([item for chunk in nested_results for item in chunk])


# -----------------------------------------------------------------------------
# Sobol index extraction
# -----------------------------------------------------------------------------

def sobol_indices_from_results(
    sobol_problem: Mapping[str, object],
    sobol_results_df: pd.DataFrame,
    functional_unit: str,
    method: tuple,
    calc_second_order: bool = False,
    print_to_console: bool = False,
) -> dict:
    """
    Compute Sobol indices for one (functional unit, method) combination.

    Returns the raw SALib Si dict (S1, S1_conf, ST, ST_conf, etc.).

    Notes
    -----
    The 'Functional unit' column in sobol_results_df may have a trailing
    ``__N`` suffix (e.g. "P. Bio-phenol assembly__0") appended by MultiLCA.
    This function tolerates both exact matches and prefix matches.
    """
    # Try exact match first
    fu_mask = sobol_results_df["Functional unit"] == functional_unit
    # If no exact match, try prefix match (handles "__0", "__1", ... suffix)
    if not fu_mask.any():
        prefix = functional_unit + "__"
        fu_mask = sobol_results_df["Functional unit"].str.startswith(prefix)

    sub = sobol_results_df[
        fu_mask
        & (sobol_results_df["Method full"].apply(lambda x: x == method))
    ].copy()

    sub = sub.sort_values("Sample index")
    Y = sub["LCA Score"].to_numpy(dtype=float)

    if len(Y) == 0:
        raise ValueError(
            f"No Sobol results found for functional unit '{functional_unit}' and method {method}"
        )

    Si = sobol_analyze.analyze(
        sobol_problem,
        Y,
        calc_second_order=calc_second_order,
        print_to_console=print_to_console,
    )
    return Si


def sobol_indices_to_dataframe(
    Si: Mapping[str, np.ndarray],
    sobol_problem: Mapping[str, object],
) -> pd.DataFrame:
    """Convert SALib Sobol index output to a tidy DataFrame sorted by ST."""
    rows = []
    names = list(sobol_problem["names"])

    for i, name in enumerate(names):
        rows.append(
            {
                "Parameter": name,
                "S1": float(Si["S1"][i]) if Si["S1"][i] is not None else np.nan,
                "S1_conf": float(Si["S1_conf"][i]) if Si["S1_conf"][i] is not None else np.nan,
                "ST": float(Si["ST"][i]) if Si["ST"][i] is not None else np.nan,
                "ST_conf": float(Si["ST_conf"][i]) if Si["ST_conf"][i] is not None else np.nan,
            }
        )

    return pd.DataFrame(rows).sort_values("ST", ascending=False).reset_index(drop=True)


# -----------------------------------------------------------------------------
# Convenience: full Sobol workflow
# -----------------------------------------------------------------------------

def run_full_sobol_workflow(
    project_name: str,
    group_name: str,
    functional_units: Sequence,
    methods: Sequence[tuple],
    N: int = 128,
    calc_second_order: bool = False,
    n_workers: Optional[int] = None,
    scramble: bool = True,
    seed: Optional[int] = 42,
    verbose: bool = False,
) -> Dict[str, object]:
    """
    Convenience notebook-facing Sobol workflow:
    1. build problem
    2. generate Sobol samples
    3. evaluate in parallel

    Returns dict with:
        - sobol_problem
        - sobol_parameters
        - sobol_samples
        - sobol_results
    """
    sobol_problem, sobol_parameters = build_sobol_problem(group_name)
    sobol_samples = generate_sobol_samples(
        sobol_problem=sobol_problem,
        N=N,
        calc_second_order=calc_second_order,
        scramble=scramble,
        seed=seed,
    )
    sobol_results = run_parallel_sobol_from_samples(
        project_name=project_name,
        group_name=group_name,
        functional_units=functional_units,
        methods=methods,
        sobol_problem=sobol_problem,
        sobol_samples=sobol_samples,
        n_workers=n_workers,
        verbose=verbose,
    )

    return {
        "sobol_problem": sobol_problem,
        "sobol_parameters": sobol_parameters,
        "sobol_samples": sobol_samples,
        "sobol_results": sobol_results,
    }
