"""
lca_fast.py -- FAST (Fourier Amplitude Sensitivity Test) global sensitivity analysis

This module provides:
- build_fast_problem(): Build a SALib FAST problem from Brightway ActivityParameters
- sanitize_fast_problem(): Validate and repair FAST problem parameters
- generate_fast_samples(): Generate FAST samples using SALib
- expected_fast_sample_rows(): Expected number of FAST sample rows
- init_worker_fast(): FAST multiprocessing worker initialization
- fast_worker(): FAST multiprocessing worker processing function
- recalculate_scores_for_fast_sample(): Re-evaluate LCA scores for a single sample (FAST version)
- run_parallel_fast_from_samples(): Parallel evaluation of FAST samples
- fast_indices_from_results(): Compute FAST indices (S1, ST) from results
- fast_indices_to_dataframe(): Convert FAST indices to a tidy DataFrame
- run_full_fast_workflow(): One-click FAST workflow (build -> sample -> evaluate -> analyze)

Notes
-----
FAST (Fourier Amplitude Sensitivity Test) decomposes the output variance into
contributions from each input parameter using Fourier analysis of the model
output sampled along a specially designed frequency-based space-filling curve.

Two variants are supported via SALib:
- FAST    : classic Fourier Amplitude Sensitivity Test -- S1 (first-order) only
- RBD-FAST: Random Balance Designs FAST -- S1 + ST (total-order) with fewer samples

Supported Brightway uncertainty types: 2 (lognormal), 3 (normal),
                                         4 (uniform), 5 (triangular).

SALib distribution format
-------------------------
- unif    : bounds = [low, high]
- norm    : bounds = [mean, std]
- lognorm : bounds = [ln_mean, ln_std]
- triang  : bounds = [low, high, c],  c = relative peak in [0, 1)
"""

from __future__ import annotations

import math
import multiprocessing as mp
import warnings
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

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

# ---------------------------------------------------------------------------
# Worker-level global for FAST parameter names (separate from Sobol globals)
# ---------------------------------------------------------------------------
G_FAST_PARAMETER_NAMES: Optional[List[str]] = None


# =============================================================================
# FAST problem construction
# =============================================================================

def build_fast_problem(group_name: str) -> Tuple[dict, list]:
    """
    Build a SALib FAST/RBD-FAST problem from Brightway ActivityParameters.

    Supported uncertainty types: 2 (lognormal), 3 (normal), 4 (uniform), 5 (triangular).

    Returns
    -------
    (problem, parameters_used)
        problem : SALib problem dictionary
        parameters_used : ordered list of Brightway parameters included
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
            _log(f"FAST skipped parameter {p.name}: uncertainty type {ut} (no uncertainty)")
            continue

        if ut == 2:  # lognormal
            loc = meta.get("loc")
            scale = meta.get("scale")
            if (
                loc is None or scale is None
                or not np.isfinite(loc) or not np.isfinite(scale)
                or scale <= 0
            ):
                warnings.warn(
                    f"[FAST] skipping parameter '{p.name}' "
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
            if (
                loc is None or scale is None
                or not np.isfinite(loc) or not np.isfinite(scale)
                or scale <= 0
            ):
                warnings.warn(
                    f"[FAST] skipping parameter '{p.name}' "
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
                    f"[FAST] skipping parameter '{p.name}' "
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
            mode = meta.get("loc")

            if (
                minimum is None or maximum is None
                or not np.isfinite(minimum) or not np.isfinite(maximum)
                or maximum <= minimum
            ):
                warnings.warn(
                    f"[FAST] skipping parameter '{p.name}' "
                    f"(invalid triangular bounds: minimum={minimum}, maximum={maximum})"
                )
                continue

            # SALib triang requires maximum >= 0
            if maximum < 0:
                warnings.warn(
                    f"[FAST] parameter '{p.name}' has negative triangular "
                    f"bounds [{minimum}, {maximum}]; falling back to uniform."
                )
                names.append(p.name)
                bounds.append([float(minimum), float(maximum)])
                dists.append("unif")
                used_params.append(p)
                continue

            minimum = float(minimum)
            maximum = float(maximum)
            width = maximum - minimum

            if mode is None or not np.isfinite(mode):
                mode = 0.5 * (minimum + maximum)
            mode = float(mode)

            mode = float(np.clip(mode, minimum, maximum))
            c = (mode - minimum) / width
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
                f"[FAST] skipping parameter '{p.name}' "
                f"(unsupported uncertainty type {ut} for FAST)"
            )

    problem = {
        "num_vars": len(names),
        "names": names,
        "bounds": bounds,
        "dists": dists,
    }

    print(
        f"FAST problem setup complete: "
        f"{problem['num_vars']} parameters, distributions={sorted(set(dists))}"
    )
    return problem, used_params


# =============================================================================
# FAST problem sanitization
# =============================================================================

def sanitize_fast_problem(problem: Mapping[str, object]) -> dict:
    """
    Validate and repair a SALib problem before FAST sampling.

    SALib uses the same problem format for FAST and Sobol, so the same
    validation rules apply:
        unif    -> [low, high]
        norm    -> [mean, std]
        lognorm -> [ln_mean, ln_std]
        triang  -> [low, high, c],  c in [0, 1)

    Invalid parameters are replaced with a very narrow legal uniform interval
    so that sampling does not abort.
    """
    names = list(problem.get("names", []))
    bounds = list(problem.get("bounds", []))
    dists = list(problem.get("dists", ["unif"] * len(bounds)))

    if len(names) != len(bounds):
        raise ValueError(
            f"Invalid FAST problem: len(names)={len(names)} != len(bounds)={len(bounds)}"
        )
    if len(dists) != len(bounds):
        raise ValueError(
            f"Invalid FAST problem: len(dists)={len(dists)} != len(bounds)={len(bounds)}"
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
                low, high, c = float(b[0]), float(b[1]), float(b[2])
                if not all(np.isfinite(v) for v in (low, high, c)):
                    raise ValueError("triang bounds contain non-finite values")
                if high <= low:
                    raise ValueError("triang requires high > low")
                if high < 0:
                    raise ValueError("triang requires high >= 0 in present SALib")
                if c < 0.0:
                    raise ValueError("triang requires c >= 0")
                c = min(c, 1.0 - EPS)
                new_names.append(pname)
                new_bounds.append([low, high, c])
                new_dists.append("triang")

            elif dist == "unif":
                if not isinstance(b, (list, tuple)) or len(b) != 2:
                    raise ValueError("unif bounds must be [low, high]")
                low, high = float(b[0]), float(b[1])
                if not (np.isfinite(low) and np.isfinite(high)):
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
                if not (np.isfinite(loc) and np.isfinite(scale)):
                    raise ValueError(f"{dist} bounds contain non-finite values")
                if scale <= 0:
                    raise ValueError(f"{dist} requires scale > 0")
                new_names.append(pname)
                new_bounds.append([loc, scale])
                new_dists.append(dist)

            else:
                raise ValueError(f"unsupported SALib dist '{dist}'")

        except Exception as e:
            try:
                center = float(np.mean([x for x in list(b)[:2] if np.isfinite(x)]))
            except Exception:
                center = 0.0
            fixed_bounds, fixed_dist = _tiny_legal_uniform(center)
            warnings.warn(
                f"[FAST sanitize] Parameter '{pname}' invalid ({e}); "
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


# =============================================================================
# FAST sampling
# =============================================================================

def generate_fast_samples(
    fast_problem: Mapping[str, object],
    M: int = 4,
    seed: Optional[int] = 42,
    method: str = "fast",
) -> np.ndarray:
    """
    Generate FAST samples safely.

    Parameters
    ----------
    fast_problem : dict
        SALib problem dictionary (output of ``build_fast_problem``).
    M : int, default 4
        Interference factor for FAST / RBD-FAST. The number of harmonics.
        Total sample size = ``4 * M^2 * num_vars + 1`` (FAST) or a multiple
        of ``num_vars`` (RBD-FAST, set by ``N`` parameter which defaults to
        ``65`` multiples inside SALib). For RBD-FAST M is the number of
        harmonics for the total-effect calculation.
    seed : int, optional
        Random seed for reproducibility (RBD-FAST only; classic FAST is
        deterministic).
    method : {'fast', 'rbd_fast'}
        Sampling method:
        - ``'fast'``    : classic FAST (S1 only, deterministic sample).
                          Uses ``SALib.sample.fast_sampler``.
        - ``'rbd_fast'``: Random Balance Designs FAST (S1 + ST, random sample).
                          Uses ``SALib.sample.latin``.

    Returns
    -------
    np.ndarray
        Sample matrix of shape (n_samples, num_vars).

    Notes
    -----
    Classic FAST sample size: ``N_per_var * num_vars`` where
    ``N_per_var = 4 * M**2 + 1`` (the minimum recommended by Saltelli et al.).
    RBD-FAST sample size: arbitrary (default ≈ 1000 in SALib; users can pass
    ``N`` to the underlying sampler by using ``generate_rbd_fast_samples``).
    """
    fast_problem_clean = sanitize_fast_problem(fast_problem)

    if fast_problem_clean["num_vars"] == 0:
        raise ValueError("No valid parameters available for FAST sampling.")

    method_lower = method.lower().replace("-", "_")

    if method_lower == "fast":
        from SALib.sample.fast_sampler import sample as fast_sample_fn
        samples = fast_sample_fn(fast_problem_clean, M=M)

    elif method_lower == "rbd_fast":
        # RBD-FAST uses a Latin Hypercube sample (random)
        from SALib.sample.latin import sample as lhs_sample_fn
        # SALib RBD-FAST recommend N >= 65 per variable; use 4*M^2+1 per var
        N_per_var = max(65, 4 * M ** 2 + 1)
        n_total = N_per_var * fast_problem_clean["num_vars"]
        rng = np.random.default_rng(seed)
        # SALib latin sampler does not accept seed directly in older versions;
        # we set the global numpy seed as fallback.
        np.random.seed(rng.integers(0, 2**31))
        samples = lhs_sample_fn(fast_problem_clean, N=n_total)

    else:
        raise ValueError(
            f"Unknown FAST method '{method}'. Choose 'fast' or 'rbd_fast'."
        )

    return samples


def expected_fast_sample_rows(
    problem: Mapping[str, object],
    M: int = 4,
    method: str = "fast",
    N_rbd: Optional[int] = None,
) -> int:
    """
    Expected number of FAST sample rows for a SALib problem.

    Parameters
    ----------
    problem : dict
        SALib problem dictionary.
    M : int
        Interference factor.
    method : {'fast', 'rbd_fast'}
        Sampling method.
    N_rbd : int, optional
        For RBD-FAST, the total number of samples (overrides auto-calculation).

    Returns
    -------
    int
        Expected number of sample rows.
    """
    D = int(problem["num_vars"])
    method_lower = method.lower().replace("-", "_")
    if method_lower == "fast":
        # SALib FAST: (4*M^2 + 1) * D
        return (4 * M ** 2 + 1) * D
    elif method_lower == "rbd_fast":
        if N_rbd is not None:
            return int(N_rbd)
        return max(65, 4 * M ** 2 + 1) * D
    else:
        raise ValueError(f"Unknown FAST method '{method}'.")


# =============================================================================
# FAST multiprocessing workers
# =============================================================================

def init_worker_fast(
    project_name: str,
    group_name: str,
    fus_keys: Sequence[Tuple[str, str]],
    methods: Sequence[tuple],
    fast_parameter_names: Sequence[str],
    verbose: bool = False,
) -> None:
    """Initialize a FAST multiprocessing worker."""
    global VERBOSE
    global G_GROUP
    global G_PARAMETER_FORMULA_MAP
    global G_EXCHANGE_FORMULA_MAP
    global G_MULTILCA
    global G_BASE_TECH
    global G_BASE_BIO
    global G_METHODS
    global G_FAST_PARAMETER_NAMES
    global G_BASELINE_PARAM_CONTEXT
    global G_BASELINE_UPDATES
    global G_BASELINE_GROUPED_UPDATES

    VERBOSE = verbose

    with Timer("fast_worker_init"):
        set_project(project_name)
        G_GROUP = group_name
        G_METHODS = list(methods)
        G_FAST_PARAMETER_NAMES = list(fast_parameter_names)

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
            f"initialized FAST worker group={group_name} "
            f"fus={len(fus)} methods={len(methods)} "
            f"param_formulas={len(G_PARAMETER_FORMULA_MAP)} "
            f"exchange_formulas={len(G_EXCHANGE_FORMULA_MAP)} "
            f"fast_params={len(G_FAST_PARAMETER_NAMES)}"
        )


def recalculate_scores_for_fast_sample(
    sample_index: int,
    sampled_values: Mapping[str, float],
) -> List[dict]:
    """Evaluate one FAST sample row and return tidy score rows."""
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


def fast_worker(sample_chunk: Sequence[Tuple[int, np.ndarray]]) -> List[dict]:
    """Multiprocessing worker for FAST sample chunks."""
    rows = []
    with Timer(f"fast chunk size={len(sample_chunk)}"):
        for sample_index, sample_row in sample_chunk:
            sampled_values = {
                G_FAST_PARAMETER_NAMES[j]: float(sample_row[j])
                for j in range(len(G_FAST_PARAMETER_NAMES))
            }
            rows.extend(
                recalculate_scores_for_fast_sample(
                    sample_index=sample_index,
                    sampled_values=sampled_values,
                )
            )
    return rows


# =============================================================================
# Parallel FAST evaluation
# =============================================================================

def run_parallel_fast_from_samples(
    project_name: str,
    group_name: str,
    functional_units: Sequence,
    methods: Sequence[tuple],
    fast_problem: Mapping[str, object],
    fast_samples: np.ndarray,
    n_workers: Optional[int] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Evaluate an existing FAST sample matrix in parallel.

    Parameters
    ----------
    project_name : str
        Brightway project name.
    group_name : str
        ActivityParameter group name.
    functional_units : list
        Brightway activity objects (functional units).
    methods : list of tuple
        Impact category method tuples.
    fast_problem : dict
        SALib problem dictionary (output of ``build_fast_problem``).
    fast_samples : np.ndarray
        Sample matrix of shape (n_samples, num_vars) from
        ``generate_fast_samples``.
    n_workers : int, optional
        Number of parallel worker processes. Defaults to ``cpu_count - 1``.
    verbose : bool
        Enable per-worker verbose logging.

    Returns
    -------
    pd.DataFrame
        Tidy DataFrame with columns:
        ``Sample index, Method, Method full, Functional unit, LCA Score``
    """
    set_project(project_name)

    if fast_samples.ndim != 2:
        raise ValueError("fast_samples must be a 2D array")
    if fast_samples.shape[1] != len(fast_problem["names"]):
        raise ValueError(
            f"fast_samples has {fast_samples.shape[1]} columns but "
            f"fast_problem has {len(fast_problem['names'])} names"
        )

    indexed_samples = [(i, fast_samples[i, :]) for i in range(fast_samples.shape[0])]

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
        initializer=init_worker_fast,
        initargs=(
            project_name,
            group_name,
            fus_keys,
            methods,
            fast_problem["names"],
            verbose,
        ),
    ) as pool:
        nested_results = pool.map(fast_worker, chunks)

    return pd.DataFrame([item for chunk in nested_results for item in chunk])


# =============================================================================
# FAST index extraction
# =============================================================================

def fast_indices_from_results(
    fast_problem: Mapping[str, object],
    fast_results_df: pd.DataFrame,
    functional_unit: str,
    method: tuple,
    fast_samples: np.ndarray,
    M: int = 4,
    analysis_method: str = "fast",
    print_to_console: bool = False,
) -> dict:
    """
    Compute FAST sensitivity indices for one (functional unit, method) pair.

    Parameters
    ----------
    fast_problem : dict
        SALib problem dictionary.
    fast_results_df : pd.DataFrame
        Output of ``run_parallel_fast_from_samples``.
    functional_unit : str
        Functional unit label (exact or prefix match supported).
    method : tuple
        Full method tuple.
    fast_samples : np.ndarray
        The sample matrix used for the evaluation (needed for RBD-FAST).
    M : int
        Interference factor (number of harmonics).  Must match the value
        used in ``generate_fast_samples``.
    analysis_method : {'fast', 'rbd_fast'}
        Analysis method to use.
    print_to_console : bool
        Whether to print indices to console.

    Returns
    -------
    dict
        Raw SALib Si dict.  Classic FAST: ``{'S1', 'S1_conf'}``.
        RBD-FAST: ``{'S1', 'S1_conf', 'ST', 'ST_conf'}``.

    Notes
    -----
    The ``Functional unit`` column may carry a ``__N`` suffix from MultiLCA.
    Both exact matches and prefix matches (``name__0``, ``name__1``, ...) are
    handled automatically.
    """
    # Try exact match, then prefix match
    fu_mask = fast_results_df["Functional unit"] == functional_unit
    if not fu_mask.any():
        prefix = functional_unit + "__"
        fu_mask = fast_results_df["Functional unit"].str.startswith(prefix)

    sub = fast_results_df[
        fu_mask
        & (fast_results_df["Method full"].apply(lambda x: x == method))
    ].copy()

    sub = sub.sort_values("Sample index")
    Y = sub["LCA Score"].to_numpy(dtype=float)

    if len(Y) == 0:
        raise ValueError(
            f"No FAST results found for functional unit '{functional_unit}' "
            f"and method {method}"
        )

    method_lower = analysis_method.lower().replace("-", "_")

    if method_lower == "fast":
        from SALib.analyze.fast import analyze as fast_analyze_fn
        Si = fast_analyze_fn(
            fast_problem,
            Y,
            M=M,
            print_to_console=print_to_console,
        )

    elif method_lower == "rbd_fast":
        from SALib.analyze.rbd_fast import analyze as rbd_analyze_fn
        Si = rbd_analyze_fn(
            fast_problem,
            fast_samples,
            Y,
            M=M,
            print_to_console=print_to_console,
        )

    else:
        raise ValueError(
            f"Unknown analysis method '{analysis_method}'. "
            "Choose 'fast' or 'rbd_fast'."
        )

    return Si


def fast_indices_to_dataframe(
    Si: Mapping[str, np.ndarray],
    fast_problem: Mapping[str, object],
    analysis_method: str = "fast",
) -> pd.DataFrame:
    """
    Convert SALib FAST index output to a tidy DataFrame.

    Parameters
    ----------
    Si : dict
        SALib sensitivity index dictionary.
    fast_problem : dict
        SALib problem dictionary.
    analysis_method : {'fast', 'rbd_fast'}
        Used to determine which keys to expect in Si.

    Returns
    -------
    pd.DataFrame
        Tidy DataFrame sorted by ``S1`` descending.
        Classic FAST columns : ``Parameter, S1, S1_conf``
        RBD-FAST columns     : ``Parameter, S1, S1_conf, ST, ST_conf``
    """
    names = list(fast_problem["names"])
    rows = []

    has_total = "ST" in Si and Si["ST"] is not None

    for i, name in enumerate(names):
        row: dict = {
            "Parameter": name,
            "S1": float(Si["S1"][i]) if Si["S1"][i] is not None else np.nan,
            "S1_conf": (
                float(Si["S1_conf"][i])
                if "S1_conf" in Si and Si["S1_conf"][i] is not None
                else np.nan
            ),
        }
        if has_total:
            row["ST"] = float(Si["ST"][i]) if Si["ST"][i] is not None else np.nan
            row["ST_conf"] = (
                float(Si["ST_conf"][i])
                if "ST_conf" in Si and Si["ST_conf"][i] is not None
                else np.nan
            )
        rows.append(row)

    df = pd.DataFrame(rows)
    sort_col = "ST" if has_total and "ST" in df.columns else "S1"
    return df.sort_values(sort_col, ascending=False).reset_index(drop=True)


# =============================================================================
# Convenience: full FAST workflow
# =============================================================================

def run_full_fast_workflow(
    project_name: str,
    group_name: str,
    functional_units: Sequence,
    methods: Sequence[tuple],
    M: int = 4,
    n_workers: Optional[int] = None,
    seed: Optional[int] = 42,
    verbose: bool = False,
    method: str = "fast",
) -> Dict[str, object]:
    """
    Convenience one-click FAST GSA workflow:
        1. build problem
        2. generate FAST samples
        3. evaluate in parallel
        4. (indices must be computed separately per FU/method via
           ``fast_indices_from_results``)

    Parameters
    ----------
    project_name : str
        Brightway project name.
    group_name : str
        ActivityParameter group name.
    functional_units : list
        Brightway activity objects.
    methods : list of tuple
        Impact category method tuples.
    M : int, default 4
        FAST interference factor / number of harmonics.
    n_workers : int, optional
        Parallel worker count.
    seed : int, optional
        Random seed (for RBD-FAST only).
    verbose : bool
        Worker verbose logging.
    method : {'fast', 'rbd_fast'}
        Sampling and analysis method.

    Returns
    -------
    dict with keys:
        - ``fast_problem``    : SALib problem dict
        - ``fast_parameters`` : list of Brightway parameters used
        - ``fast_samples``    : np.ndarray sample matrix
        - ``fast_results``    : pd.DataFrame of LCA scores
        - ``method``          : method string used
    """
    fast_problem, fast_parameters = build_fast_problem(group_name)
    fast_samples = generate_fast_samples(
        fast_problem=fast_problem,
        M=M,
        seed=seed,
        method=method,
    )
    fast_results = run_parallel_fast_from_samples(
        project_name=project_name,
        group_name=group_name,
        functional_units=functional_units,
        methods=methods,
        fast_problem=fast_problem,
        fast_samples=fast_samples,
        n_workers=n_workers,
        verbose=verbose,
    )

    return {
        "fast_problem": fast_problem,
        "fast_parameters": fast_parameters,
        "fast_samples": fast_samples,
        "fast_results": fast_results,
        "method": method,
    }
