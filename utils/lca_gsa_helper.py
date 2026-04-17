"""
lca_gsa_helper.py -- Backward-compatible re-export layer

The original monolithic file has been split into the following sub-modules:
    lca_config         -- Global configuration, logging, Timer, LCASetup dataclass
    lca_project_setup  -- Project initialization, MultiLCA construction, deterministic LCA
    lca_parameters     -- Parameter queries, formula resolution, validation, exchange mapping
    lca_matrices       -- Sparse matrix update utilities
    lca_contribution   -- Contribution analysis, OAT sensitivity, analytical variance
    lca_monte_carlo    -- Monte Carlo sampling and parallel evaluation
    lca_sobol          -- Sobol global sensitivity analysis

If you previously used ``from lca_gsa_helper import X``, it still works.
New code should import from the corresponding sub-module directly.
"""

# -- lca_config ----------------------------------------------------------------
from utils.lca_config import (                          # noqa: F401
    VERBOSE,
    G_GROUP,
    G_PARAMETER_FORMULA_MAP,
    G_EXCHANGE_FORMULA_MAP,
    G_MULTILCA,
    G_BASE_TECH,
    G_BASE_BIO,
    G_METHODS,
    G_ACTIVITY_PARAMETER_NAMES,
    G_SOBOL_PARAMETER_NAMES,
    G_BASELINE_PARAM_CONTEXT,
    G_BASELINE_UPDATES,
    G_BASELINE_GROUPED_UPDATES,
    LCASetup,
    _log,
    Timer,
)

# -- lca_project_setup ---------------------------------------------------------
from utils.lca_project_setup import (                      # noqa: F401
    set_project,
    build_fu_labels,
    build_demands,
    initialize_multilca,
    flatten_multilca_scores,
    deterministic_lca,
)

# -- lca_parameters ------------------------------------------------------------
from utils.lca_parameters import (                         # noqa: F401
    get_activity_parameters,
    build_parameter_formula_map,
    resolve_parameter_context,
    validate_exchange_formulas,
    add_database_exchanges_to_group,
    build_exchange_formula_map,
    evaluate_exchange_formulas,
)

# -- lca_matrices --------------------------------------------------------------
from utils.lca_matrices import (                           # noqa: F401
    aggregate_updates_by_cell,
    apply_updates_to_matrices,
    apply_grouped_updates_to_matrices,
)

# -- lca_contribution ----------------------------------------------------------
from utils.lca_contribution import (                       # noqa: F401
    contribution_analysis,
    oat_sensitivity,
    analytical_variance_from_parameter,
    parameter_metadata,
    combine_oat_and_analytical_variance,
    _reset_and_recalculate_mlca,
    _recompute_scores_with_param_context,
)

# -- lca_monte_carlo -----------------------------------------------------------
from utils.lca_monte_carlo import (                        # noqa: F401
    mc_sample,
    chunkify_rows,
    _parse_score_key,
    recalculate_scores_for_sample,
    init_worker,
    worker,
    run_parallel_monte_carlo,
)

# -- lca_sobol -----------------------------------------------------------------
from utils.lca_sobol import (                              # noqa: F401
    build_sobol_problem,
    sanitize_sobol_problem,
    generate_sobol_samples,
    expected_sobol_sample_rows,
    init_worker_sobol,
    sobol_worker,
    recalculate_scores_for_sobol_sample,
    run_parallel_sobol_from_samples,
    sobol_indices_from_results,
    sobol_indices_to_dataframe,
    run_full_sobol_workflow,
)

# -- Standard library re-exports (for compatibility) ---------------------------
from __future__ import annotations
import math                                         # noqa: F401
import multiprocessing as mp                         # noqa: F401
import os                                           # noqa: F401
import time                                         # noqa: F401
from dataclasses import dataclass                    # noqa: F401
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple  # noqa: F401
from collections import defaultdict                  # noqa: F401

import numpy as np                                  # noqa: F401
import pandas as pd                                  # noqa: F401
from scipy import stats                              # noqa: F401
import warnings                                     # noqa: F401

import bw2analyzer as bwa                            # noqa: F401
import bw2calc as bc                                 # noqa: F401
import bw2data as bd                                 # noqa: F401
from bw2data import labels                           # noqa: F401
from asteval import Interpreter                      # noqa: F401
from bw2calc.method_config import MethodConfig       # noqa: F401
from bw2calc.multi_lca import MultiLCA               # noqa: F401
from bw2data.backends.schema import ExchangeDataset  # noqa: F401
from bw2data.parameters import (                     # noqa: F401
    ActivityParameter, ParameterizedExchange, parameters,
)

# SALib re-exports
from SALib.analyze import sobol as sobol_analyze     # noqa: F401
from SALib.sample import sobol                       # noqa: F401
