"""
utils — LCA Global Sensitivity Analysis (GSA) Toolkit

Modular toolkit for parameter-perturbation Life Cycle Assessment
using Brightway2.5, SALib, and multiprocessing.

Sub-modules
-----------
lca_config          Global configuration, logging, timing, LCASetup dataclass
lca_project_setup   Brightway project initialization, MultiLCA construction
lca_parameters      Activity parameter queries, formula resolution
lca_matrices        Sparse matrix update utilities
lca_contribution    Contribution analysis, OAT sensitivity, analytical variance
lca_monte_carlo     Monte Carlo sampling and parallel evaluation
lca_sobol           Sobol global sensitivity analysis
lca_fast            FAST / RBD-FAST global sensitivity analysis
lca_gsa_helper      Backward-compatible re-export layer
"""

from utils.lca_config import (                         # noqa: F401
    VERBOSE,
    LCASetup,
    _log,
    Timer,
)

from utils.lca_project_setup import (                   # noqa: F401
    set_project,
    build_fu_labels,
    build_demands,
    initialize_multilca,
    flatten_multilca_scores,
    deterministic_lca,
)

from utils.lca_parameters import (                      # noqa: F401
    get_activity_parameters,
    build_parameter_formula_map,
    resolve_parameter_context,
    validate_exchange_formulas,
    add_database_exchanges_to_group,
    build_exchange_formula_map,
    evaluate_exchange_formulas,
)

from utils.lca_matrices import (                        # noqa: F401
    aggregate_updates_by_cell,
    apply_updates_to_matrices,
    apply_grouped_updates_to_matrices,
)

from utils.lca_contribution import (                    # noqa: F401
    contribution_analysis,
    oat_sensitivity,
    analytical_variance_from_parameter,
    parameter_metadata,
    combine_oat_and_analytical_variance,
)

from utils.lca_monte_carlo import (                     # noqa: F401
    mc_sample,
    chunkify_rows,
    _parse_score_key,
    recalculate_scores_for_sample,
    init_worker,
    worker,
    run_parallel_monte_carlo,
)

from utils.lca_sobol import (                           # noqa: F401
    build_sobol_problem,
    sanitize_sobol_problem,
    generate_sobol_samples,
    expected_sobol_sample_rows,
    run_parallel_sobol_from_samples,
    sobol_indices_from_results,
    sobol_indices_to_dataframe,
    run_full_sobol_workflow,
)

from utils.lca_fast import (                            # noqa: F401
    build_fast_problem,
    sanitize_fast_problem,
    generate_fast_samples,
    expected_fast_sample_rows,
    run_parallel_fast_from_samples,
    fast_indices_from_results,
    fast_indices_to_dataframe,
    run_full_fast_workflow,
)

from utils.lca_gsa_helper import *                      # noqa: F401,F403
