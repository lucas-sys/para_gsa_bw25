"""
test_full_workflow.py -- Run the full GSA workflow (Sections 3-9)

Before running, update the configuration section below with your own values.
"""
import sys
import os
import warnings

warnings.filterwarnings('ignore')

import traceback
import multiprocessing as mp

import bw2data as bd
from utils import lca_project_setup
from utils import lca_parameters
from utils import lca_contribution
from utils import lca_monte_carlo
from utils import lca_sobol


# ============================================================
# Configuration -- UPDATE THESE VALUES BEFORE RUNNING
# ============================================================
PROJECT_NAME = 'your_project_name'
FOREGROUND_DB = 'your_foreground_db'
GROUP_NAME = 'your_parameter_group'
FU_CODE = 'your_fu_code'  # Activity code for the functional unit

METHODS = [
    # Update with your impact methods, e.g.:
    # ('your_ecoinvent_db', 'EF v3.1', 'climate change',
    #  'global warming potential (GWP100)'),
]


def run_step(name, fn):
    print(f"\n{'=' * 60}", flush=True)
    print(f"[STEP] {name}", flush=True)
    print(f"{'=' * 60}", flush=True)
    try:
        fn()
        print(f"[OK] {name}", flush=True)
        return True
    except Exception as e:
        print(f"[FAIL] {name}: {e}", flush=True)
        traceback.print_exc()
        return False


# ============================================================
# Steps
# ============================================================
def step_setup():
    lca_project_setup.set_project(PROJECT_NAME)
    print(f"  Project: {bd.projects.current}")
    print(f"  Databases: {list(bd.databases)}")
    db = bd.Database(FOREGROUND_DB)
    fu = db.get(code=FU_CODE)
    print(f"  FU: {fu['name']} (code={fu['code'][:20]}...)")
    globals()['FUS'] = [fu]


def step_register_exchanges():
    lca_parameters.add_database_exchanges_to_group(FOREGROUND_DB, GROUP_NAME)
    df = lca_parameters.validate_exchange_formulas(FOREGROUND_DB, GROUP_NAME)
    print(f"  Formula check: {df.shape[0]} exchanges, {df.shape[1]} columns")
    errors = df[df.get('valid', True) == False] if 'valid' in df.columns else df[
        df.iloc[:, -1].astype(str).str.contains('ERROR', na=False)]
    if len(errors) > 0:
        print(f"  WARNING: {len(errors)} exchanges have formula errors")
    else:
        print(f"  All formulas valid")
    df.to_excel('results/deterministic/formula_check_results.xlsx', index=False)


def step_deterministic():
    df = lca_project_setup.deterministic_lca(globals()['FUS'], METHODS)
    print(f"  Results: {df.shape}")
    print(df.to_string())
    df.to_excel('results/deterministic/deterministic_results.xlsx', index=False)
    globals()['deterministic_df'] = df


def step_contribution():
    df = lca_contribution.contribution_analysis(globals()['FUS'], METHODS, max_level=1,
                                                cutoff=0.001)
    print(f"  Contribution results: {df.shape}")
    df.to_excel('results/contribution/contribution_analysis_results.xlsx', index=False)
    print(df.head(10).to_string())


def step_oat():
    df = lca_contribution.oat_sensitivity(globals()['FUS'], METHODS, GROUP_NAME, rel_step=0.1)
    print(f"  OAT results: {df.shape}")
    df.to_excel('results/sensitivity_oat/oat_sensitivity_results.xlsx', index=False)
    print(df.head(5).to_string())
    globals()['sensitivity_df'] = df


def step_analytical():
    param_meta = lca_contribution.parameter_metadata(GROUP_NAME)
    print(f"  Parameter metadata: {param_meta.shape}")
    df = lca_contribution.combine_oat_and_analytical_variance(globals()['sensitivity_df'],
                                                              param_meta)
    print(f"  Analytical results: {df.shape}")
    df.to_excel('results/sensitivity_analytical/analytical_sensitivity_results.xlsx',
                index=False)
    print(df.head(5).to_string())


def step_monte_carlo():
    df = lca_monte_carlo.run_parallel_monte_carlo(
        project_name=PROJECT_NAME,
        group_name=GROUP_NAME,
        functional_units=globals()['FUS'],
        methods=METHODS,
        n_samples=100,
        n_workers=min(6, mp.cpu_count()),
        random_seed=42,
        export_parameter_samples='results/monte_carlo/parameter_samples.xlsx',
    )
    print(f"  MC results: {df.shape}")
    df.to_excel('results/monte_carlo/monte_carlo_results.xlsx', index=False)
    print(df.head(5).to_string())


def step_sobol():
    problem, params_used = lca_sobol.build_sobol_problem(GROUP_NAME)
    print(f"  Sobol problem: {len(problem['names'])} parameters")

    samples = lca_sobol.generate_sobol_samples(problem, N=16)
    print(f"  Samples shape: {samples.shape}")

    results_df = lca_sobol.run_parallel_sobol_from_samples(
        project_name=PROJECT_NAME,
        group_name=GROUP_NAME,
        functional_units=globals()['FUS'],
        methods=METHODS,
        sobol_problem=problem,
        sobol_samples=samples,
        n_workers=2,
    )
    print(f"  Sobol results shape: {results_df.shape}")

    all_dfs = []
    for fu in globals()['FUS']:
        fu_label = fu['name']
        for method in METHODS:
            try:
                Si = lca_sobol.sobol_indices_from_results(
                    sobol_problem=problem,
                    sobol_results_df=results_df,
                    functional_unit=fu_label,
                    method=method,
                )
                df_si = lca_sobol.sobol_indices_to_dataframe(Si, problem)
                df_si['Functional unit'] = fu_label
                df_si['Method'] = method[2] if len(method) > 2 else str(method)
                all_dfs.append(df_si)
                print(f"  [OK] Sobol indices for '{fu_label}' x '{method[2]}'")
            except Exception as exc:
                print(f"  [WARN] Sobol indices failed for '{fu_label}' x '{method[2]}': {exc}")

    if all_dfs:
        import pandas as pd
        final_df = pd.concat(all_dfs, ignore_index=True)
        print(f"  Final Sobol DataFrame: {final_df.shape}")
        print(final_df.head(10).to_string())
        final_df.to_excel('results/sobol/sobol_results.xlsx', index=False)
    else:
        print("  WARNING: No Sobol indices computed")


if __name__ == '__main__':
    steps = [
        ("Setup & FU selection", step_setup),
        ("Register exchanges", step_register_exchanges),
        ("Deterministic LCA", step_deterministic),
        ("Contribution analysis", step_contribution),
        ("OAT sensitivity", step_oat),
        ("Analytical variance", step_analytical),
        ("Monte Carlo", step_monte_carlo),
        ("Sobol GSA", step_sobol),
    ]

    results = []
    for name, fn in steps:
        r = run_step(name, fn)
        results.append(r)
        if not r:
            print(f"\nStopping after failure in: {name}")
            break

    print(f"\n{'=' * 60}")
    print(f"FINAL SUMMARY: {sum(results)}/{len(results)} steps passed")
    if all(results):
        print("ALL TESTS PASSED!")
    else:
        print("Some steps failed.")
        sys.exit(1)
