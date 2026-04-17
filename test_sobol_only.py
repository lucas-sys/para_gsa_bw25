"""
test_sobol_only.py -- Test only the Sobol GSA step

Before running, update the configuration section below with your own values.
"""
import sys
import os
import datetime
import warnings

warnings.filterwarnings('ignore')

_log_fh = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_sobol_result.txt'),
               'w', encoding='utf-8')


def log(msg):
    ts = datetime.datetime.now().strftime('%H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    _log_fh.write(line + '\n')
    _log_fh.flush()


log("=== Sobol-only test started ===")
log(f"Python: {sys.executable}")

import traceback
import multiprocessing as mp

log("Importing bw2data...")
import bw2data as bd

log("Importing lca modules...")
from utils import lca_project_setup
from utils import lca_parameters
from utils import lca_sobol
import pandas as pd

log("All imports done")


# ============================================================
# Configuration -- UPDATE THESE VALUES BEFORE RUNNING
# ============================================================
PROJECT_NAME = 'your_project_name'
FOREGROUND_DB = 'your_foreground_db'
GROUP_NAME = 'your_parameter_group'
FU_CODE = 'your_fu_code'

METHODS = [
    # Update with your impact methods, e.g.:
    # ('your_ecoinvent_db', 'EF v3.1', 'climate change',
    #  'global warming potential (GWP100)'),
]


if __name__ == '__main__':
    try:
        log("Setting up project...")
        lca_project_setup.set_project(PROJECT_NAME)
        db = bd.Database(FOREGROUND_DB)
        fu = db.get(code=FU_CODE)
        FUS = [fu]
        log(f"FU: {fu['name']}")

        log("Building Sobol problem...")
        problem, params_used = lca_sobol.build_sobol_problem(GROUP_NAME)
        log(f"Sobol problem: {len(problem['names'])} parameters")
        log(f"Distributions: {sorted(set(problem['dists']))}")

        N = 16
        log(f"Generating Sobol samples (N={N})...")
        samples = lca_sobol.generate_sobol_samples(problem, N=N)
        log(f"Samples shape: {samples.shape}")

        log("Running parallel Sobol (n_workers=2)...")
        results_df = lca_sobol.run_parallel_sobol_from_samples(
            project_name=PROJECT_NAME,
            group_name=GROUP_NAME,
            functional_units=FUS,
            methods=METHODS,
            sobol_problem=problem,
            sobol_samples=samples,
            n_workers=2,
        )
        log(f"Sobol results shape: {results_df.shape}")

        raw_pkl = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'results/sobol/sobol_raw_results.pkl')
        results_df.to_pickle(raw_pkl)
        log(f"Saved raw results to {raw_pkl}")

        log("Computing Sobol indices...")
        all_dfs = []
        for fu in FUS:
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
                    log(f"  [OK] '{fu_label}' x '{method[2]}'")
                except Exception as exc:
                    log(f"  [WARN] '{fu_label}' x '{method[2]}': {exc}")

        if all_dfs:
            final_df = pd.concat(all_dfs, ignore_index=True)
            log(f"Final Sobol DataFrame: {final_df.shape}")
            out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'results/sobol/sobol_results.xlsx')
            final_df.to_excel(out, index=False)
            log(f"Saved to {out}")
            log("Top 10 results:")
            _log_fh.write(final_df.head(10).to_string() + '\n')
            _log_fh.flush()
            log("=== Sobol test PASSED ===")
        else:
            log("=== Sobol test FAILED: no indices computed ===")

    except Exception as e:
        log(f"=== Sobol test FAILED: {e} ===")
        traceback.print_exc()
        _log_fh.write(traceback.format_exc())
        _log_fh.flush()
    finally:
        _log_fh.close()
