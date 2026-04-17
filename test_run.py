"""
test_run.py -- Quick validation: project setup + deterministic LCA

Before running, update the configuration section below with your own values.
"""
import sys
import os
import warnings

warnings.filterwarnings('ignore')

import traceback
import bw2data as bd
from utils import lca_project_setup


# ============================================================
# Configuration -- UPDATE THESE VALUES BEFORE RUNNING
# ============================================================
PROJECT_NAME = 'your_project_name'
FOREGROUND_DB = 'your_foreground_db'
ECOINVENT_DB_PREFIX = 'ecoinvent'  # Used to filter impact methods


def run_step(name, fn):
    print(f"\n{'=' * 60}")
    print(f"[STEP] {name}")
    print(f"{'=' * 60}", flush=True)
    try:
        fn()
        print(f"[OK] {name} completed successfully", flush=True)
        return True
    except Exception as e:
        print(f"[FAIL] {name}: {e}", flush=True)
        traceback.print_exc()
        return False


def step_setup():
    lca_project_setup.set_project(PROJECT_NAME)
    print(f"  Project: {bd.projects.current}")
    print(f"  Databases: {list(bd.databases)}")


def step_methods():
    all_methods = list(bd.methods)
    methods = []
    for m in all_methods:
        if 'gwp100' in str(m[-1]).lower() and m[0].startswith(ECOINVENT_DB_PREFIX):
            methods.append(m)
    methods.sort()
    print(f"  Found {len(methods)} GWP100 methods")
    for m in methods[:3]:
        print(f"    {m}")
    print(f"    ... and {len(methods) - 3} more")
    globals()['METHODS'] = methods


def step_fu():
    db = bd.Database(FOREGROUND_DB)
    fu = [db.get(code=a['code']) for a in db]
    print(f"  Found {len(fu)} activities")
    globals()['FU'] = fu


def step_deterministic_lca():
    methods = globals()['METHODS']
    fu = globals()['FU']
    print(f"  Running LCA with {len(fu)} FUs and {len(methods)} methods...")
    df = lca_project_setup.deterministic_lca(fu, methods)
    print(f"  Results shape: {df.shape}")
    top = df.reindex(df.Score.abs().sort_values(ascending=False).index).head(5)
    print("  Top 5 scores:")
    for _, row in top.iterrows():
        print(f"    {row['Functional unit'][:35]:35s} | {row['Score']:.4e}")


if __name__ == '__main__':
    results = []
    results.append(run_step("Setup project", step_setup))
    if results[-1]:
        results.append(run_step("Get GWP100 methods", step_methods))
    if results[-1]:
        results.append(run_step("Select functional units", step_fu))
    if results[-1]:
        results.append(run_step("Deterministic LCA", step_deterministic_lca))

    print(f"\n{'=' * 60}")
    print(f"SUMMARY: {sum(results)}/{len(results)} steps passed")
    if all(results):
        print("All sections passed!")
    else:
        print("Some steps failed.")
        sys.exit(1)
