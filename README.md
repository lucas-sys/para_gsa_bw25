## Contributors

<a href="https://github.com/lucas-sys/para_gsa_bw25/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=lucas-sys/para_gsa_bw25" />
</a>

# LCA Global Sensitivity Analysis (GSA) Toolkit

A modular Python toolkit for performing **parameter-perturbation Life Cycle Assessment (LCA)** with **Global Sensitivity Analysis (GSA)**, built on [Brightway2.5](https://brightway.dev/), [SALib](https://salib.readthedocs.io/), and Python multiprocessing.

This toolkit was developed based on previous work of Rebecca Belfiore, Christhel Andrade Diaz, and Weier Liu. Code X was used to improve the code efficiency, and WorkBuddy was used to organize all files. 

## Features

- **Deterministic LCA** -- Compute baseline environmental impact scores across multiple functional units and impact categories
- **Contribution Analysis** -- Recursive Tier-2 contribution analysis to identify key process contributions
- **Formula Validation** -- Verify stored exchange amounts against parameterized formulas
- **OAT Sensitivity** -- One-at-a-time finite-difference sensitivity analysis with elasticity coefficients
- **Analytical Variance** -- Variance decomposition from uncertainty metadata (13 distribution types)
- **Monte Carlo Simulation** -- Parallel stochastic sampling with configurable worker count
- **Sobol GSA** -- First-order (S1) and total-order (ST) Sobol indices via SALib with parallel evaluation

## Project Structure

```
.
├── utils/                            # Core LCA modules
│   ├── __init__.py                   # Public API re-exports
│   ├── lca_config.py                 # Configuration, logging, Timer, LCASetup
│   ├── lca_project_setup.py          # Brightway project init, MultiLCA
│   ├── lca_parameters.py             # Parameter queries, formula resolution
│   ├── lca_matrices.py               # Sparse matrix update utilities
│   ├── lca_contribution.py           # Contribution analysis, OAT, analytical variance
│   ├── lca_monte_carlo.py            # Monte Carlo sampling & parallel evaluation
│   ├── lca_sobol.py                  # Sobol GSA (build, sample, evaluate, analyze)
│   └── lca_gsa_helper.py             # Backward-compatible re-export layer
├── data/                             # Input databases (user-provided)
│   └── ...                           # Place your Excel database files here
├── results/                          # Output files (organized by analysis type)
│   ├── deterministic/
│   ├── contribution/
│   ├── sensitivity_oat/
│   ├── sensitivity_analytical/
│   ├── monte_carlo/
│   └── sobol/
├── gsa_workflow.ipynb                # Interactive notebook (full workflow)
├── init_databases.py                 # Database initialization script
├── test_full_workflow.py             # End-to-end test
├── test_run.py                       # Quick validation test
├── test_sobol_only.py                # Sobol-only test
├── requirements.txt
├── LICENSE
└── README.md
```

## Prerequisites

- Python >= 3.10
- [Brightway2.5](https://brightway.dev/) (`bw2data`, `bw2calc`, `bw2io`, `bw2analyzer`)
- An ecoinvent database (requires ecoinvent license)

## Installation

```bash
# Create and activate conda environment
conda create -n bw25 python=3.10
conda activate bw25

# Install core dependencies
pip install -r requirements.txt
```

## Configuration

Before running any analysis, you need to configure the following parameters:

| Parameter                     | Description                                    | Where to set                                 |
| ----------------------------- | ---------------------------------------------- | -------------------------------------------- |
| `PROJECT_NAME`                | Your Brightway project name                    | `init_databases.py`, test scripts, notebook  |
| `FOREGROUND_DB`               | Name of your foreground database               | Test scripts, notebook                       |
| `GROUP_NAME`                  | ActivityParameter group name                   | Test scripts, notebook                       |
| `FU_CODE`                     | Activity code for the functional unit          | Test scripts, notebook                       |
| `METHODS`                     | List of impact method tuples                   | Test scripts, notebook                       |
| `ECOINVENT_USERNAME/PASSWORD` | ecoinvent credentials                          | `init_databases.py` or environment variables |
| `CUSTOM_DATABASES`            | List of (db_name, xlsx_path, match_dbs) tuples | `init_databases.py`                          |

All parameters use `your_*` placeholders in the template files. Replace them with your actual values before running.

## Database Initialization

```bash
# Set ecoinvent credentials (recommended: use environment variables)
export ECOINVENT_USERNAME="your_username"
export ECOINVENT_PASSWORD="your_password"

# Or edit init_databases.py directly and update the configuration section

# Run the initialization script
python init_databases.py
```

Place your custom database Excel files in the `data/` directory and add entries to the `CUSTOM_DATABASES` list in `init_databases.py`.

## Quick Start

### In Python

```python
from utils import (
    set_project, initialize_multilca, deterministic_lca,
    oat_sensitivity, run_parallel_monte_carlo,
    build_sobol_problem, generate_sobol_samples,
    run_parallel_sobol_from_samples, sobol_indices_from_results,
)

import bw2data as bd

# Setup
PROJECT_NAME = 'your_project_name'
FOREGROUND_DB = 'your_foreground_db'
GROUP_NAME = 'your_parameter_group'

set_project(PROJECT_NAME)
db = bd.Database(FOREGROUND_DB)
fu = [db.get(code='your_fu_code')]

methods = [
    ('your_ecoinvent_db', 'EF v3.1', 'climate change',
     'global warming potential (GWP100)'),
]

# Deterministic LCA
df = deterministic_lca(fu, methods)

# Monte Carlo (1000 iterations, 2 workers)
mc_results = run_parallel_monte_carlo(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    n_samples=1000,
    n_workers=2,
)

# Sobol GSA
problem, params = build_sobol_problem(GROUP_NAME)
samples = generate_sobol_samples(problem, N=512)
sobol_results = run_parallel_sobol_from_samples(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    sobol_problem=problem,
    sobol_samples=samples,
    n_workers=2,
)
```

### In Jupyter Notebook

Open `gsa_workflow.ipynb` and run cells sequentially. The notebook covers the complete workflow from database setup through Sobol GSA. Before running, update the configuration cells with your project-specific values.

## Workflow Sections

| Section | Description               | Output Directory                  |
| ------- | ------------------------- | --------------------------------- |
| 0       | Brightway project setup   | --                                |
| 1       | Impact method selection   | --                                |
| 2       | Functional unit selection | --                                |
| 3       | Deterministic LCA         | `results/deterministic/`          |
| 4       | Contribution analysis     | `results/contribution/`           |
| 5       | OAT sensitivity           | `results/sensitivity_oat/`        |
| 6       | Analytical variance       | `results/sensitivity_analytical/` |
| 7       | Monte Carlo simulation    | `results/monte_carlo/`            |
| 8       | Sobol GSA                 | `results/sobol/`                  |

## Supported Uncertainty Distributions

| Code | Distribution     | Parameters                               |
| ---- | ---------------- | ---------------------------------------- |
| 0, 1 | No uncertainty   | --                                       |
| 2    | Lognormal        | loc (ln_mean), scale (ln_std)            |
| 3    | Normal           | loc (mean), scale (std)                  |
| 4    | Uniform          | minimum, maximum                         |
| 5    | Triangular       | minimum, maximum, loc (mode)             |
| 6    | Bernoulli        | loc (probability)                        |
| 7    | Discrete uniform | minimum, maximum                         |
| 8    | Weibull          | shape, loc, scale                        |
| 9    | Gamma            | shape (alpha), loc, scale                |
| 10   | Beta             | shape (alpha), shape2 (beta), loc, scale |
| 11   | Gumbel (GEV)     | shape (xi), loc, scale                   |
| 12   | Student-t        | shape (df), loc, scale                   |

## Multiprocessing

Both Monte Carlo and Sobol evaluations use Python's `spawn` multiprocessing context for Windows compatibility. Each worker independently loads the Brightway project and databases, so memory usage scales with worker count. Recommended: 2-4 workers.

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

- [Brightway2.5](https://brightway.dev/) -- LCA framework
- [SALib](https://salib.readthedocs.io/) -- Sensitivity Analysis Library
- [ecoinvent](https://ecoinvent.org/) -- Life cycle inventory database
