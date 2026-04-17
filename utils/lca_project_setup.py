"""
lca_project_setup.py -- Brightway project initialization, MultiLCA construction, and deterministic LCA

This module provides:
- set_project(): Switch the active Brightway project
- build_fu_labels(): Generate stable unique labels for functional units
- build_demands(): Build MultiLCA demand dictionaries
- initialize_multilca(): Create LCASetup with LCI/LCIA computed and base matrices cached
- flatten_multilca_scores(): Flatten mlca.scores to a uniform key format
- deterministic_lca(): Run deterministic LCA and return a tidy DataFrame
"""

from __future__ import annotations

from typing import Dict, Mapping, Sequence, Tuple

import pandas as pd

import bw2data as bd
from bw2calc.method_config import MethodConfig
from bw2calc.multi_lca import MultiLCA

from utils.lca_config import LCASetup


def set_project(project_name: str) -> None:
    """Set Brightway project and force database metadata load."""
    bd.projects.set_current(project_name)
    _ = bd.databases


def build_fu_labels(functional_units: Sequence) -> Dict[object, str]:
    """Build stable, unique labels for FUs even when names repeat."""
    return {fu: f"{fu.get('name', 'FU')}__{i}" for i, fu in enumerate(functional_units)}


def build_demands(fu_labels: Mapping[object, str]) -> Dict[str, Dict[int, float]]:
    return {label: {fu.id: 1.0} for fu, label in fu_labels.items()}


def initialize_multilca(functional_units: Sequence, methods: Sequence[tuple]) -> LCASetup:
    """Create a fully initialized LCASetup with LCI/LCIA computed and base matrices cached."""
    fu_labels = build_fu_labels(functional_units)
    demands = build_demands(fu_labels)

    method_config = MethodConfig(impact_categories=list(methods))
    method_config_dict = {
        key: value for key, value in method_config.model_dump().items() if value is not None
    }
    data_objs = bd.get_multilca_data_objs(
        functional_units=demands,
        method_config=method_config_dict,
    )

    mlca = MultiLCA(
        demands=demands,
        method_config=method_config_dict,
        data_objs=data_objs,
    )
    mlca.lci()
    mlca.lcia()

    return LCASetup(
        project_name=bd.projects.current,
        group_name="",
        functional_units=functional_units,
        methods=list(methods),
        fu_labels=fu_labels,
        demands=demands,
        method_config_dict=method_config_dict,
        data_objs=data_objs,
        mlca=mlca,
        base_tech=mlca.technosphere_matrix.copy().tocsr(),
        base_bio=mlca.biosphere_matrix.copy().tocsr(),
    )


def flatten_multilca_scores(scores_dict: Mapping) -> Dict[Tuple[str, tuple], float]:
    """Normalize mlca.scores to {(fu_label, method_tuple): score}."""
    flat = {}
    for key, score in scores_dict.items():
        fu_label = None
        method = None
        if isinstance(key, tuple) and len(key) == 2:
            a, b = key
            if isinstance(a, tuple) and isinstance(b, str):
                method, fu_label = a, b
            elif isinstance(a, str) and isinstance(b, tuple):
                fu_label, method = a, b
        if fu_label is None:
            fu_label = "unknown_fu"
            method = key if isinstance(key, tuple) else (str(key),)
        flat[(fu_label, method)] = float(score)
    return flat


def deterministic_lca(functional_units: Sequence, methods: Sequence[tuple]) -> pd.DataFrame:
    """Run deterministic LCA and return a tidy DataFrame of scores."""
    setup = initialize_multilca(functional_units, methods)
    rows = []
    for (fu_label, method), score in flatten_multilca_scores(setup.mlca.scores).items():
        rows.append(
            {
                "Functional unit": fu_label,
                "Method": method[2] if isinstance(method, tuple) and len(method) > 2 else str(method),
                "Method full": method,
                "Score": float(score),
            }
        )
    return pd.DataFrame(rows).sort_values(["Functional unit", "Method"]).reset_index(drop=True)
