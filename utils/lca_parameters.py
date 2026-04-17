"""
lca_parameters.py -- Brightway activity parameter and formula utilities

This module provides:
- get_activity_parameters(): Query activity parameters (optionally by group)
- build_parameter_formula_map(): Build parameter name -> formula mapping
- resolve_parameter_context(): Iteratively resolve parameter contexts with formulas
- validate_exchange_formulas(): Validate stored exchange amounts against formulas
- add_database_exchanges_to_group(): Register database exchanges into a parameter group
- build_exchange_formula_map(): Map parameterized exchanges to matrix coordinates
- evaluate_exchange_formulas(): Evaluate all exchange formulas with a parameter context
"""

from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

import bw2data as bd
from bw2data import labels
from bw2data.backends.schema import ExchangeDataset
from bw2data.parameters import ActivityParameter, ParameterizedExchange, parameters
from asteval import Interpreter

from utils.lca_config import _log


# -----------------------------------------------------------------------------
# Activity parameter queries
# -----------------------------------------------------------------------------

def get_activity_parameters(group_name: Optional[str] = None) -> List[ActivityParameter]:
    """Query activity parameters, optionally filtered by group."""
    query = ActivityParameter.select()
    if group_name is not None:
        query = query.where(ActivityParameter.group == group_name)
    return list(query)


def build_parameter_formula_map(group_name: str) -> List[dict]:
    """Return list of {'name': ..., 'formula': ...} for parameters that have formulas."""
    rows = get_activity_parameters(group_name)
    return [
        {"name": row.name, "formula": getattr(row, "formula", None)}
        for row in rows
        if getattr(row, "formula", None)
    ]


# -----------------------------------------------------------------------------
# Parameter context resolution (iterative formula evaluation)
# -----------------------------------------------------------------------------

def resolve_parameter_context(
    group_name: str,
    sampled_values_by_name: Mapping[str, float],
    parameter_formula_map: Optional[Sequence[dict]] = None,
    max_iter: int = 50,
    tol: float = 1e-15,
) -> Dict[str, float]:
    """
    Resolve a consistent parameter context including derived formulas.

    This is safer than evaluating exchange formulas against raw parameter amounts only.
    """
    context = ActivityParameter.static(group_name, full=True).copy()
    context.update(dict(sampled_values_by_name))
    formulas = list(parameter_formula_map or build_parameter_formula_map(group_name))

    if not formulas:
        return context

    interpreter = Interpreter()
    interpreter.symtable.update(context)

    for _ in range(max_iter):
        changed = 0
        for item in formulas:
            value = interpreter(item["formula"])
            if interpreter.error:
                errors = "; ".join(err.get_error()[1] for err in interpreter.error)
                raise ValueError(f"Parameter formula eval failed for {item['name']}: {errors}")
            interpreter.error = []
            new = float(value)
            old = context.get(item["name"])
            if old is None or abs(new - old) > tol:
                context[item["name"]] = new
                interpreter.symtable[item["name"]] = new
                changed += 1
        if changed == 0:
            return context

    raise RuntimeError(
        f"Parameter formulas for group '{group_name}' did not converge after {max_iter} iterations"
    )


# -----------------------------------------------------------------------------
# Exchange formula validation and registration
# -----------------------------------------------------------------------------

def validate_exchange_formulas(database_name: str, group_name: str, tol: float = 1e-6) -> pd.DataFrame:
    """Validate stored exchange amounts against formulas using resolved parameter context."""
    db = bd.Database(database_name)
    param_context = resolve_parameter_context(group_name, sampled_values_by_name={})
    interpreter = Interpreter()
    interpreter.symtable.update(param_context)

    rows = []
    for act in db:
        for exc in act.exchanges():
            formula = exc.get("formula")
            if not formula:
                continue
            stored_amount = exc["amount"]
            try:
                calculated = float(interpreter(formula))
                if interpreter.error:
                    errors = "; ".join(err.get_error()[1] for err in interpreter.error)
                    raise ValueError(errors)
                interpreter.error = []
                diff = calculated - stored_amount
                status = "OK" if abs(diff) <= tol else "MISMATCH"
            except Exception as exc_err:
                calculated = np.nan
                diff = np.nan
                status = f"ERROR: {exc_err}"

            rows.append(
                {
                    "Activity name": act.get("name", ""),
                    "Activity code": act["code"],
                    "Exchange input": str(exc.input),
                    "Stored amount": stored_amount,
                    "Formula": formula,
                    "Calculated from formula": calculated,
                    "Difference": diff,
                    "Status": status,
                }
            )
    return pd.DataFrame(rows)


def add_database_exchanges_to_group(database_name: str, group_name: str) -> None:
    """Register all exchanges of a database into a parameter group and recalculate them."""
    for act in bd.Database(database_name):
        parameters.add_exchanges_to_group(group_name, act)
    ActivityParameter.recalculate_exchanges(group_name)


# -----------------------------------------------------------------------------
# Exchange → matrix coordinate mapping
# -----------------------------------------------------------------------------

def build_exchange_formula_map(group_name: str, reference_mlca) -> List[dict]:
    """
    Map parameterized exchanges to matrix coordinates using Brightway label semantics.

    Returns list of dicts with keys:
        exchange_id, formula, matrix, exchange_type, row_index, col_index, sign
    """
    exchange_formula_map = []

    pos_types = set(getattr(labels, "technosphere_positive_edge_types", []))
    neg_types = set(getattr(labels, "technosphere_negative_edge_types", []))
    bio_types = set(getattr(labels, "biosphere_edge_types", ["biosphere"]))

    for pe in ParameterizedExchange.select().where(ParameterizedExchange.group == group_name):
        exc_ds = ExchangeDataset.get(id=pe.exchange)
        exc = exc_ds.data

        input_node = bd.get_node(database=exc["input"][0], code=exc["input"][1])
        output_node = bd.get_node(database=exc["output"][0], code=exc["output"][1])
        exc_type = exc.get("type")

        if exc_type in bio_types:
            matrix = "biosphere_matrix"
            row_index = reference_mlca.dicts.biosphere.get(input_node.id)
            sign = 1.0

        elif exc_type in pos_types or exc_type in neg_types:
            matrix = "technosphere_matrix"
            row_index = reference_mlca.dicts.product.get(input_node.id)
            sign = 1.0 if exc_type in pos_types else -1.0

        else:
            continue

        col_index = reference_mlca.dicts.activity.get(output_node.id)

        if row_index is None or col_index is None:
            continue

        exchange_formula_map.append(
            {
                "exchange_id": pe.exchange,
                "formula": pe.formula,
                "matrix": matrix,
                "exchange_type": exc_type,
                "row_index": int(row_index),
                "col_index": int(col_index),
                "sign": float(sign),
            }
        )

    return exchange_formula_map


def evaluate_exchange_formulas(
    exchange_formula_map: Sequence[dict],
    param_context: Mapping[str, float],
) -> List[dict]:
    """
    Evaluate all exchange formulas with the given parameter context.

    Each returned dict extends the input with a 'matrix_value' key.
    """
    interpreter = Interpreter()
    interpreter.symtable.update(dict(param_context))

    updates = []
    for item in exchange_formula_map:
        value = interpreter(item["formula"])
        if interpreter.error:
            errors = "; ".join(err.get_error()[1] for err in interpreter.error)
            raise ValueError(f"Formula eval failed for exchange {item['exchange_id']}: {errors}")
        interpreter.error = []
        updates.append({**item, "matrix_value": float(value) * item["sign"]})
    return updates
