# Copyright 2026 Asmani Yamin
# SPDX-License-Identifier: LicenseRef-Asmani-Yamin-Noncommercial-Research-1.0

# chem/dummy/edit.py
from __future__ import annotations
from typing import Any
import numpy as np
from core.structs import MolGraph
from chem.dummy.query import dummy_by_label

def _ensure_attachment_label_array(g: MolGraph) -> None:
    if g.arrays.attachment_label is None:
        n = int(g.arrays.atomic_num.shape[0])
        g.arrays.attachment_label = np.empty(n, dtype=object)
        g.arrays.attachment_label[:] = None

def set_dummy_label(g: MolGraph, idx: int, label: Any) -> MolGraph:
    """Assign a canonical attachment label to a dummy atom in place.

    Args:
        g: MolGraph to update. A missing attachment_label array is initialized
            as a length-N object array containing None.
        idx: Atom index, converted to int and interpreted using NumPy indexing.
        label: Value to store; uniqueness is not checked.

    Returns:
        The same graph object. Derived attachment views are not rebuilt.

    Raises:
        ValueError: If the selected atom is not a dummy (atomic_num != 0).
        IndexError: If the index falls outside the NumPy array bounds.
    """
    idx = int(idx)
    if int(g.arrays.atomic_num[idx]) != 0:
        raise ValueError(f"Atom {idx} is not a dummy (atomic_num != 0)")
    _ensure_attachment_label_array(g)
    g.arrays.attachment_label[idx] = label
    return g

def clear_dummy_label(g: MolGraph, idx: int) -> MolGraph:
    """Clear an atom's attachment label in place.

    Args:
        g: MolGraph to update. If no attachment_label array exists, nothing changes.
        idx: Atom index, converted to int and interpreted using NumPy indexing.
            This helper does not require the atom to be a dummy.

    Returns:
        The same graph object, with the selected label set to None when present.
        Derived attachment views are not rebuilt.

    Raises:
        IndexError: If the label array exists and the index is outside its bounds.
    """
    idx = int(idx)
    if g.arrays.attachment_label is None:
        return g
    g.arrays.attachment_label[idx] = None
    return g

def relabel_dummy(g: MolGraph, old_label: Any, new_label: Any) -> MolGraph:
    """Replace an existing dummy label with an unused label in place.

    Args:
        g: MolGraph containing the labeled dummy atom.
        old_label: Existing label resolved through dummy_by_label, which compares
            labels as strings and selects the first matching dummy.
        new_label: Replacement label, which must not already identify a dummy.

    Returns:
        The same graph object with the canonical label updated.

    Raises:
        ValueError: If old_label is absent or new_label already exists, including
            when old_label and new_label are the same. Label-resolution errors
            from dummy_by_label also propagate.
    """
    idx = dummy_by_label(g, old_label)
    if idx is None:
        raise ValueError(f"No dummy found with label '{old_label}'")
    # optional: prevent duplicates among existing dummy labels
    if dummy_by_label(g, new_label) is not None:
        raise ValueError(f"Label '{new_label}' already exists")
    return set_dummy_label(g, idx, new_label)

def enforce_label_invariant(g: MolGraph) -> MolGraph:
    """
    Enforce: non-dummy atoms must have attachment_label = None.
    This is a canonical graph invariant.
    """
    if g.arrays.attachment_label is None:
        return g
    non_dummy = g.arrays.atomic_num != 0
    g.arrays.attachment_label[non_dummy] = None
    return g
