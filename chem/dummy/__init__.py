# Copyright 2026 Asmani Yamin
# SPDX-License-Identifier: LicenseRef-Asmani-Yamin-Noncommercial-Research-1.0

# chem/dummy/__init__.py
"""Manage dummy atoms used as molecular attachment points.

A dummy atom has arrays.atomic_num equal to zero. Its canonical label lives in
arrays.attachment_label; attachment views and metadata label maps are derived
and can be rebuilt. The query, edit, derive, and clean modules provide lookup,
label mutation, view construction, and validation respectively.
"""

from .query import (
    dummy_indices,
    dummy_labels,
    dummy_by_label,
    neighbors,
    dummy_target,
    is_insertion_dummy,
)

from .edit import (
    set_dummy_label,
    clear_dummy_label,
    relabel_dummy,
    enforce_label_invariant,
)

# derived view builder
try:
    from .derive import build_attachment_points  # noqa: F401
except Exception:  # pragma: no cover
    build_attachment_points = None  # type: ignore

# cleaning/validation helpers
try:
    from .clean import (
        validate_dummy_invariants,  # noqa: F401
        remove_orphan_dummies_mark_only,  # noqa: F401
    )
except Exception:  # pragma: no cover
    validate_dummy_invariants = None  # type: ignore
    remove_orphan_dummies_mark_only = None  # type: ignore


__all__ = [
    # query
    "dummy_indices",
    "dummy_labels",
    "dummy_by_label",
    "neighbors",
    "dummy_target",
    "is_insertion_dummy",
    # edit
    "set_dummy_label",
    "clear_dummy_label",
    "relabel_dummy",
    "enforce_label_invariant",
    # derive
    "build_attachment_points",
    # clean/validate
    "validate_dummy_invariants",
    "remove_orphan_dummies_mark_only",
]
