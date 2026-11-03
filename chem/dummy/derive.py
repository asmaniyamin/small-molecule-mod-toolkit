# Copyright 2026 Asmani Yamin
# SPDX-License-Identifier: LicenseRef-Asmani-Yamin-Noncommercial-Research-1.0

# chem/dummy/derive.py
from __future__ import annotations
from typing import Optional
import numpy as np
from core.structs import MolGraph, AttachmentPointArray, AttachmentKind
from chem.dummy.query import dummy_indices, dummy_target, dummy_labels, neighbors

def build_attachment_points(g: MolGraph) -> Optional[AttachmentPointArray]:
    """Build a derived attachment-point view from canonical dummy atoms.

    Args:
        g: MolGraph with dummy atoms identified by atomic_num == 0 and optional
            labels stored in arrays.attachment_label.

    Returns:
        An AttachmentPointArray, or None if no dummies exist. Entries follow dummy
        atom order and use int32 atom indices. A target is the sole non-dummy
        neighbor, or -1 when ambiguous or unavailable. is_insertion marks dummies
        with at least two neighbors. label_id is None if all labels are absent;
        insertion_anchors is left unset.

    The graph is not changed and the result is not assigned to g.attachments.
    """
    d = dummy_indices(g)
    if d.size == 0:
        return None

    idx = d.astype(np.int32, copy=False)
    kind = np.full(d.size, AttachmentKind.DUMMY, dtype=np.int8)

    # target: -1 if ambiguous/unavailable
    tgt = np.empty(d.size, dtype=np.int32)
    for i, di in enumerate(d):
        t = dummy_target(g, int(di))
        tgt[i] = int(t) if t is not None else -1

    labels = None
    if g.arrays.attachment_label is not None:
        lab = dummy_labels(g)
        labels = lab if np.any(lab != None) else None  # noqa: E711

    # insertion flags (optional)
    is_ins = np.array([neighbors(g, int(di)).size >= 2 for di in d], dtype=bool)

    return AttachmentPointArray(
        idx=idx,
        kind=kind,
        target=tgt,
        label_id=labels,
        is_insertion=is_ins,
        insertion_anchors=None,
    )
