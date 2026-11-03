# Copyright 2026 Asmani Yamin
# SPDX-License-Identifier: LicenseRef-Asmani-Yamin-Noncommercial-Research-1.0

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple

from core.structs import MolGraph
from chem.ops.base import Transform, OpInfo, ChangeSummary
from core.molgraph.validate import validate_molgraph, MolGraphValidationError


@dataclass
class Validate(Transform):
    strict: bool = True
    name: str = "validate"

    def apply(self, mg: MolGraph) -> Tuple[MolGraph, OpInfo]:
        """Check graph invariants using the configured strictness.

        Args:
            mg: MolGraph to process using this transform's configured fields.

        Returns:
            A (graph, OpInfo) pair. The graph structure is unchanged. Non-strict results are recorded in
            OpInfo.warnings; strict validation fails instead of returning a result.

        Call apply_transform(graph, op) to append the returned operation record to
        the graph log; calling this method directly does not log this outer operation.

        Raises:
            MolGraphValidationError: If strict validation detects an invalid graph.
        """
        n0 = int(mg.arrays.atomic_num.shape[0])
        m0 = int(mg.arrays.bonds.shape[0])

        warnings = []
        if self.strict:
            # will raise MolGraphValidationError if invalid
            validate_molgraph(mg, strict=True)
        else:
            warnings = validate_molgraph(mg, strict=False)

        info = OpInfo(
            op=self.name,
            params={"strict": bool(self.strict)},
            warnings=list(warnings),
            changes=ChangeSummary(
                n_atoms_before=n0,
                n_atoms_after=n0,
                n_bonds_before=m0,
                n_bonds_after=m0,
                notes=["validation only (no structural changes)"],
            ),
        )
        return mg, info
