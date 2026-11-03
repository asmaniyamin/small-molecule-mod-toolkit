# Copyright 2026 Asmani Yamin
# SPDX-License-Identifier: LicenseRef-Asmani-Yamin-Noncommercial-Research-1.0

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple, Optional

from core.structs import MolGraph
from chem.ops.base import Transform, OpInfo, ChangeSummary
from core.molgraph.hash import hash_molgraph


@dataclass
class HashMolGraph(Transform):
    mode: str = "stereo"          # "stereo" or "topology"
    write_to_meta: bool = True
    meta_key: Optional[str] = None
    name: str = "hash_molgraph"

    def apply(self, mg: MolGraph) -> Tuple[MolGraph, OpInfo]:
        """Hash an already canonicalized graph using the configured mode.

        Args:
            mg: MolGraph to process using this transform's configured fields.

        Returns:
            A (graph, OpInfo) pair. The same graph is returned. OpInfo.params["hash"] contains the digest.
            If write_to_meta is true, the digest is also stored under meta_key or
            hash_<mode>. No canonicalization is performed.

        Call apply_transform(graph, op) to append the returned operation record to
        the graph log; calling this method directly does not log this outer operation.

        Raises:
            ValueError: If mode is unsupported.
        """
        n0 = int(mg.arrays.atomic_num.shape[0])
        m0 = int(mg.arrays.bonds.shape[0])

        key = hash_molgraph(mg, mode=str(self.mode))

        if self.write_to_meta:
            mg.meta = dict(getattr(mg, "meta", {}) or {})
            k = self.meta_key or f"hash_{self.mode}"
            mg.meta[k] = key

        info = OpInfo(
            op=self.name,
            params={"mode": str(self.mode), "hash": key, "write_to_meta": bool(self.write_to_meta)},
            changes=ChangeSummary(
                n_atoms_before=n0,
                n_atoms_after=n0,
                n_bonds_before=m0,
                n_bonds_after=m0,
                notes=["hash computed (no structural changes)"],
            ),
        )
        return mg, info
