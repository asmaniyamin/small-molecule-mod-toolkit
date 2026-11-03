# Copyright 2026 Asmani Yamin
# SPDX-License-Identifier: LicenseRef-Asmani-Yamin-Noncommercial-Research-1.0

# chem/functional_groups/__init__.py
"""Expose AccFG functional-group detection for MolGraph structures.

The package provides group lookup, atom-index mappings, hierarchical
relationships, and visualization. Detection accepts SMILES or MolGraph input;
lite mode omits detailed subcategories.
"""

from chem.functional_groups.match import (
    FunctionalGroupMatch,
    FunctionalGroupResult,
)

from chem.functional_groups.detect import (
    # Core detection
    detect_functional_groups,
    detect_from_molgraph,
    detect_from_mol,
    compare_molecules,
    # Library inspection
    get_accfg_instance,
    get_smarts_for_fg,
    list_all_functional_groups,
    search_functional_groups,
    # Hierarchy visualization
    print_fg_tree,
    # Molecule visualization
    draw_mol_with_fgs,
    draw_compare_mols,
    # Convenience functions
    find_carboxylic_acids,
    find_amines,
    find_halides,
    find_alcohols,
)

__all__ = [
    # Core classes
    "FunctionalGroupMatch",
    "FunctionalGroupResult",
    # Detection functions
    "detect_functional_groups",
    "detect_from_molgraph",
    "detect_from_mol",
    "compare_molecules",
    # Library inspection
    "get_accfg_instance",
    "get_smarts_for_fg",
    "list_all_functional_groups",
    "search_functional_groups",
    # Visualization
    "print_fg_tree",
    "draw_mol_with_fgs",
    "draw_compare_mols",
    # Convenience functions
    "find_carboxylic_acids",
    "find_amines",
    "find_halides",
    "find_alcohols",
]