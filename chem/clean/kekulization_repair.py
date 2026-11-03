# Copyright 2026 Asmani Yamin
# SPDX-License-Identifier: LicenseRef-Asmani-Yamin-Noncommercial-Research-1.0

import traceback
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import ResonanceMolSupplier
from chem.clean.fatal_repair_tools import clear_nonring_aromaticity

# --- Logging wrapper that passes kwargs ------------------------------------------------
def step_with_logging(step_fn, smiles, verbose, **kwargs):
    """Run one repair strategy, logging its outcome when requested.

    Args:
        step_fn: Callable accepting a SMILES string and keyword arguments.
        smiles: Input SMILES passed to the strategy.
        verbose: Print success messages or failure tracebacks when true.
        **kwargs: Keyword arguments forwarded unchanged to step_fn.

    Returns:
        The strategy's molecule when truthy, otherwise None. Exceptions derived
        from Exception are caught and converted to None so later steps can run.
    """
    try:
        mol = step_fn(smiles, **kwargs)
        if mol:
            if verbose:
                print(f"[PASS] Success in {step_fn.__name__}")
            return mol
    except Exception as e:
        if verbose:
            print(f"[FAIL] Failed in {step_fn.__name__}:\n{traceback.format_exc()}")
    return None

# --- Step 1: Sanitize all except kekulize, then attempt Kekulize -----------------------
def try_sanitize_then_kekulize_preserving_aromatic(smiles, **kwargs):
    """Sanitize a SMILES molecule before attempting kekulization.

    Parsing skips sanitization initially. All sanitization operations except
    kekulization are then applied, non-ring aromatic flags are cleared, and
    kekulization runs without clearing the remaining aromatic flags.

    Args:
        smiles: Input SMILES string.
        **kwargs: Unused options accepted for compatibility with the strategy runner.

    Returns:
        The newly parsed or repaired RDKit molecule.

    Parsing or repair exceptions propagate to the caller. Use step_with_logging
    to convert a strategy failure into None.
    """
    params = Chem.SmilesParserParams()
    params.sanitize = False
    mol = Chem.MolFromSmiles(smiles, params)
    flags = Chem.SanitizeFlags.SANITIZE_ALL & ~Chem.SanitizeFlags.SANITIZE_KEKULIZE
    Chem.SanitizeMol(mol, flags)
    clear_nonring_aromaticity(mol)
    Chem.Kekulize(mol, clearAromaticFlags=False)
    return mol

# --- Step 2: Sanitize skipping aromaticity perception (aromatizeIfPossible=False) ------
def try_sanitize_without_aromatization(smiles, **kwargs):
    """Attempt sanitization with the requested aromaticity option.

    The aromatizeIfPossible keyword is forwarded to the RDKit sanitization call.
    Compatibility with that keyword depends on the installed RDKit API.

    Args:
        smiles: Input SMILES string.
        **kwargs: Additional strategy options. aromatizeIfPossible is read; other keys
            are ignored.

    Returns:
        The newly parsed or repaired RDKit molecule.

    Parsing or repair exceptions propagate to the caller. Use step_with_logging
    to convert a strategy failure into None.
    """
    aromatize = kwargs.get("aromatizeIfPossible", False)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    flags = Chem.SanitizeFlags.SANITIZE_ALL
    Chem.rdmolops.SanitizeMol(mol, flags, aromatizeIfPossible=aromatize)
    return mol

# --- Step 3: Clear aromaticity flags and sanitize --------------------------------------
def try_clear_aromaticity_and_sanitize(smiles, **kwargs):
    """Clear all aromatic flags before sanitizing a SMILES molecule.

    Both atom and bond aromatic flags are cleared on the newly parsed molecule.
    This fallback can alter the aromatic representation.

    Args:
        smiles: Input SMILES string.
        **kwargs: Unused options accepted for compatibility with the strategy runner.

    Returns:
        The newly parsed or repaired RDKit molecule.

    Parsing or repair exceptions propagate to the caller. Use step_with_logging
    to convert a strategy failure into None.
    """
    params = Chem.SmilesParserParams()
    params.sanitize = False
    mol = Chem.MolFromSmiles(smiles, params)
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        bond.SetIsAromatic(False)
    Chem.SanitizeMol(mol)
    return mol

# --- Step 4: Add explicit hydrogens ----------------------------------------------------
def try_add_hs_and_sanitize(smiles, **kwargs):
    """Add explicit hydrogens before sanitizing a SMILES molecule.

    The molecule is parsed without sanitization before hydrogen addition.

    Args:
        smiles: Input SMILES string.
        **kwargs: Unused options accepted for compatibility with the strategy runner.

    Returns:
        The newly parsed or repaired RDKit molecule.

    Parsing or repair exceptions propagate to the caller. Use step_with_logging
    to convert a strategy failure into None.
    """
    params = Chem.SmilesParserParams()
    params.sanitize = False
    mol = Chem.MolFromSmiles(smiles, params)
    mol = Chem.AddHs(mol)
    Chem.SanitizeMol(mol)
    return mol

# --- Step 5: Normalize charges / tautomers ---------------------------------------------
def try_normalize_and_sanitize(smiles, **kwargs):
    """Normalize a SMILES molecule before sanitizing it.

    Uses RDKit MolStandardize.Normalizer on an unsanitized input molecule.

    Args:
        smiles: Input SMILES string.
        **kwargs: Unused options accepted for compatibility with the strategy runner.

    Returns:
        The newly parsed or repaired RDKit molecule.

    Parsing or repair exceptions propagate to the caller. Use step_with_logging
    to convert a strategy failure into None.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)
    Chem.SanitizeMol(mol)
    return mol

# --- Step 6: Try resonance forms and pick a good one -----------------------------------
def try_resonance_forms(smiles, **kwargs):
    """Return the first resonance form that passes RDKit sanitization.

    Initial SMILES parsing uses normal sanitization. Individual resonance-form
    sanitization failures are skipped; initial parsing and supplier errors propagate.

    Args:
        smiles: Input SMILES string.
        **kwargs: Unused options accepted for compatibility with the strategy runner.

    Returns:
        An RDKit molecule for the first successful form, or None if none pass.

    Parsing or repair exceptions propagate to the caller. Use step_with_logging
    to convert a strategy failure into None.
    """
    mol = Chem.MolFromSmiles(smiles)
    resonance_forms = ResonanceMolSupplier(mol, Chem.KEKULE_ALL)
    for res_mol in resonance_forms:
        try:
            Chem.SanitizeMol(res_mol)
            return res_mol
        except:
            continue
    return None

# --- Main controller -------------------------------------------------------------------
def fix_kekulization_modular(smiles, verbose=True, extreme=False, aromatizeIfPossible=False):
    """Try successive RDKit strategies to repair a SMILES kekulization failure.

    Args:
        smiles: SMILES string reparsed independently by each strategy.
        verbose: Print strategy successes and failure tracebacks.
        extreme: Append a final attempt that clears all aromatic flags.
        aromatizeIfPossible: Forwarded to the sanitization strategy that reads
            this option; other strategies ignore it.

    Returns:
        The first successfully repaired RDKit molecule. Strategies attempt
        sanitization and kekulization, alternate aromaticity handling, hydrogen
        addition, normalization, and resonance forms in that order.

    Raises:
        ValueError: If all strategies fail or return no molecule.

    Individual strategy exceptions are caught by step_with_logging. Repairs may
    change the molecular representation and should be checked by the caller.
    """
    steps = [
        try_sanitize_then_kekulize_preserving_aromatic,
        try_sanitize_without_aromatization,
        try_add_hs_and_sanitize,
        try_normalize_and_sanitize,
        try_resonance_forms,
    ]
    
    if extreme:
        steps.append(try_clear_aromaticity_and_sanitize)

    for step in steps:
        mol = step_with_logging(step, smiles, verbose, aromatizeIfPossible=aromatizeIfPossible)
        if mol:
            return mol

    raise ValueError(f"[FAIL] Kekulization failed for SMILES: {smiles}")
