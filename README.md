# Small Molecule Modification Toolkit

Python utilities for building, inspecting, and modifying small molecules as NumPy-based graphs. The toolkit connects RDKit molecules and SMILES with an editable `MolGraph` representation, allows for atom and bond-level edits, supports fragment attachment and merges, and supports functional-group level modifications. Its array-based structure supports straightforward conversion of atom, bond, and connectivity data into tensors for machine-learning applications.

## Capabilities

| Area | What it provides |
| --- | --- |
| Molecule conversion | Convert SMILES or RDKit `Mol` objects to `MolGraph`, then reconstruct RDKit molecules or canonical/isomeric SMILES. |
| Molecular features | Atom and bond arrays for element identities, formal and partial charges, aromaticity, stereochemistry, hydrogen counts, isotopes, atom maps, and optional coordinates. |
| Molecule preparation | Parsing, sanitization, Gasteiger charge calculation, and a rule-based cleanup workflow with repair logs. Cleanup can change the input structure. |
| Coordinate generation | RDKit ETKDG conformer generation, with optional MMFF optimization and single-conformer or ensemble options. |
| Graph editing | Add, delete, replace, or reorder atoms; add, delete, or edit bonds; update graph indices and invalidate derived annotations where required. |
| Fragment attachment | Query and label dummy atoms, identify their neighboring anchors, and merge fragments at labeled sites. Supports substituent attachment and insertion sites. |
| Functional groups | Detect groups from SMILES, RDKit molecules, or `MolGraph`; inspect matched atom indices, query SMARTS and group names, compare molecules, and draw highlighted structures. |
| Operation pipelines | Compose transforms, repeat operations, choose operations randomly, apply conditions, and collect operation metadata and warnings. |
| Graph utilities | Check structural invariants, derive deterministic atom orders, and compute topology or stereochemistry-aware graph hashes. |

## Project layout

```text
chem/
  build/               SMILES/RDKit conversion, molecular features, coordinates
  clean/               Molecule detection, sanitization, and repair helpers
  dummy/               Dummy-atom labels, attachment queries, and invariants
  edit/                Atom and bond editing
  functional_groups/   AccFG integration, queries, comparisons, and visualization
  merge/               Labeled fragment merging
  ops/                 Transform classes, composition, and operation logging
core/
  structs.py           MolGraph, MolArrays, and related data structures
  structs_helper.py    Array validation helpers
  molgraph/            Graph validation, ordering, hashing, and transform wrappers
integration/           Package support for integration code
utils/                 Logging utilities
run_all_tests.py        Aggregate test-suite runner
```

Test modules are located alongside their components in `tests/` or `testing/` directories. This repository contains the small-molecule toolkit.

## Molecular graph representation

[`core/structs.py`](core/structs.py) defines the shared data structures. `MolArrays` is the molecular state; the surrounding `MolGraph` holds annotations, derived caches, and provenance.

| Component | Contents |
| --- | --- |
| `MolArrays` | Required `atomic_num` and `formal_charge` vectors of length N, `bonds` with shape `(M, 2)`, and `bond_type` of length M; optional atom, bond, and coordinate arrays. Bonds represent undirected edges with normalized endpoint order. |
| `AttachmentPointArray` | Derived dummy-atom indices, labels, neighboring anchors, and insertion-site information. The canonical dummy labels live in `arrays.attachment_label`. |
| `Fragment` and `ResonanceSystem` | Atom/bond subsets and annotations for fragments or delocalized systems. |
| `Editability` | Per-atom and per-bond mutation scores supplied by the caller. |
| `DerivedCaches` and `DirtyFlags` | Directed edge indices, feature arrays, bond lookup tables, and flags indicating when cached data needs rebuilding. |
| `MolGraph.meta` | Input SMILES, repair history, operation records, merge history, and other workflow metadata. |

A single conformer's coordinates use `(N, 3)`; conformer ensembles use `(K, N, 3)`, with corresponding validity masks. The directed `edge_index` helper emits both directions of each undirected bond, giving shape `(2, 2M)` for graph-learning workflows.

Structural mutations clear affected derived caches and mark features dirty. Operations that delete or reorder atoms can return `old_to_new`, `new_to_old`, and bond-permutation mappings. Deleted atoms use `-1` in the old-to-new map. Use the returned mappings rather than reusing atom indices from before an edit; mapping availability depends on the operation.

## Setup

Use Python 3.11 and install the dependencies in a virtual environment:

```sh
git clone https://github.com/asmaniyamin/small-molecule-mod-toolkit.git
cd small-molecule-mod-toolkit
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Run the examples from the repository root so that `chem`, `core`, and `utils` are importable. The project does not currently include an installable package configuration. When calling it from another project, add this repository's root to that project's Python import path.

The requirements file pins the scientific dependencies, including RDKit, NumPy, NetworkX, and AccFG. Functional-group detection requires AccFG; its detector is loaded when first used.

## Examples

### Convert SMILES to a graph and back

```python
from chem.build.create_molgraph import smiles_to_molgraph
from chem.build.molgraph_to_mol import molgraph_to_smiles

graph = smiles_to_molgraph("CCO", compute_charges=False)
if graph is None:
    raise ValueError("Molecule preparation failed")

print(graph.arrays.atomic_num)
print(graph.arrays.bonds)
print(molgraph_to_smiles(graph))
```

`MolGraph.arrays` holds the atom and bond data. `graph.meta` carries provenance and operation metadata. Atom indices refer to the current graph and may change after deletions, permutations, or molecule preparation.

By default, conversion does not generate coordinates. Pass `coords="canonical"` for a single conformer or `coords="ensemble", num_confs=10` for multiple conformers. Check the returned coordinate shapes against the requirements of downstream editing and validation routines.

### Merge fragments directly

Call `chem.merge.merge.merge_by_labels` directly to attach a fragment to a molecule. This is the API used by [molecular-fragment-assembly](https://github.com/asmaniyamin/molecular-fragment-assembly/blob/main/env/env.py) for fragment assembly.

```python
from copy import deepcopy
from chem.build.create_molgraph import smiles_to_molgraph
from chem.build.molgraph_to_mol import molgraph_to_smiles
from chem.merge.merge import merge_by_labels
from core.structs import BondType

left = smiles_to_molgraph("[*:1]C", compute_charges=False)
right = smiles_to_molgraph("[*:2]C", compute_charges=False)
if left is None or right is None:
    raise ValueError("Fragment preparation failed")

merged = merge_by_labels(
    deepcopy(left),
    deepcopy(right),
    label_a=1,
    label_b=2,
    bond_type_code=BondType.SINGLE,
    validate_with_rdkit=True,
    on_rdkit_fail="raise",
)
print(molgraph_to_smiles(merged))  # CC
```

Here `[*:1]` and `[*:2]` are atom-mapped dummy atoms. Conversion stores their map numbers in `arrays.attachment_label`, which `label_a` and `label_b` select. The function returns the merged `MolGraph`; the selected dummy atoms are removed and their neighboring atoms are connected. The example copies its inputs because merging can modify them.

For iterative assembly, use the returned graph as the growing molecule and select one of its remaining dummy sites for the next merge. Use `chem.dummy.query.dummy_indices` and `dummy_by_label` to inspect sites, and `chem.dummy.edit.set_dummy_label` to label a selected site. When several dummies share a label, assign a unique label to the intended site before calling `merge_by_labels`, as the selected-site assembly workflow in molecular-fragment-assembly does.

Substituent attachment uses one heavy-atom neighbor at each site. Insertion uses a two-neighbor site and a one-neighbor site; merging two insertion sites is unsupported. `bond_type_code` controls the new bond type. `validate_with_rdkit=True` enables an RDKit reconstruction check, and `on_rdkit_fail="raise"` makes a failed check raise an exception. Coordinate generation can be requested through `coords`, `num_confs`, and `optimize_coords`.

Selected dummy labels are consumed, derived attachment data is invalidated, and `meta["merge_log"]` records the operation. An isotope label such as `[1*]` does not by itself populate the canonical attachment label used by this function. The optional `chem.ops.merge.MergeByLabels` wrapper exposes the same merge through the transform/logging interface; two-input merge transforms do not provide a single combined atom-index mapping.

### Apply an edit and inspect its operation log

```python
from chem.build.create_molgraph import smiles_to_molgraph
from chem.ops.atoms import ReplaceAtom
from chem.ops.base import apply_transform

graph = smiles_to_molgraph("CCO", compute_charges=False)
if graph is None:
    raise ValueError("Molecule preparation failed")

# Replace the oxygen with sulfur in this three-atom graph.
graph, info = apply_transform(graph, ReplaceAtom(idx=2, atomic_num=16))
print(graph.arrays.atomic_num)
print(info.op, info.warnings)
print(graph.meta["op_log"][-1])
```

Low-level edits are also available in `chem.edit.atoms` and `chem.edit.bonds`. Use the `chem.ops` layer when you need operation records or composition through `Compose`, `Repeat`, `Maybe`, `RandomChoice`, or `Conditional`. Edits can modify their input graph; do not assume the original graph is preserved. Structural edits do not automatically recompute all chemical properties, such as partial charges.

### Detect and visualize functional groups

```python
from pathlib import Path
from chem.functional_groups import detect_functional_groups, draw_mol_with_fgs

smiles = "CC(=O)O"
result = detect_functional_groups(smiles, include_hierarchy=True)
for match in result.matches:
    print(match.name, match.atom_indices)

output = Path("outputs")
output.mkdir(exist_ok=True)
(output / "functional_groups.png").write_bytes(draw_mol_with_fgs(smiles))
```

The result also supports named-group queries and atom-to-group mappings. `lite=True` omits detailed subcategories. Custom definitions can be supplied with `user_defined_fgs`; the hierarchy is available when requested. Matched atom indices correspond to the analyzed molecule, which may be canonicalized during detection.

### Compose operations with validation

```python
from chem.build.create_molgraph import smiles_to_molgraph
from chem.ops.atoms import ReplaceAtom
from chem.ops.base import apply_transform
from chem.ops.pipeline import Compose
from core.molgraph.ops.canonicalize import Canonicalize
from core.molgraph.ops.validate import Validate

graph = smiles_to_molgraph("CCO", compute_charges=False)
if graph is None:
    raise ValueError("Molecule preparation failed")

pipeline = Compose([
    ReplaceAtom(idx=2, atomic_num=16),
    Canonicalize(stereo=False),
    Validate(strict=True),
])
graph, info = apply_transform(graph, pipeline)
print(info.warnings)
print([entry["op"] for entry in graph.meta["op_log"]])
```

`apply_transform` returns `(graph, OpInfo)`. Records contain operation parameters, warnings, atom/bond count changes, and index mappings where available. `Compose` logs each nested operation and a pipeline summary; it does not compose the nested index mappings into one overall mapping.

`Repeat`, `Maybe`, `RandomChoice`, and `Conditional` support repeated or conditional transformations. The random combinators accept seeds. `StopOnWarning` returns the original graph reference when warnings occur; it does not undo mutations already made to that object or catch exceptions. Copy a graph before editing when you need to preserve it.

### Normalize graph order and calculate a hash

```python
from chem.build.create_molgraph import smiles_to_molgraph
from core.molgraph.canonicalize import canonicalize
from core.molgraph.hash import hash_molgraph

graph = smiles_to_molgraph("CCO", compute_charges=False)
if graph is None:
    raise ValueError("Molecule preparation failed")

graph, old_to_new, new_to_old, bond_perm = canonicalize(graph, stereo=False)
print(hash_molgraph(graph, mode="topology"))
```

`canonicalize` uses Weisfeiler-Lehman refinement to order atoms and remaps bonds and related annotations. Hashes serialize ordered atom/bond fields, so normalize ordering before comparisons. Use matching `stereo=False`/`mode="topology"` settings for topology comparisons, or `stereo=True`/`mode="stereo"` when including stereochemical fields. Hashes exclude coordinates and workflow metadata and are not a substitute for a full molecular-equivalence check.

## Validation and preparation behavior

- `core.structs_helper` checks array and attachment shapes; `core.molgraph.validate.validate_molgraph` also checks graph invariants such as endpoint bounds, self-bonds, normalized bonds, and numeric values. Strict graph validation raises on errors; non-strict validation returns messages.
- These checks are not comprehensive valence, stereochemical, geometric, or physical validation. RDKit reconstruction and sanitization provide additional checks, subject to the configured conversion options.
- `MoleculeCleaner` can retain the largest connected component and apply charge, hydrogen, tautomer, stereochemistry, kekulization, or hypervalence repairs. Inspect `meta["repair_log"]` after SMILES conversion when preserving the original input chemistry matters.
- Dummy-atom partial charges are handled using temporary hydrogen capping during Gasteiger calculations. Low-level edits do not automatically recalculate charges after changing the graph.

## Testing

Run the chemistry and core test modules from the repository root:

```sh
python -m pytest chem core
```

The transform demonstration is available separately:

```sh
python -m chem.ops.demo_ops_pipeline
```

`run_all_tests.py` will run all test suites available.

## Citation

If you use this toolkit in research, cite *Small Molecule Modification Toolkit* by Asmani Yamin. See [CITATION.cff](CITATION.cff) for software citation metadata. Acknowledge the scientific dependencies used in the reported workflow, including RDKit and AccFG where applicable.

## Research use

Computational results require independent scientific review and validation
before use in synthesis, experiments, drug discovery, or development. A
computed molecular structure or transformation does not establish compound
safety or efficacy. See [research-use limitations](RESEARCH_USE.md).

## License and attribution

Copyright 2026 Asmani Yamin.

Unmodified noncommercial use and private modifications for noncommercial
research are permitted subject to [LICENSE](LICENSE). Commercial use and
sharing modified versions require prior written permission from Asmani Yamin.
Retain [NOTICE](NOTICE) and use [CITATION.cff](CITATION.cff) when citing the work.
