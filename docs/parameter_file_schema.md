# meshtools3d Parameter File Schema

Reference for the parameter (`.par`) file consumed by `meshtools3d` via its
`-f <data_file>` argument.

The **source of truth** for this schema is the bundled `parfile_builder`
application (`applications/parfile_builder.cpp`), which emits a complete file
with the defaults documented here. If you only need a starting file, run:

```bash
parfile_builder -o heart.par
```

This document is for anyone writing their own generator (Python, shell, Matlab,
GUI, etc.) or hand-editing a file without taking a dependency on
`parfile_builder`.

---

## 1. File format

The file is an INI-style block format read by GetPot: `[section]` headers group
`key = value` entries, and `meshtools3d` looks up values as `section/key`.

- One `[section]` header per group of parameters.
- `key = value` entries beneath each section.
- Keys are **case-sensitive** — `rescaleFactor` and `dimKrilovSp` must stay
  camelCase.
- All values are written as strings; numerics are parsed by the C++ side.
- Section order does not matter, but key names must match exactly.
- `meshtools3d` does **not** fill in missing keys — emit all five sections with
  every key.

A complete valid file (identical to `parfile_builder` output with defaults):

```ini
[segmentation]
seg_dir = ./
seg_name = seg_final_smooth_corrected.inr
mesh_from_segmentation = 1
boundary_relabeling = 0

[meshing]
facet_angle = 30
facet_size = 0.8
facet_distance = 4
cell_rad_edge_ratio = 2.0
cell_size = 0.8
rescaleFactor = 1000

[laplacesolver]
abs_toll = 1e-6
rel_toll = 1e-6
itr_max = 700
dimKrilovSp = 500
verbose = 1

[others]
eval_thickness = 0

[output]
outdir = ./myocardium_OUT
name = heart_mesh
out_medit = 0
out_carp = 1
out_carp_binary = 0
out_vtk = 1
out_vtk_binary = 0
out_potential = 0
```

### CLI overrides

A few keys can be overridden on the `meshtools3d` command line; these take
precedence over the data file (useful inside scripts):
`-seg_dir`, `-seg_name`, `-mesh_dir`, `-mesh_name`, `-out_dir`, `-out_name`.

---

## 2. Sections and keys

All keys below are **required**. Defaults match `parfile_builder`'s built-in
schema.

### 2.1 `[segmentation]`

Inputs describing the source segmentation.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `seg_dir` | path | `./` | Directory containing the segmentation file. |
| `seg_name` | filename | `seg_final_smooth_corrected.inr` | INR segmentation file name. |
| `mesh_from_segmentation` | bool (0/1) | `1` | If `1`, build the mesh from the segmentation. |
| `boundary_relabeling` | bool (0/1) | `0` | Apply boundary relabeling pass. |

### 2.2 `[meshing]`

CGAL meshing criteria. Tune `facet_size` and `cell_size` for resolution.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `facet_angle` | float (deg) | `30` | Lower bound on triangle angles. |
| `facet_size` | float | `0.8` | Upper bound on facet circumradius. **Resolution knob.** |
| `facet_distance` | float | `4` | Max distance facet ↔ surface. |
| `cell_rad_edge_ratio` | float | `2.0` | Tetra radius/edge ratio bound. |
| `cell_size` | float | `0.8` | Upper bound on tetra circumradius. **Resolution knob.** |
| `rescaleFactor` | float | `1000` | Scaling applied to CARP / VTK output (note camelCase). |

### 2.3 `[laplacesolver]`

Linear solver settings used by the Laplace pass (e.g. thickness / potential).
Shared with the standalone `laplace_solver` application.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `abs_toll` | float | `1e-6` | Absolute tolerance (note spelling: `toll`). |
| `rel_toll` | float | `1e-6` | Relative tolerance. |
| `itr_max` | int | `700` | Max iterations. |
| `dimKrilovSp` | int | `500` | Krylov subspace dimension (note camelCase). |
| `verbose` | bool (0/1) | `1` | Solver verbosity. |

### 2.4 `[others]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `eval_thickness` | bool (0/1) | `0` | If `1`, run wall-thickness evaluation. Also controls VTK contents: `1` writes potential+thickness, `0` writes region_labels. |

### 2.5 `[output]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `outdir` | path | `./myocardium_OUT` | Output directory (created if missing). |
| `name` | string | `heart_mesh` | Base name for output files. |
| `out_medit` | bool (0/1) | `0` | Write MEDIT `.mesh`. |
| `out_carp` | bool (0/1) | `1` | Write CARP ASCII (`.elem`, `.pts`, `.lon`). |
| `out_carp_binary` | bool (0/1) | `0` | Write CARP binary (untested against a real openCARP consumer — see Roadmap). |
| `out_vtk` | bool (0/1) | `1` | Write VTK ASCII. |
| `out_vtk_binary` | bool (0/1) | `0` | Write VTK binary. |
| `out_potential` | bool (0/1) | `0` | Write Laplace potential field. |

---

## 3. Generating files with `parfile_builder`

`parfile_builder` is stdlib-only (no CGAL, no `m3d/` link) and bakes in the
defaults above. Output goes to stdout unless `-o` is given.

```
Usage: parfile_builder [OPTIONS]

  -o, --out PATH            Output file path
  --seg_dir DIR             Set segmentation.seg_dir
  --seg_name FILE           Set segmentation.seg_name
  --out_dir DIR             Set output.outdir
  --out_name NAME           Set output.name
  --set SECTION.KEY=VAL     Override any key (repeatable)
  -h, --help                Show this help
```

Examples:

```bash
parfile_builder > heart.par
parfile_builder -o heart.par --seg_dir /data --out_dir /data/out
parfile_builder --set meshing.facet_size=0.5 --set output.out_vtk=1
```

`--set SECTION.KEY=VALUE` validates both the section and the key against the
schema and exits with an error on an unknown name.

---

## 4. Reference Python writer (no dependencies)

For pipelines that prefer to emit the file from Python, this self-contained
writer reproduces the same schema. No third-party imports.

```python
import configparser

DEFAULTS = {
    "segmentation": {
        "seg_dir": "./",
        "seg_name": "seg_final_smooth_corrected.inr",
        "mesh_from_segmentation": "1",
        "boundary_relabeling": "0",
    },
    "meshing": {
        "facet_angle": "30",
        "facet_size": "0.8",
        "facet_distance": "4",
        "cell_rad_edge_ratio": "2.0",
        "cell_size": "0.8",
        "rescaleFactor": "1000",
    },
    "laplacesolver": {
        "abs_toll": "1e-6",
        "rel_toll": "1e-6",
        "itr_max": "700",
        "dimKrilovSp": "500",
        "verbose": "1",
    },
    "others": {
        "eval_thickness": "0",
    },
    "output": {
        "outdir": "./myocardium_OUT",
        "name": "heart_mesh",
        "out_medit": "0",
        "out_carp": "1",
        "out_carp_binary": "0",
        "out_vtk": "1",
        "out_vtk_binary": "0",
        "out_potential": "0",
    },
}

def write_par(path: str, overrides: dict | None = None) -> None:
    cfg = configparser.ConfigParser()
    cfg.optionxform = str  # preserve case (rescaleFactor, dimKrilovSp)
    cfg.read_dict(DEFAULTS)
    for section, kvs in (overrides or {}).items():
        if section not in cfg:
            raise KeyError(f"unknown section: {section}")
        for k, v in kvs.items():
            if k not in cfg[section]:
                raise KeyError(f"unknown key: [{section}] {k}")
            cfg[section][k] = str(v)
    with open(path, "w") as fh:
        cfg.write(fh)
```

Usage:

```python
write_par(
    "heart.par",
    overrides={
        "segmentation": {"seg_dir": "/data/case01", "seg_name": "seg.inr"},
        "meshing":      {"facet_size": 0.5, "cell_size": 0.5},
        "output":       {"outdir": "/data/case01/mesh", "name": "case01"},
    },
)
```

---

## 5. Authoring checklist

When writing your own generator:

1. Emit **all five sections** with **every key** from §2 — `meshtools3d` does
   not fill in missing keys.
2. Preserve case for `rescaleFactor` and `dimKrilovSp`.
3. Keep tolerance keys spelled `abs_toll` / `rel_toll` (the existing files use
   the `toll` spelling).
4. Booleans are integers `0` / `1`, not `true` / `false`.
5. Paths can be relative; `meshtools3d` resolves them against its own CWD, so
   prefer absolute paths in batch pipelines.
6. Write values as strings; do not quote them.
7. After writing, round-trip with an INI parser (or diff against
   `parfile_builder` output) to confirm the file parses.
