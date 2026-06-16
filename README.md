# MeshTools3D

MeshTools3D generates 3D tetrahedral meshes from `.inr` image segmentations
using CGAL 6.x, with optional Laplace harmonic extension and wall-thickness
computation. It is built for cardiac volumetric models but works on any
labelled `.inr` segmentation.

The build produces three command-line tools:

- **`meshtools3d`** — generate a tetrahedral mesh from a segmentation.
- **`laplace_solver`** — run the harmonic-extension + thickness pipeline on an
  existing CARP mesh (no CGAL required).
- **`parfile_builder`** — generate a parameter file with sensible defaults.

Outputs: 
+ CARP (`.elem` / `.pts` + per-region `.vtx`), 
+ VTK, and 
+ INRIA `.mesh` (medit)

See [`docs/parameter_file_schema.md`](docs/parameter_file_schema.md) for the full parameter reference.

---

## Install (pre-built binaries)

Download the archive for your platform from the
[Releases page](https://github.com/OpenHeartDevelopers/meshtools3d/releases),
then extract and run:

```bash
tar xzf meshtools3d-<version>-<os>-<arch>.tar.gz
./meshtools3d-<version>-<os>-<arch>/bin/meshtools3d --help
```

Each release also publishes a `.sha256` sidecar you can verify with
`sha256sum -c <file>.sha256`.

If **macOS** blocks an archive downloaded via a browser, try [these instructions](docs/macOS_bug_fix_v2.0-beta.md)


---

## Quick start

**1. Create a parameter file.** Use `parfile_builder` (recommended):

```bash
parfile_builder -o heart.par --seg_dir /data/case01 --seg_name seg.inr --out_dir /data/case01/mesh
```

Override any individual key with `--set SECTION.KEY=VALUE`, e.g.
`--set meshing.facet_size=0.5`. The complete schema (sections, keys, defaults)
is documented in
[`docs/parameter_file_schema.md`](docs/parameter_file_schema.md).

**2. Generate the mesh:**

```bash
meshtools3d -f heart.par
```

The `-seg_dir`, `-seg_name`, `-out_dir`, and `-out_name` flags override the
matching values in the parameter file (handy in scripts). Run
`meshtools3d --help` for the full list, including `--read_the_mesh` (mesh input
instead of segmentation) and `--thickness-algorithm`.

**3. (Optional) Run the standalone Laplace solver** on an existing CARP mesh:

```bash
laplace_solver -mesh_dir /data/case01/mesh -mesh_name case01 -out_dir out -out_name case01 --vtk
```

Each `meshtools3d` run also writes `<out_name>_params.data` and
`<out_name>_invocation.sh` next to the outputs, so any run can be reproduced
exactly.

---

## Parallel runs (TBB)

If Intel TBB is available, meshing is parallelised. The thread count defaults to
**1** (so the tool doesn't grab every core) and is controlled by the
`TBB_NUM_THREADS` environment variable:

```bash
export TBB_NUM_THREADS=4
```

> `TBB_NUM_THREADS` is specific to MeshTools3D — it is not a standard TBB
> variable.

---

## Build from source

Quick start (Linux / macOS, CGAL 6.x via system packages):

```bash
# Dependencies (macOS):  brew install cmake cgal boost tbb
# Dependencies (Linux):  apt install cmake libcgal-dev libboost-dev libtbb-dev zlib1g-dev libeigen3-dev

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The three tools land in `build/`. Full dependency details, version requirements,
and platform notes are in [`docs/install.md`](docs/install.md). To reproduce the
legacy CGAL 4.x environment for v1.0 results, see
[`docs/install_legacy.md`](docs/install_legacy.md).

---

## Documentation

- [`docs/install.md`](docs/install.md) — full build instructions (Linux / macOS).
- [`docs/install_legacy.md`](docs/install_legacy.md) — legacy CGAL 4.x build for
  reproducing v1.0 results.
- [`docs/parameter_file_schema.md`](docs/parameter_file_schema.md) — parameter
  file reference.
- [`docs/Roadmap.md`](docs/Roadmap.md) — planned work and known gaps.
- `examples/` — sample inputs (`sphereCoarse`, `sphereCoarse4Thickness`,
  `sphereMultilabel`, …).

A Python wrapper, **`pycemrg-meshing`**, is maintained separately — it authors
parameter files, fetches released binaries, and runs the tools from Python.

---

## Limitations

These are tracked in [`docs/Roadmap.md`](docs/Roadmap.md):

- **Windows is best-effort** — the release build may fail (GetPot /
  `boost::filesystem` issue); Linux and macOS are the supported platforms.
- CARP **binary** output (`out_carp_binary = 1`) is implemented but untested
  against a real openCARP consumer.
- `.mesh` re-emission from an already-loaded `Mesh` object, boundary-triangle
  surface export, and outward-normal triangle re-orientation are not yet wired
  up for general use.

---

## Changelog

Per-release notes are on the
[Releases page](https://github.com/OpenHeartDevelopers/meshtools3d/releases).

## Contributors

- Cesare Corrado — original author
- Jose Alonso Solis-Lemus
