# meshtools3d

C++ tool for generating 3D tetrahedral meshes from segmented medical images (.inr format),
built on CGAL 4.x. Targets cardiac mesh generation (atria, ventricles) with region labeling,
boundary extraction, and wall-thickness evaluation via a Laplace solver.

## Build

Requires: CGAL (4.6+), Boost, optional Intel TBB for parallelism.

```sh
cmake -DCGAL_DIR=<path_to_CGAL> .
make
```

Debug build (disables rounding-math checks required for valgrind):
```sh
cmake -DCMAKE_BUILD_TYPE=Debug -DCGAL_DONT_OVERRIDE_CMAKE_FLAGS=FALSE -DCGAL_DISABLE_ROUNDING_MATH_CHECK=ON .
```

To disable TBB: comment out the `if( TBB_FOUND )` block in `CMakeLists.txt`.

## Run

```sh
cp data-template data
# edit data with paths and meshing parameters
./meshtools3d -f data
```

Key CLI overrides (supersede values in the data file):
```
-seg_dir / -seg_name      input segmentation directory and filename
-out_dir / -out_name      output directory and filename prefix
--read_the_mesh / -meshr  read an existing mesh instead of generating one
--thickness-algorithm N   N=1 (Martin Bishop) or N=2 (Cesare Corrado)
-v                        verbose output
```

Parallelism: set `TBB_NUM_THREADS` env var before running. The Docker entrypoint
(`dockerM3D.sh`) handles this automatically; pass `-a` to auto-detect cores.

## Architecture

- `applications/main.cpp` — orchestrator: reads the data file via GetPot, sets up CGAL
  mesh criteria, drives the pipeline, writes outputs. All file I/O lives here.
- `m3d/include/` + `m3d/src/` — library classes:
  - `Mesh` — internal mesh representation (points, triangles, tetrahedra with region labels);
    handles CARP output, boundary extraction, rescaling, VTX lists.
  - `INRreader` — reads `.inr` segmentation volumes; supports bounding-box label queries
    for the manual-segmentation (re-labeling) workflow.
  - `LaplaceSolver` — GMRES-based solver for the harmonic extension used in thickness.
  - `ThicknessEvaluation` — extends `LaplaceSolver`; two algorithms (Bishop/Corrado).
  - `VtkWriter` — writes legacy VTK unstructured grid files (ASCII or binary).
  - `CGALDataType.hpp` — all CGAL typedef aliases; selects parallel or sequential
    triangulation types based on `CGAL_LINKED_WITH_TBB`.
  - `GetPot` — third-party header-only config/CLI parser.

## Two Meshing Modes

1. **`mesh_from_segmentation = 1`** (default): CGAL reads the `.inr` directly via
   `Labeled_image_mesh_domain_3`. Region labels come from the image.
2. **`mesh_from_segmentation = 0`** (manual / re-labeling mode): `INRreader` loads the
   segmentation; meshing uses `My_Labeled_image_mesh_domain_3` with no labels (all tetras
   get label=1); triangle labels are assigned post-hoc by centroid lookup against the INR
   bounding boxes. Use this when CGAL's direct label meshing produces artifacts.

## Domain Terminology

- **INR / `.inr`**: INRIA image format — the required input segmentation format.
- **CARP**: cardiac simulation file format (`.elem`, `.pts`, `.vtx`) used by the openCARP
  and CARP solvers.
- **VTX files**: point-index lists identifying anatomical surface regions (endocardium,
  epicardium, mitral valves, etc.).
- **Endo/Epi**: endocardium (inner wall) and epicardium (outer wall) — determined
  automatically as the two largest boundary regions.
- **rescaleFactor**: applied at output time for CARP/VTK; does not affect `.mesh` output.
  CARP expects coordinates in micrometers (factor ~1000 from mm).

## Gotchas

- `NDEBUG` controls whether array accesses use `.at()` (bounds-checked) or `[]` (unchecked).
  Release builds silently skip bounds checks; use Debug builds when debugging indexing issues.
- `boundary_relabeling = 0` skips the endo/epi detection step entirely — VTX files will
  not be written.
- Binary CARP output (`out_carp_binary = 1`) is implemented but untested.
- `thicknessOld/` is legacy code — do not modify or build it.
- `matlab_utilities/` provides MATLAB helpers for converting NRRD segmentations to INR.
