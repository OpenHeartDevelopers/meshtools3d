# MeshTools3D Roadmap

This document tracks planned refactors and feature work. Phases are the structured
architectural work the library is moving through; the feature backlog at the end
captures smaller items carried over from the original README and to be picked up
once the phases are complete.

---

## Completed

### Phase 1 + 2a — Laplace pipeline extraction + standalone solver
- Extracted Laplace setup/solve/output from `main.cpp` into `m3d/LaplacePipeline`
- New `applications/laplace_solver.cpp` standalone app (no CGAL dependency)
- `LaplaceBCConfig` accepts both VTX file paths and in-memory node sets
- VTK output split: `eval_thickness=true` writes potential+thickness; `eval_thickness=false`
  writes region_labels. Region labels always written as `.vtx` files regardless.
- **Status:** merged to `master`

### Phase 4 — CGAL 4.x to 6.x upgrade
- Modernised `CMakeLists.txt` to target-based linking (`CGAL::CGAL`, `TBB::tbb`)
- Updated `CGALDataType.hpp` and the custom domain wrappers for the CGAL 5/6 API
- TBB migrated from `task_scheduler_init` to `tbb::global_control`
- Dockerfiles updated; legacy CGAL 4.14.3 from-source install preserved in
  `Dockerfile.legacy` and `docs/install_legacy.md`
- **Status:** merged to `master` via `feature/cgal_upgrade`

### Phase 2b — Meshing-side refactor of `main.cpp`
- `MeshingParams` struct extracted, loaded via `loadMeshingParams()`
- TBB setup isolated into `configureTbbThreads()` helper
- Post-meshing boilerplate templated in `m3d/include/MeshingPipeline.hpp`
  (`writeMeditFile`, `validateTriangulation`, `populateCarpMeshFromC3t3`)
- Labeled vs. manual-segmentation branches collapsed via traits bundles
  (`LabeledMeshingTraits`, `ManualSegMeshingTraits`) and a single
  `runCGALMeshing<Traits>()` template
- Run reproducibility: `snapshotRunInputs()` writes `_params.data` + `_invocation.sh`
- **Status:** on `development`, not yet merged to `master`

### Phase 5a — Release CI + packaging
- `.github/workflows/release.yml`: Linux/macOS/Windows build matrix, CGAL 6.1.1,
  triggered on `release: published`
- `CMakeLists.txt`: CPack, install rpath, Homebrew dylib bundling (macOS),
  vcpkg DLL bundling (Windows), system-lib deps documented for Linux
- **Status:** on `development`, not yet exercised (no tag cut yet)

---

## Remaining Phases

### Phase 2c — Move `MeshingParams` out of CGAL
> To be considered as an option. Not decided!
- Change numeric fields from `FaceNumericalType`/`CellNumericalType` to `double`
- Move the struct from `CGALDataType.hpp` to `m3d/include/` alongside `LaplaceParams`
- `main.cpp` passes `double` values to CGAL criteria constructors directly (they accept double)

The only risk is if a future CGAL version makes the criteria constructors finicky about
`double` vs. their own `FT` — unlikely but worth noting.

### Phase 3 — Python wrapper
Replace the user's repeated ad-hoc parameter-file scripts with a proper Python API.

- Subprocess runner that invokes `meshtools3d` / `laplace_solver`
- GetPot data-file builder from Python dicts / dataclasses
- CARP / VTX output parsers (read `.pts`, `.elem`, `.vtx` back into NumPy arrays)
- Packaging: pip-installable, separate repo or `python/` subdirectory TBD
- Integration with the broader pycemrg suite where it makes sense

**Depends on:** Phase 2b (done — parameter contract is stable).
**Blocks:** `v2.0.0` release tag — no merge to `master` until the parameter-file
builder is in place.

### Phase 5b — CI smoke test + Docker publish
Extensions to Phase 5a, deferred until a small fixture exists.

- Smoke test job: build + run `examples/sphereCoarse4Thickness` (or a smaller
  fixture) and assert non-empty outputs
- Publish updated Docker image on push to `master`

### Phase 6 — Push/PR CI (Linux-only)
`.github/workflows/ci.yml`: Linux build on push to `development`/`master` and on
PRs to either. Reuses the `release.yml` dependency install steps.

- Build-only, no packaging — fail-fast signal for refactors
- Open-source-friendly: external PRs get an automatic green/red check
- Future hook: wire in the Phase 5b smoke test once the fixture is ready

---

## Feature Backlog (post-phase)

These are smaller features carried over from the original `README.md` "TO DO" list
and items flagged as untested in the codebase. They are scoped for after the
phase work is complete.

### Output formats
- **Verify binary CARP output.** `out_carp_binary = 1` is implemented but has never
  been validated against a real openCARP consumer. Needs a round-trip test: write
  binary, read back in openCARP, confirm mesh identity.
- **`.mesh` writer inside `Mesh` class.** Currently only CGAL's `c3t3.output_to_medit`
  can produce `.mesh` output, which means the library cannot re-emit a `.mesh` file
  from an already-loaded `Mesh` object (useful when `--read_the_mesh` is used).
- **Triangle (boundary element) output.** Export the extracted boundary triangles as
  a standalone surface mesh (e.g. `.vtp`, `.obj`, or CARP surface `.elem`). Currently
  boundary information is consumed internally for VTX generation but never written
  out as geometry.

### Boundary processing
- **Triangle re-orientation with outward normals.** Mentioned in the README as
  implemented-but-unused. Needs a sanity pass to confirm it still works against the
  current `Mesh::extractBoundary()` output, and should become a CLI flag or a
  default for surface output.

---

## Phase Ordering (quick reference)

```
Phase 1+2a  DONE  Laplace pipeline extraction + laplace_solver app
Phase 4     DONE  CGAL 4.x to 6.x upgrade
Phase 2b    DONE  Meshing-side refactor (on development)
Phase 5a    DONE  Release CI + CPack packaging (on development, no tag yet)
Phase 6     NEXT  Push/PR CI, Linux-only
Phase 3           Python wrapper (blocks v2.0.0 tag)
Phase 2c          Optional: MeshingParams → m3d/, double fields
Phase 5b          CI smoke test + Docker publish on master
Backlog           Feature items from README TODOs
```
