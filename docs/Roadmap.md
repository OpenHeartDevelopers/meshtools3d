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

### Phase 3a — C++ `parfile_builder` (in this repo)
- Lean stdlib-only application (`applications/parfile_builder.cpp`) that emits a
  complete meshtools3d `.par` file with documented defaults. No CGAL, no `m3d/` link.
- CLI: defaults always baked in; `--seg_dir` / `--seg_name` / `--out_dir` /
  `--out_name` quality-of-life flags; `--set SECTION.KEY=VALUE` for arbitrary
  overrides; output to stdout by default, `-o PATH` to write a file.
- Bundled into the CPack tarballs alongside `meshtools3d` and `laplace_solver`.
- **Status:** on `development`, ships in the v2.0.0-beta.2 binaries.
- Schema documented in `docs/parameter_file_schema.md` (`parfile_builder` is the
  source of truth).

### Phase 3b — `pycemrg-meshing` (out of repo) — LIVE
Python package wrapping the meshtools3d binaries; lives in a separate repo,
imports `pycemrg`, uses `ModelManager` to fetch versioned binaries from the
GitHub Releases page, and vendors the parameter-file writer in Python.
**Already released and in use** — out of scope for this repo and not a v2.0.0
blocker. Spec retained at `.claude/ticket_pycemrg_meshing.md`.

### Phase 4 — CGAL 4.x to 6.x upgrade
- Modernised `CMakeLists.txt` to target-based linking (`CGAL::CGAL`, `TBB::tbb`)
- Updated `CGALDataType.hpp` and the custom domain wrappers for the CGAL 5/6 API
- TBB migrated from `task_scheduler_init` to `tbb::global_control`
- Dockerfiles updated; legacy CGAL 4.14.3 from-source install preserved in
  `Dockerfile.legacy` and `docs/install_legacy.md`
- **Status:** merged to `master` via `feature/cgal_upgrade`

### Phase 5a — Release CI + packaging
- `.github/workflows/release.yml`: Linux/macOS/Windows build matrix, CGAL 6.1.1,
  triggered on `release: published`
- `CMakeLists.txt`: CPack, install rpath, Homebrew dylib bundling (macOS),
  vcpkg DLL bundling (Windows), system-lib deps documented for Linux
- SHA256 checksum sidecars (`.sha256`) computed and uploaded per artifact on all
  platforms, for `pycemrg.ModelManager` verification (was tracked separately in
  `.claude/ticket_release_ci_for_model_manager.md` — now done)
- **macOS code-signing fix:** the `install_name_tool` rpath rewrites invalidate
  each bundled dylib's code signature; on Apple Silicon AMFI then SIGKILLs the
  process on launch ("killed: 9"). The install step now re-signs all dylibs and
  executables ad-hoc (`codesign --force --sign -`) after the rewrites and after
  CPack stripping. Workaround for already-published beta artifacts documented in
  `docs/macOS_bug_fix_v2.0-beta.md`.
- **Status:** on `development`, exercised via the `v2.0.0-beta.2` tag.

### Phase 6 — Push/PR CI (Linux-only)
- `.github/workflows/ci.yml`: Linux build on push to `development`/`master` and on
  PRs to either; reuses the `release.yml` dependency install steps
- Build-only, no packaging — fail-fast signal for refactors
- **Status:** done

### Documentation
- Promoted `m3d_python_params.md` to `docs/parameter_file_schema.md` — the public
  parameter-file schema reference (the `.par` contract shared by `parfile_builder`
  and downstream consumers). `parfile_builder` is the source of truth; the doc
  matches its built-in defaults.
- Added `docs/macOS_bug_fix_v2.0-beta.md` — first-run workaround for the
  already-published beta.2 macOS binaries.

---

## Remaining before v2.0.0 (non-beta)

The C++/CLI and library surface is feature-complete for v2.0. The Python wrapper
is **out of scope for this repo** (separate `pycemrg-meshing` package) and does
**not** block the v2.0.0 tag,  it ships and versions independently.

### Release housekeeping
- Merge `development` → `master` (carries Phase 2b, 3a, 5a, 6, and the macOS
  code-signing fix).
- Write v2.0.0 final release notes (the macOS section can drop the manual
  `xattr` / `codesign` workaround now that the CMake fix ships in-build).

### Windows build — deferred
- The Windows job is `continue-on-error` / experimental pending a GetPot /
  `boost::filesystem` fix (see `.claude/ticket_windows_build_fix.md`).
- **Decision (JA):** deferred. Only one Windows 10 machine is available and it
  holds several golden artifacts, so we return to this when it's feasible.
  v2.0.0 ships with **Windows documented as best-effort**; Linux and macOS are
  the supported platforms.

---

## Remaining Phases (post-v2.0 / optional)

### Phase 2c — Move `MeshingParams` out of CGAL
> To be considered as an option. Not decided!
- Change numeric fields from `FaceNumericalType`/`CellNumericalType` to `double`
- Move the struct from `CGALDataType.hpp` to `m3d/include/` alongside `LaplaceParams`
- `main.cpp` passes `double` values to CGAL criteria constructors directly (they accept double)

The only risk is if a future CGAL version makes the criteria constructors finicky about
`double` vs. their own `FT` — unlikely but worth noting.

### Phase 5b — CI smoke test + Docker publish
Extensions to Phase 5a, deferred until a small fixture exists.

- Smoke test job: build + run `examples/sphereCoarse4Thickness` (or a smaller
  fixture) and assert non-empty outputs
- Publish updated Docker image on push to `master`

---

## Ticketed for `pycemrg-meshing` (out of repo)
- **CARP / VTX → NumPy parsers.** Read `.pts`, `.elem`, `.vtx` back into NumPy
  arrays from Python. Belongs to `pycemrg-meshing`, not this repo. Spec drafted
  at `.claude/ticket_carp_vtx_numpy_parsers.md` for JA to move into the
  pycemrg-meshing repo (targets `pycemrg-meshing v0.2`).

---

## Feature Backlog (post-phase)

These are smaller features carried over from the original `README.md` "TO DO" list
and items flagged as untested in the codebase. They are scoped for after the
phase work is complete.

### Output formats
- **Verify binary CARP output.** `out_carp_binary = 1` is implemented but has never
  been validated against a real openCARP consumer. Needs a round-trip test: write
  binary, read back in openCARP, confirm mesh identity.
  **Priority: low (JA).**
- **`.mesh` writer inside `Mesh` class.** Currently only CGAL's `c3t3.output_to_medit`
  can produce `.mesh` output, which means the library cannot re-emit a `.mesh` file
  from an already-loaded `Mesh` object (useful when `--read_the_mesh` is used).
- **Triangle (boundary element) output.** Export the extracted boundary triangles as
  a standalone surface mesh (e.g. `.vtp`, `.obj`, or CARP surface `.elem`). Currently
  boundary information is consumed internally for VTX generation but never written
  out as geometry.
  **Priority: useful (JA)** — good near-term pickup.

### Boundary processing
- **Triangle re-orientation with outward normals.** Mentioned in the README as
  implemented-but-unused. Needs a sanity pass to confirm it still works against the
  current `Mesh::extractBoundary()` output, and should become a CLI flag or a
  default for surface output.

---

## Phase Ordering (quick reference)

```
Phase 1+2a  DONE  Laplace pipeline extraction + laplace_solver app (master)
Phase 4     DONE  CGAL 4.x to 6.x upgrade (master)
Phase 2b    DONE  Meshing-side refactor (development)
Phase 3a    DONE  C++ parfile_builder (development, ships in beta.2)
Phase 3b    LIVE  pycemrg-meshing Python package (out of repo, released)
Phase 5a    DONE  Release CI + CPack + SHA256 + macOS codesign fix (development)
Phase 6     DONE  Push/PR CI, Linux-only
Docs        DONE  docs/parameter_file_schema.md

── v2.0.0 (non-beta) ──────────────────────────────────────────────
  Merge development → master; final release notes
  Windows: deferred — ship best-effort (Linux + macOS supported)
  (Python wrapper is out-of-repo and does NOT block this tag)

Phase 2c          Optional: MeshingParams → m3d/, double fields
Phase 5b          CI smoke test + Docker publish on master
pycemrg           CARP / VTX → NumPy parsers (ticketed, out of repo)
Backlog           Feature items from README TODOs
```
