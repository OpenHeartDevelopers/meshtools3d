# Installation

`meshtools3d` is a C++ CLI built on CGAL 5.6+ / 6.x. This document covers:

1. [Requirements](#requirements)
2. [Build on macOS](#build-on-macos-homebrew)
3. [Build on Ubuntu 22.04](#build-on-ubuntu-2204)
4. [Portable binary (rpath bundle)](#portable-binary-rpath-bundle)
5. [Runtime](#runtime)
6. [Troubleshooting](#troubleshooting)

---

## Requirements

| Component | Minimum | Notes |
|---|---|---|
| C++ compiler | C++17 | AppleClang 15+, GCC 9+, Clang 10+ |
| CMake | 3.14 | Homebrew / Ubuntu 22.04 ship newer |
| CGAL | 5.6 | 6.x works too; 5.x before 5.6 may miss `create_labeled_image_mesh_domain` |
| Boost | 1.70 | Header-only usage plus `filesystem` on Windows |
| oneTBB | 2021 | Optional but recommended; enables parallel meshing |
| zlib | — | Pulled in by CGAL's `ImageIO` component |

CGAL 6 dropped `Labeled_image_mesh_domain_3`; the build will not succeed against CGAL 4.x.

---

## Build on macOS (Homebrew)

```sh
brew install cmake cgal boost tbb

cd <repo-root>
cmake -S . -B ../build/MESHTOOLS3D_release -DCMAKE_BUILD_TYPE=Release
cmake --build ../build/MESHTOOLS3D_release -j
```

Produces `../build/MESHTOOLS3D_release/meshtools3d` and `.../laplace_solver`.

Debug build:

```sh
cmake -S . -B ../build/MESHTOOLS3D_debug \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCGAL_DISABLE_ROUNDING_MATH_CHECK=ON
cmake --build ../build/MESHTOOLS3D_debug -j
```

---

## Build on Ubuntu 22.04

```sh
sudo apt-get update
sudo apt-get install -y \
    build-essential cmake \
    libcgal-dev libboost-all-dev libtbb-dev zlib1g-dev

cd <repo-root>
cmake -S . -B ../build/MESHTOOLS3D_release -DCMAKE_BUILD_TYPE=Release
cmake --build ../build/MESHTOOLS3D_release -j
```

Ubuntu 22.04 ships CGAL 5.4, which **predates** the `create_labeled_image_mesh_domain` factory. If apt's CGAL is too old, install from source into `/opt/cgal` and pass `-DCGAL_DIR=/opt/cgal/lib/cmake/CGAL`. Ubuntu 24.04 (CGAL 5.6) works out of the box.

---

## Portable binary (rpath bundle)

The binary is meant to be shipped inside another application that invokes it as a subprocess. To produce a self-contained tree that runs on a target machine without Homebrew / apt, use the install + bundle recipe:

```sh
# 1. build as usual
cmake -S . -B ../build/MESHTOOLS3D_release -DCMAKE_BUILD_TYPE=Release
cmake --build ../build/MESHTOOLS3D_release -j

# 2. install into a staging prefix
cmake --install ../build/MESHTOOLS3D_release --prefix meshtools3d-dist

# 3. copy non-system shared libs into meshtools3d-dist/lib and rewrite rpaths
scripts/bundle.sh meshtools3d-dist
```

Resulting layout:

```
meshtools3d-dist/
├── bin/
│   ├── meshtools3d
│   └── laplace_solver
└── lib/
    ├── libtbb.12.dylib       (Linux: libtbb.so.12)
    ├── libtbbmalloc.2.dylib
    ├── libCGAL_ImageIO.*
    ├── libgmp.*
    └── libmpfr.*
```

### How the bundling works

- `CMakeLists.txt` sets `CMAKE_INSTALL_RPATH` to `@loader_path/../lib` on macOS and `$ORIGIN/../lib` on Linux. The installed binary looks for shared libraries there first.
- `scripts/bundle.sh` walks the binary's dynamic dependencies, filters out system libraries (`/usr/lib`, `/System`, libc/libm/libstdc++), copies the rest into `meshtools3d-dist/lib/`, and rewrites references:
  - **macOS** uses `otool -L` to inspect and `install_name_tool -change` / `-id` to rewrite each reference to `@rpath/<basename>`.
  - **Linux** uses `ldd` to get the transitive closure and `patchelf --set-rpath` to point each bundled lib at `$ORIGIN` (so they find each other).
- System libraries (`libc++` on macOS, `libstdc++`/`libc`/`libm` on Linux) are **not** bundled — they are ABI-stable across recent OS releases and bundling them usually causes more problems than it solves.

### Verify the bundle

```sh
# macOS
otool -L meshtools3d-dist/bin/meshtools3d

# Linux
ldd   meshtools3d-dist/bin/meshtools3d
```

Every non-system reference should be `@rpath/lib*.dylib` (macOS) or a path resolved under `meshtools3d-dist/lib/` via `$ORIGIN` (Linux). If you see an absolute path into `/opt/homebrew` or `/usr/local`, the bundler missed something — add that basename to the dependency walk in `scripts/bundle.sh`.

### Relocation test

```sh
mv meshtools3d-dist /tmp/ship-test
/tmp/ship-test/bin/meshtools3d -h
```

If the binary runs, the bundle is relocatable.

### Known limitations

- Built on one CPU architecture, runs on that architecture only. To ship both `arm64` and `x86_64` macOS binaries, build twice and combine with `lipo` (out of scope here).
- On Linux, the bundle is glibc-forward-compatible: building on Ubuntu 22.04 (glibc 2.35) runs on 22.04+ but not on older distros. To target older systems, build inside a container based on that system.
- macOS code-signing / notarization is not automated; sign the binary if shipping outside a containing app that is already signed.

---

## Runtime

Inputs: an INR-format segmentation plus a GetPot `data` file (see `data-template`).

```sh
cp data-template data
meshtools3d -f data
```

Override single options without editing `data`:

```sh
meshtools3d -f data -seg_dir /path/to/seg -seg_name heart.inr \
            -out_dir /path/to/out -out_name heart \
            --thickness-algorithm 2 -v
```

Set TBB thread count via environment:

```sh
export TBB_NUM_THREADS=8
meshtools3d -f data
```

The separate `laplace_solver` binary solves Laplace BCs on a pre-existing mesh (no CGAL dependency):

```sh
laplace_solver -f data --zero-bc vtx/endo.vtx --one-bc vtx/epi.vtx
```

---

## Troubleshooting

**`dyld: Library not loaded: @rpath/libtbb.12.dylib`** — the bundle is missing that library or the binary's rpath is wrong. Re-run `scripts/bundle.sh`; check `otool -l bin/meshtools3d | grep -A2 RPATH`.

**`error while loading shared libraries: libCGAL_ImageIO.so.14: cannot open shared object file`** — same story on Linux. Check `readelf -d bin/meshtools3d | grep RPATH` and `ls lib/`.

**`Failed to find Boost`** under CMake 4.x with a recent Boost — the legacy `FindBoost` module was removed. The project sets `CMP0167 NEW` to opt into `BoostConfig.cmake`; if your system Boost doesn't provide `BoostConfig.cmake`, install a newer Boost (Homebrew and Ubuntu 24.04 do; Ubuntu 22.04's apt Boost does too via `libboost-all-dev`).

**`call to 'lloyd' is ambiguous`** — CGAL 6 exposes two `parameters::lloyd()` overloads that collide on zero-arg calls. `main.cpp` calls `parameters::lloyd<bool>()` (and siblings) to disambiguate. Do not remove the explicit template argument.

**`TBB_NUM_THREADS not set; nb of threads is: 1 (default)`** — informational, not an error. Export the variable if you want parallelism.

**Build succeeds, run crashes with a segfault inside mesh generation** — most often a Debug/Release ABI mismatch. `#define NDEBUG` toggles `std::vector::at` vs `operator[]` behavior in the library; don't mix Debug and Release objects in one binary.
