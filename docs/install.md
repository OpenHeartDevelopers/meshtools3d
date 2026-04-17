# Installing MeshTools3D

This guide covers building MeshTools3D on macOS and Linux using system
package managers. For the legacy CGAL 4.x from-source build on older
systems, see [install_legacy.md](install_legacy.md).

## Dependencies

- CMake 3.14+
- C++17 compiler (GCC 9+, Clang 12+, Apple Clang 14+)
- CGAL 6.x (header-only; pulls in GMP and MPFR transitively)
- Boost
- zlib
- TBB (optional, for parallel mesh generation)

---

## macOS

Tested on macOS 13+ (Ventura) with Homebrew. Works on both Apple Silicon and Intel.

### 1. Install dependencies

```bash
brew install cmake cgal boost tbb
```

### 2. Build

```bash
git clone <repo-url> meshtools3d && cd meshtools3d
mkdir build && cd build
cmake ..
make -j$(sysctl -n hw.logicalcpu)
```

If CMake cannot find CGAL (rare on Apple Silicon), pass the prefix explicitly:

```bash
cmake -DCGAL_DIR=$(brew --prefix cgal)/lib/cmake/CGAL ..
```

### 3. Verify

```bash
./meshtools3d -h
./laplace_solver -h
```

---

## Linux (Ubuntu 22.04+ / Debian 12+)

### 1. Install dependencies

```bash
sudo apt-get update && sudo apt-get install -y \
    build-essential cmake \
    libcgal-dev libboost-all-dev libtbb-dev zlib1g-dev
```

### 2. Build

```bash
git clone <repo-url> meshtools3d && cd meshtools3d
mkdir build && cd build
cmake ..
make -j$(nproc)
```

---

## Linux (Fedora 38+ / RHEL 9+)

### 1. Install dependencies

```bash
sudo dnf install -y \
    gcc-c++ cmake CGAL-devel boost-devel tbb-devel zlib-devel
```

### 2. Build

Same as Ubuntu from step 2 onward.

---

## Build Variants

**Debug build** (enables bounds-checked array access via `.at()`):

```bash
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DCGAL_DONT_OVERRIDE_CMAKE_FLAGS=FALSE \
      -DCGAL_DISABLE_ROUNDING_MATH_CHECK=ON ..
```

**Without TBB** (single-threaded mesh generation):

```bash
cmake -DCMAKE_DISABLE_FIND_PACKAGE_TBB=ON ..
```

---

## Parallelism

Control thread count at runtime via the `TBB_NUM_THREADS` environment variable.
Defaults to 1 thread if unset.

```bash
export TBB_NUM_THREADS=4
./meshtools3d -f data
```

---

## Docker

```bash
docker build -t meshtools3d .
docker run -v $(pwd)/my_data:/data meshtools3d -f /data/data
```

Mount your input data directory to `/data` inside the container.

---

## Troubleshooting

**CMake cannot find CGAL** -- Pass `-DCGAL_DIR=<path>` pointing to the directory
containing `CGALConfig.cmake`. Common locations: `/opt/homebrew/lib/cmake/CGAL`
(macOS ARM), `/usr/local/lib/cmake/CGAL` (macOS Intel),
`/usr/lib/x86_64-linux-gnu/cmake/CGAL` (Ubuntu).

**Rounding math errors at runtime** -- Add `-DCGAL_DISABLE_ROUNDING_MATH_CHECK=ON`
to your CMake invocation. Release builds set `-frounding-math` automatically.

**Linker errors about GMP or MPFR** -- These are transitive CGAL dependencies.
Install explicitly if needed: `brew install gmp mpfr` (macOS) or
`sudo apt install libgmp-dev libmpfr-dev` (Ubuntu).
