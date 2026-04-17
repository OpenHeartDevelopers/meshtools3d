# Build MeshTools3D on CGAL 4.x (Legacy Build)

This guide documents how to build MeshTools3D against CGAL 4.14.3 on Ubuntu 20.04.
The steps mirror the `Dockerfile.legacy-deps` and `Dockerfile.legacy` build process.

All libraries are installed under `$HOME/installs/` to keep them isolated from system packages.

---

## 1. System packages

Install the base toolchain and Boost from the Ubuntu package repository:

```bash
sudo apt-get update && sudo apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    ca-certificates \
    m4 \
    python3 \
    libssl-dev \
    libboost-all-dev
```

---

## 2. Environment variables

Set these before building. Add them to your shell profile or export them in your build session:

```bash
export GMP_PREFIX=$HOME/installs/gmp-6.2.1
export MPFR_PREFIX=$HOME/installs/mpfr-4.1.0
export TBBROOT=$HOME/installs/tbb-2020.3
export CGAL_DIR=$HOME/installs/cgal-4.14.3/lib/cmake/CGAL

export LD_LIBRARY_PATH="${GMP_PREFIX}/lib:${MPFR_PREFIX}/lib:${TBBROOT}/lib/intel64/gcc4.8"
export PATH="$HOME/installs/cmake-3.16.9/bin:${PATH}"
```

---

## 3. CMake 3.16.9 (pre-compiled binary)

```bash
mkdir -p $HOME/installs/cmake-3.16.9
wget https://cmake.org/files/v3.16/cmake-3.16.9-Linux-x86_64.tar.gz -O /tmp/cmake.tar.gz
tar -xzf /tmp/cmake.tar.gz -C $HOME/installs/cmake-3.16.9 --strip-components=1
rm /tmp/cmake.tar.gz
```

Verify the `PATH` export above is active; `cmake --version` should report `3.16.9`.

---

## 4. zlib 1.2.13

```bash
mkdir -p $HOME/installs/zlib-src
wget https://github.com/madler/zlib/releases/download/v1.2.13/zlib-1.2.13.tar.gz \
    -O /tmp/zlib.tar.gz
tar -xzf /tmp/zlib.tar.gz -C $HOME/installs/zlib-src --strip-components=1
cd $HOME/installs/zlib-src
./configure --prefix=$HOME/installs/zlib-1.2.13
make -j$(nproc) && make install
cd / && rm -rf $HOME/installs/zlib-src /tmp/zlib.tar.gz
```

---

## 5. GMP 6.2.1

```bash
mkdir -p $HOME/installs/gmp-src
wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz -O /tmp/gmp.tar.xz
tar -xJf /tmp/gmp.tar.xz -C $HOME/installs/gmp-src --strip-components=1
cd $HOME/installs/gmp-src
./configure --prefix=${GMP_PREFIX} --enable-cxx
make -j$(nproc) && make install
cd / && rm -rf $HOME/installs/gmp-src /tmp/gmp.tar.xz
```

---

## 6. MPFR 4.1.0

Depends on GMP (step 5 must be complete first).

```bash
mkdir -p $HOME/installs/mpfr-src
wget https://www.mpfr.org/mpfr-4.1.0/mpfr-4.1.0.tar.xz -O /tmp/mpfr.tar.xz
tar -xJf /tmp/mpfr.tar.xz -C $HOME/installs/mpfr-src --strip-components=1
cd $HOME/installs/mpfr-src
./configure --prefix=${MPFR_PREFIX} --with-gmp=${GMP_PREFIX}
make -j$(nproc) && make install
cd / && rm -rf $HOME/installs/mpfr-src /tmp/mpfr.tar.xz
```

---

## 7. Eigen 3.3.9 (header-only)

```bash
mkdir -p $HOME/installs/eigen-3.3.9
wget https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz \
    -O /tmp/eigen.tar.gz
tar -xzf /tmp/eigen.tar.gz -C $HOME/installs/eigen-3.3.9 --strip-components=1
rm /tmp/eigen.tar.gz
```

---

## 8. TBB 2020.3 (pre-compiled binary)

```bash
mkdir -p $HOME/installs/tbb-2020.3
wget https://github.com/uxlfoundation/oneTBB/releases/download/v2020.3/tbb-2020.3-lin.tgz \
    -O /tmp/tbb.tar.gz
tar -xzf /tmp/tbb.tar.gz -C $HOME/installs/tbb-2020.3 --strip-components=1
rm /tmp/tbb.tar.gz
```

Libraries are located at `${TBBROOT}/lib/intel64/gcc4.8`, which is already covered by the
`LD_LIBRARY_PATH` export in step 2.

---

## 9. CGAL 4.14.3

Depends on GMP, MPFR, Eigen, zlib, TBB, and Boost.

```bash
mkdir -p $HOME/installs/cgal-src $HOME/installs/cgal-4.14.3-build
wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.14.3/CGAL-4.14.3.tar.xz \
    -O /tmp/cgal.tar.xz
tar -xJf /tmp/cgal.tar.xz -C $HOME/installs/cgal-src --strip-components=1
cd $HOME/installs/cgal-4.14.3-build
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$HOME/installs/cgal-4.14.3 \
    -DGMP_INCLUDE_DIR=${GMP_PREFIX}/include \
    -DGMP_LIBRARIES=${GMP_PREFIX}/lib/libgmp.so \
    -DMPFR_INCLUDE_DIR=${MPFR_PREFIX}/include \
    -DMPFR_LIBRARIES=${MPFR_PREFIX}/lib/libmpfr.so \
    -DEIGEN3_INCLUDE_DIR=$HOME/installs/eigen-3.3.9 \
    -DZLIB_ROOT=$HOME/installs/zlib-1.2.13 \
    -DTBB_ROOT=${TBBROOT} \
    -DTBB_ARCH_PLATFORM=intel64/gcc4.8 \
    $HOME/installs/cgal-src
make -j$(nproc) && make install
cd / && rm -rf $HOME/installs/cgal-src $HOME/installs/cgal-4.14.3-build /tmp/cgal.tar.xz
```

---

## 10. Build MeshTools3D

Clone or copy the MeshTools3D source, then:

```bash
mkdir -p $HOME/installs/meshtools3D_build
cd $HOME/installs/meshtools3D_build
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCGAL_DIR=${CGAL_DIR} \
    -DGMP_INCLUDE_DIR=${GMP_PREFIX}/include \
    -DGMP_LIBRARIES=${GMP_PREFIX}/lib/libgmp.so \
    -DMPFR_INCLUDE_DIR=${MPFR_PREFIX}/include \
    -DMPFR_LIBRARIES=${MPFR_PREFIX}/lib/libmpfr.so \
    -DTBB_ROOT=${TBBROOT} \
    -DTBB_ARCH_PLATFORM=intel64/gcc4.8 \
    /path/to/meshtools3d
make -j$(nproc)
```

Executables are placed in `$HOME/installs/meshtools3D_build/`.

---

## Notes

- **Boost** is sourced from the Ubuntu system package (`libboost-all-dev`) rather than
  built from source. If you need a specific Boost version, build it from
  [boost.org](https://www.boost.org/) and add `-DBOOST_ROOT=<prefix>` to the CMake invocations.
- **LEDA** is not required for this build configuration.
- All paths above use `$HOME/installs/` as a prefix. Adjust to match your environment, keeping the
  environment variables in step 2 consistent.
