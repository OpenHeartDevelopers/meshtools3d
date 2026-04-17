#!/usr/bin/env bash
# bundle.sh — turn a meshtools3d build directory into a relocatable bundle.
#
# Given a build directory where cmake was run (binaries live at the root),
# walk each binary's dynamic dependencies, copy non-system libraries into
# <build-dir>/lib, and rewrite install-names / rpaths so the binary runs on
# another machine with only system libc + libstdc++/libc++ available.
#
# Usage: scripts/bundle.sh <build-dir>
#
# Platforms: macOS (otool + install_name_tool) and Linux (ldd + patchelf).

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <build-dir>" >&2
  exit 2
fi

BUILD_DIR="$1"
BIN_DIR="$BUILD_DIR"
LIB_DIR="$BUILD_DIR/lib"

if [[ ! -d "$BUILD_DIR" ]]; then
  echo "error: $BUILD_DIR does not exist" >&2
  exit 1
fi

mkdir -p "$LIB_DIR"

OS="$(uname -s)"

# Paths whose libraries we consider "system" and therefore skip.
is_system_lib_macos() {
  case "$1" in
    /usr/lib/*|/System/*|@rpath/*|@loader_path/*|@executable_path/*) return 0 ;;
    *) return 1 ;;
  esac
}

# Library basenames treated as "system" on Linux regardless of path.
is_system_lib_linux() {
  case "$(basename "$1")" in
    libc.so.*|libm.so.*|libdl.so.*|libpthread.so.*|librt.so.*|\
    libstdc++.so.*|libgcc_s.so.*|ld-linux-*.so.*|libnsl.so.*|\
    libresolv.so.*|libutil.so.*|linux-vdso.so.*|libcrypt.so.*)
      return 0 ;;
    *) return 1 ;;
  esac
}

bundle_macos() {
  local target="$1"
  local queue=("$target")
  local seen_file
  declare -A seen=()

  # ensure binary has an rpath entry pointing at ../lib
  install_name_tool -add_rpath "@loader_path/../lib" "$target" 2>/dev/null || true

  while [[ ${#queue[@]} -gt 0 ]]; do
    local current="${queue[0]}"
    queue=("${queue[@]:1}")

    # otool -L prints the file itself as first line, deps after
    while IFS= read -r line; do
      # strip leading tabs + trim "(compatibility ...)" suffix
      local dep
      dep="$(echo "$line" | sed -E 's/^\t//; s/ \(compatibility.*//')"
      [[ -z "$dep" ]] && continue
      [[ "$dep" == "$current:" ]] && continue
      is_system_lib_macos "$dep" && continue

      local base dest
      base="$(basename "$dep")"
      dest="$LIB_DIR/$base"

      # rewrite the reference in 'current' to @rpath/<base>
      install_name_tool -change "$dep" "@rpath/$base" "$current" 2>/dev/null || true

      if [[ -z "${seen[$base]+x}" ]]; then
        seen[$base]=1
        if [[ ! -f "$dest" ]]; then
          # resolve symlinks so we copy the actual file
          local real
          real="$(readlink -f "$dep" 2>/dev/null || echo "$dep")"
          cp "$real" "$dest"
          chmod u+w "$dest"
          install_name_tool -id "@rpath/$base" "$dest" 2>/dev/null || true
        fi
        queue+=("$dest")
      fi
    done < <(otool -L "$current" | tail -n +2)
  done
}

bundle_linux() {
  local target="$1"

  # Set rpath on the binary itself: look first in ../lib relative to bin/
  patchelf --set-rpath '$ORIGIN/../lib' "$target"

  # ldd gives the full transitive closure — walk it once, filter, copy.
  local line dep base dest
  while IFS= read -r line; do
    # lines look like: "<soname> => <path> (0x...)"  OR  "<path> (0x...)"
    if [[ "$line" =~ \=\>[[:space:]]+(/[^[:space:]]+) ]]; then
      dep="${BASH_REMATCH[1]}"
    elif [[ "$line" =~ ^[[:space:]]*(/[^[:space:]]+) ]]; then
      dep="${BASH_REMATCH[1]}"
    else
      continue
    fi
    is_system_lib_linux "$dep" && continue

    base="$(basename "$dep")"
    dest="$LIB_DIR/$base"
    if [[ ! -f "$dest" ]]; then
      cp "$(readlink -f "$dep")" "$dest"
      chmod u+w "$dest"
      # bundled libs look in the same directory for each other
      patchelf --set-rpath '$ORIGIN' "$dest" 2>/dev/null || true
    fi
  done < <(ldd "$target")
}

for exe in meshtools3d laplace_solver; do
  exe_path="$BIN_DIR/$exe"
  [[ -x "$exe_path" ]] || continue
  echo "bundling: $exe_path"
  case "$OS" in
    Darwin) bundle_macos "$exe_path" ;;
    Linux)  bundle_linux  "$exe_path" ;;
    *) echo "error: unsupported OS: $OS" >&2; exit 1 ;;
  esac
done

echo
echo "bundle complete. layout:"
echo "  $BIN_DIR/{meshtools3d,laplace_solver}"
echo "  $LIB_DIR"
ls -1 "$LIB_DIR" 2>/dev/null | sed 's/^/    /' || echo "    (no bundled shared libs)"
echo
case "$OS" in
  Darwin) echo "verify with: otool -L $BUILD_DIR/meshtools3d" ;;
  Linux)  echo "verify with: ldd   $BUILD_DIR/meshtools3d" ;;
esac
