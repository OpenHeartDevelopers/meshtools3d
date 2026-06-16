#!/usr/bin/env bash
# macos-fix-download.sh — make a downloaded macOS build runnable.
#
# When the release archive is downloaded through a browser, macOS tags it with
# the com.apple.quarantine attribute. Combined with ad-hoc (un-notarized)
# signatures, Gatekeeper blocks the tools and the process is terminated with
# "zsh: killed" (killed: 9) before printing anything.
#
# This script, run once on the extracted folder, removes the quarantine flag
# and re-signs the bundled libraries and executables ad-hoc so they launch.
#
# Usage:
#   ./macos-fix-download.sh <extracted-dir>
#
# Pass the path to the extracted folder explicitly, e.g.:
#   ./macos-fix-download.sh ~/installs/meshtools3d-<version>-macos-arm64

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <extracted-dir>" >&2
  echo "  e.g. $0 ~/installs/meshtools3d-<version>-macos-arm64" >&2
  exit 1
fi

DIR="$1"

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "This script only applies to macOS." >&2
  exit 1
fi

BIN_DIR="$DIR/bin"
LIB_DIR="$DIR/lib"

if [[ ! -d "$BIN_DIR" ]]; then
  echo "error: '$BIN_DIR' not found." >&2
  echo "Pass the extracted folder as an argument, e.g.:" >&2
  echo "  $0 ~/installs/meshtools3d-<version>-macos-arm64" >&2
  exit 1
fi

echo "Fixing macOS bundle at: $DIR"

# 1. Remove the "downloaded from the internet" quarantine flag.
echo "  - removing quarantine attribute"
xattr -dr com.apple.quarantine "$DIR" 2>/dev/null || true

# 2. Re-sign the bundled libraries first (skip symlinks).
if [[ -d "$LIB_DIR" ]]; then
  for f in "$LIB_DIR"/*.dylib; do
    [[ -e "$f" ]] || continue        # no dylibs present
    [[ -L "$f" ]] && continue        # skip symlinks
    echo "  - signing lib: $(basename "$f")"
    codesign --force --sign - "$f"
  done
fi

# 3. Re-sign the executables (must come after the libraries).
for b in "$BIN_DIR"/*; do
  [[ -f "$b" ]] || continue
  [[ -L "$b" ]] && continue
  echo "  - signing bin: $(basename "$b")"
  codesign --force --sign - "$b"
done

# 4. Verify every signed file; report a clear pass/fail summary.
echo "Verifying signatures..."
failed=0
for f in "$LIB_DIR"/*.dylib "$BIN_DIR"/*; do
  [[ -f "$f" ]] || continue
  [[ -L "$f" ]] && continue
  if ! codesign --verify --strict "$f" 2>/dev/null; then
    echo "  FAILED: $f"
    failed=1
  fi
done

if [[ "$failed" -ne 0 ]]; then
  echo "One or more files failed verification. See messages above." >&2
  exit 1
fi

echo "Done. The tools should now run, e.g.:"
echo "  $BIN_DIR/meshtools3d --help"
