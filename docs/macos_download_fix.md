# macOS (Apple Silicon): first-run setup for downloaded releases

The macOS `arm64` binaries are **ad-hoc signed but not notarized**. When you
download the release archive through a browser, macOS tags it with the
`com.apple.quarantine` attribute, and Gatekeeper then blocks the tools: the
process is silently terminated with `zsh: killed` (`killed: 9`) and no error
message.

This affects browser downloads from the
[Releases page](https://github.com/OpenHeartDevelopers/meshtools3d/releases).
It is a one-time setup per download.

## Fix (recommended)

Each macOS archive ships a `macos-fix-download.sh` helper at its root. After
extracting, run it and pass the path to the extracted folder:

```bash
DIR=~/installs/meshtools3d-<version>-macos-arm64
"$DIR/macos-fix-download.sh" "$DIR"
```

It clears the quarantine flag, re-signs the bundled libraries and executables
(libraries first), and verifies the result.

## Fix (manual fallback)

If you don't have the helper script, run the equivalent steps by hand. Set
`DIR` to wherever you extracted the archive:

```bash
# Path to the extracted folder
DIR=~/installs/meshtools3d-<version>-macos-arm64

# 1. Remove the "downloaded from the internet" quarantine flag
xattr -dr com.apple.quarantine "$DIR"

# 2. Re-sign the bundled libraries (skips symlinks)
for f in "$DIR"/lib/*.dylib; do
  [ -L "$f" ] && continue
  codesign --force --sign - "$f"
done

# 3. Re-sign the executables (must be done after the libraries)
for b in "$DIR"/bin/*; do
  [ -L "$b" ] && continue
  codesign --force --sign - "$b"
done
```

## Verify

```bash
"$DIR/bin/meshtools3d" --help
```

You should see the `MeshTools3D` usage text. The tools will now run normally.

## Why this is needed

The binaries are signed ad-hoc during packaging, so they are valid on the build
machine. Two things can still block a browser download:

- **Quarantine.** macOS adds `com.apple.quarantine` to anything downloaded via a
  browser; Gatekeeper refuses to run ad-hoc (un-notarized) code while that flag
  is set. Removing it is the main fix.
- **Broken signatures (fallback).** The bundled libraries (`libgmp`, `libmpfr`,
  `libtbb`, …) have their load paths rewritten during packaging
  (`install_name_tool`), which invalidates code signatures. The packaging step
  re-signs them, but if an archive ever loses its signatures in transit, the
  re-sign above restores them. On Apple Silicon, macOS refuses to load a library
  with an invalid signature and kills the process before it can print anything.

A future release will ship pre-signed **and notarized** binaries so these steps
are no longer necessary.

## Diagnosing (optional)

To confirm the cause on a fresh extraction:

```bash
DIR=~/installs/meshtools3d-<version>-macos-arm64

# (a) Quarantine flag present?
xattr -lr "$DIR" | grep -c quarantine

# (b) Broken dylib signatures? (the fatal issue)
for f in "$DIR"/lib/*.dylib; do
  printf "%s: " "$f"
  codesign --verify "$f" 2>&1 || true
  echo
done
```

A nonzero count in (a), or `invalid signature (code or signature have been
modified)` for a dylib in (b), confirms the issue.
