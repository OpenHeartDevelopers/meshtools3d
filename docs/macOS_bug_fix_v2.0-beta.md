# macOS (Apple Silicon): first-run setup for v2.0.0-beta

The macOS `arm64` binaries in this beta are **ad-hoc signed and not notarized**,
and the bundled libraries need to be re-signed after extraction. Without this
step macOS will silently terminate the tools with `zsh: killed` (`killed: 9`)
and no error message.

## Fix

Run the following **once** after extracting the archive. Set `DIR` to wherever
you extracted it:

```bash
# Path to the extracted folder
DIR=~/installs/meshtools3d-2.0.0-beta.2-macos-arm64

# 1. Remove the "downloaded from the internet" quarantine flag
xattr -dr com.apple.quarantine "$DIR"

# 2. Re-sign the bundled libraries (skips symlinks)
for f in "$DIR"/lib/*.dylib; do
  [ -L "$f" ] && continue
  codesign --force --sign - "$f"
done

# 3. Re-sign the executables (must be done after the libraries)
for b in meshtools3d laplace_solver parfile_builder; do
  codesign --force --sign - "$DIR/bin/$b"
done
```

## Verify

```bash
"$DIR/bin/meshtools3d" --help
```

You should see the `MeshTools3D` usage text. The tools will now run normally.

## Why this is needed

The bundled libraries (`libgmp`, `libmpfr`, `libtbb`, …) have their load paths
rewritten during packaging (`install_name_tool`), which invalidates their code
signatures. On Apple Silicon, macOS refuses to load libraries with invalid
signatures and kills the process before it can print anything. Re-signing
locally resolves it.

A future release will ship pre-signed (and notarized) binaries so these steps
are no longer necessary.

## Diagnosing (optional)

To confirm the cause on a fresh extraction:

```bash
DIR=~/installs/meshtools3d-2.0.0-beta.2-macos-arm64

# (a) Quarantine flag present?
xattr -lr "$DIR" | grep -c quarantine

# (b) Broken dylib signatures? (the fatal issue)
for f in "$DIR"/lib/*.dylib; do
  printf "%s: " "$f"
  codesign --verify "$f" 2>&1 || true
  echo
done
```

A nonzero count in (a) and `invalid signature (code or signature have been
modified)` for each dylib in (b) confirm the issue.
