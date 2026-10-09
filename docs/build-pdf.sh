#!/usr/bin/env bash
# Build Wavefront-MAT-Spec-v<version>.pdf from spec.md (pandoc + KaTeX + Chromium via Playwright).
# The version comes from the "version:" field in spec.md's front matter.
set -euo pipefail
cd "$(dirname "$0")"
version=$(sed -n 's/^version: *"\{0,1\}\([^"]*\)"\{0,1\}$/\1/p' spec.md | head -1)
[ -n "$version" ] || { echo "no version: field in spec.md" >&2; exit 1; }
[ -d node_modules/katex ] || npm install katex@0.16.11 --silent
pandoc spec.md -s --columns=1000 --katex=node_modules/katex/dist/ --highlight-style=tango -c spec.css --metadata pagetitle="Wavefront MAT" -o spec.html
SPEC_VERSION="$version" NODE_PATH="$(npm root -g)" node print.js
mv spec.pdf "Wavefront-MAT-Spec-v${version}.pdf"
rm -f spec.html
echo "built Wavefront-MAT-Spec-v${version}.pdf"
