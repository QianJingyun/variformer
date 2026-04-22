#!/usr/bin/env bash
set -euo pipefail

# Run from anywhere inside the repo.
repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

echo "[check] scanning staged additions/modifications for forbidden artifacts..."

staged_files="$(git diff --cached --diff-filter=ACMR --name-only)"

if [[ -z "${staged_files}" ]]; then
  echo "[ok] no staged additions/modifications."
  exit 0
fi

forbidden_regex='^(data/|results/|logs/)|\.(ckpt|pt|pth|parquet|vcf|vcf\.gz|bcf|fa|fasta|fna|fai|bed|bed\.gz|gct|gct\.gz|zip)$'

violations="$(printf '%s\n' "${staged_files}" | awk 'NF > 0' | grep -E "${forbidden_regex}" || true)"

if [[ -n "${violations}" ]]; then
  echo "[fail] found staged files that should not be pushed:"
  printf '%s\n' "${violations}"
  echo
  echo "Fix:"
  echo "  - unstage file: git restore --staged <path>"
  echo "  - stop tracking folder while keeping local files: git rm -r --cached <folder>"
  exit 1
fi

echo "[ok] staged files look safe for GitHub push."
