#!/bin/bash
# Safely fast-forward an existing Biowulf V2 checkout from GitHub.

set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

if [[ ! -d .git ]]; then
  echo "This directory is not a Git checkout: $repo_root" >&2
  exit 2
fi
if [[ -n "$(git status --porcelain --untracked-files=normal)" ]]; then
  echo "Tracked or untracked code files have local changes. Review 'git status' before pulling." >&2
  exit 3
fi

git fetch origin main
git pull --ff-only origin main
chmod u+x biowulf/*.sh biowulf/*.sbatch
echo "Updated V2 to $(git rev-parse --short HEAD). Ignored data/runtime files were preserved."
