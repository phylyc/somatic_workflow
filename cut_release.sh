#!/usr/bin/env bash
#
# cut_release.sh — repoint every self-referential link to this repository to a
# single git ref, in preparation for tagging a release.
#
# Three kinds of self-reference carry a git ref and are rewritten:
#   1. raw script fetches   https://github.com/phylyc/somatic_workflow/raw/<ref>/python/*.py
#                           (the WDLs download these at runtime — must be pinned
#                            to an immutable tag for a reproducible release)
#   2. Dockstore selectors  https://dockstore.org/workflows/github.com/phylyc/somatic_workflow/<Name>:<ref>
#   3. source links         https://github.com/phylyc/somatic_workflow/blob/<ref>/...  (and /tree/<ref>/...)
#
# Usage:
#   ./cut_release.sh <git-ref>
#     e.g. ./cut_release.sh v2.0.0     # freeze everything to the release tag
#          ./cut_release.sh master     # point everything back at the dev branch
#
# Typical release flow:
#   ./cut_release.sh v2.0.0
#   git commit -am "Pin self-references to v2.0.0"
#   (merge to master, then) git tag v2.0.0 && git push --tags
# The raw/blob/Dockstore links resolve as soon as the tag exists.
#
set -euo pipefail

REPO="phylyc/somatic_workflow"

ref="${1:-}"
if [[ -z "$ref" ]]; then
  echo "Usage: $0 <git-ref>   (e.g. v2.0.0, master)" >&2
  exit 1
fi

# Refs with slashes would break the single-segment URL rewrites below.
if [[ ! "$ref" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "Refusing: ref must match [A-Za-z0-9._-]+ (no slashes)." >&2
  exit 1
fi

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Any file mentioning a ref-bearing self-reference. Exclude VCS and run outputs,
# and exclude *.sh so this script (which contains the patterns as literals) is
# never itself rewritten.
find_files() {
  grep -rlE "(github\.com/${REPO}/(raw|blob|tree)/|dockstore\.org/workflows/github\.com/${REPO}/)" "$root" \
    --exclude-dir=.git --exclude-dir=cromwell-executions --exclude-dir=execute \
    --exclude='*.sh' 2>/dev/null || true
}

mapfile -t files < <(find_files)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No self-references found."
  exit 0
fi

export REF="$ref"
for f in "${files[@]}"; do
  # \x27 = single quote, kept out of the ref token's stop-set alongside the
  # other URL/markdown delimiters so refs inside 'single-quoted' strings work.
  perl -i -pe 'BEGIN { $r = $ENV{REF} }
    s{(github\.com/phylyc/somatic_workflow/(?:raw|blob|tree)/)[^/\s"\x27`)>\]]+}{$1$r}g;
    s{(dockstore\.org/workflows/github\.com/phylyc/somatic_workflow/[A-Za-z0-9_]+:)[A-Za-z0-9._-]+}{$1$r}g;
  ' "$f"
  echo "  updated ${f#"$root"/}"
done

echo
echo "Repointed all self-references to '${ref}'. Resulting refs in tree:"
grep -rhoE "github\.com/${REPO}/(raw|blob|tree)/[A-Za-z0-9._-]+" "$root" \
  --exclude-dir=.git --exclude-dir=cromwell-executions --exclude-dir=execute --exclude='*.sh' \
  | sed -E 's#/[A-Za-z0-9._-]+$#/&#' | sort | uniq -c
grep -rhoE "dockstore\.org/workflows/github\.com/${REPO}/[A-Za-z0-9_]+:[A-Za-z0-9._-]+" "$root" \
  --exclude-dir=.git --exclude-dir=cromwell-executions --exclude-dir=execute --exclude='*.sh' \
  | sort | uniq -c
