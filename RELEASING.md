# Releasing

This repo's WDLs, docs, and Dockstore links reference themselves by git ref
(e.g. `.../raw/<ref>/python/...`). A release must **pin every such ref to the
release tag** so a tagged run is reproducible — while `master` keeps pointing at
`master` for day-to-day development.

To avoid carrying the pinned refs (and a revert commit) on `master`, cut the
release on a throwaway branch and **tag it without merging back**. The release
commit's parent is the `master` commit it was cut from, so the release stays
anchored to `master`'s history even though the tag itself is not an ancestor of
`master`.

## Steps

```bash
git checkout master              # start from a clean master (refs say "master")
git pull

git checkout -b release-tmp
./cut_release.sh vX.Y.Z          # pin all self-references to the tag
# review the diff, then:
git commit -am "Release vX.Y.Z: pin self-references"
git tag vX.Y.Z

git push origin vX.Y.Z           # push the TAG ONLY — never push the branch

git checkout master
git branch -D release-tmp        # the commit survives; the tag holds it
```

`master` is now untouched (still `master` refs, no revert commit). Then create
the GitHub release from the `vX.Y.Z` tag.

## Before tagging

- Bump `version:` / `date-released:` in `CITATION.cff` and add a `CHANGELOG.md`
  entry for the new version.
- Confirm `git status` is clean of local run artifacts and real input JSONs
  (`cromwell-executions/`, `execute/`, `*.input.full.json`, `cromwell_config/`).
- Publish before **20:00 US-Eastern** so the release's UTC date matches the
  intended calendar day (GitHub stores release timestamps in UTC).

## Notes

- **`git describe` / `git log master`** will not show the tag — expected, since
  the pinned tree deliberately diverges from `master`. Use `git checkout vX.Y.Z`
  or the GitHub `compare/<prev>...vX.Y.Z` view.
- **Dockstore** ingests the tag regardless of branch ancestry; its `:vX.Y.Z`
  version selector links resolve once ingest completes (a few minutes after the
  tag is pushed).
- To point everything back at a branch for local testing, just run
  `./cut_release.sh master` in your working tree (don't commit it to `master`).
