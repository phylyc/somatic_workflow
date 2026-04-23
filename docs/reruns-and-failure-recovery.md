# Reruns and failure recovery

This page describes the most common rerun patterns.

## General guidance

The workflow emits many outputs that can be reused as cache inputs, and Terra/Cromwell call caching may help avoid recomputation. Even so, if core inputs or assumptions change, rerunning is usually safer than assuming previous intermediates remain valid.

## Main rerun patterns

### 1. Rerun after ABSOLUTE solution selection

This is the main expected rerun path.

Process:

1. run the workflow from BAMs/resources
2. inspect `absolute_acr_plot`
3. create `Cache.absolute_solution`
4. rerun with cached ABSOLUTE-related outputs plus `Cache.absolute_solution`

See [ABSOLUTE review and rerun](absolute-review-and-rerun.md).

### 2. Rerun after failed Mutect2 shards

The workflow exposes two shard-level controls for handling problematic Mutect2 shards:

- `Cache.skip_shards`
- `Cache.high_mem_shards`

These are arrays of 0-based shard indices.

Use cases:

- a specific shard silently ran out of memory and should be rerun with doubled memory
- one or a few regions consistently fail due to local BAM/pathology issues and need to be skipped

This avoids globally increasing memory or rerunning all shards unnecessarily.

## Practical advice for shard failures

- use `Cache.high_mem_shards` first when failure looks like a silent OOM
- use `Cache.skip_shards` only when a region remains persistently problematic and cannot be fixed upstream
- document any skipped shards in project notes because skipped regions imply intentionally incomplete SNV discovery in those intervals

## Cache conservatism

As a rule of thumb, if something important changes, rerun. Prefer reproducibility and interpretability over aggressive reuse of partially stale intermediates.
