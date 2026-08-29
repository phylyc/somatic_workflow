# Failure recovery

Troubleshooting for some potential workflow failures

### 1. Failed Mutect2 shards

The workflow exposes three shard-level controls for handling problematic Mutect2 shards:

- `z2_Cache.high_mem_shards`: a specific shard silently ran out of memory and should be rerun with doubled memory. This avoids globally increasing memory or rerunning all shards unnecessarily.
- `z2_Cache.huge_mem_shards`: a shard failed again at the high-memory tier and should be rerun with 4x the default memory. If a shard appears in both lists, the huge-memory tier takes precedence.
- `z2_Cache.skip_shards`: one or a few regions consistently fail due to local BAM issues and need to be skipped

**NOTE**
- use `z2_Cache.high_mem_shards` first when failure looks like a silent OOM
- use `z2_Cache.huge_mem_shards` for the subset that still fails at the high-memory tier
- use `z2_Cache.skip_shards` only when a region remains persistently problematic and cannot be fixed upstream
- document any skipped shards in project notes because skipped regions imply intentionally incomplete SNV discovery in those intervals

### 2. Failed Funcotator run

- Some long indels might cause Funcotator failure. In this case, it is best to manually remove the offending variant from the input VCF file and then re-run.

### Cache conservatism

As a rule of thumb, if something important changes, rerun. Prefer reproducibility and interpretability over aggressive reuse of partially stale intermediates.
