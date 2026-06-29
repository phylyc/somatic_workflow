# Failure recovery

Troubleshooting for some potential workflow failures

### 1. Failed Mutect2 shards

The workflow exposes two shard-level controls for handling problematic Mutect2 shards:

- `Cache.high_mem_shards`: a specific shard silently ran out of memory and should be rerun with doubled memory. This avoids globally increasing memory or rerunning all shards unnecessarily.
- `Cache.skip_shards`: one or a few regions consistently fail due to local BAM issues and need to be skipped

**NOTE**
- use `Cache.high_mem_shards` first when failure looks like a silent OOM
- use `Cache.skip_shards` only when a region remains persistently problematic and cannot be fixed upstream
- document any skipped shards in project notes because skipped regions imply intentionally incomplete SNV discovery in those intervals

### 2. Failed Funcotator run

- Some long indels might cause Funcotator failure. In this case, it is best to manually remove the offending variant from the input VCF file and then re-run.

### Cache conservatism

As a rule of thumb, if something important changes, rerun. Prefer reproducibility and interpretability over aggressive reuse of partially stale intermediates.
