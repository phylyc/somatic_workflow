"""Single source of truth for the per-task Runtime construction shape.

This drives the three mechanical bands of runtime_collection.wdl (the RuntimeCollection
struct fields, the `Runtime <name> = {...}` constructions, and the final collection map),
which are filled between `# CODEGEN:BEGIN <id>` / `# CODEGEN:END <id>` markers. The
workflow's `input {}` block — the actual tuning values, with their comments, formulas
and JUST_RUN toggles — stays hand-written; only the boilerplate that mirrors it is
generated. Adding a task = add its `mem_/time_` inputs by hand + one RT row here.

Each row records only how a task DEVIATES from the standard Runtime literal:

    Runtime <name> = {
        "docker": <docker>,
        "jar_override": gatk_override,            # only if jar=True
        "preemptible": <preemptible>,
        "max_retries": <max_retries>,
        "cpu": <cpu>,
        "machine_mem": <mem> + <overhead>,
        "command_mem": <mem>,
        "runtime_minutes": <time>,
        "disk": <disk>,
        "boot_disk_size": boot_disk_size
    }

Defaults: cpu="cpu", preemptible="preemptible", max_retries="max_retries",
mem="mem_<name>", overhead="mem_machine_overhead", time="time_startup + time_<name>",
disk="disk". `mem_derive`, if set, emits `Int <mem> = <mem_derive>` before the block
(for the num_bams-scaling tasks). All values are WDL expression strings referencing the
hand-written input variables.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class RT:
    name: str
    docker: str
    jar: bool = False
    cpu: str = "cpu"
    preemptible: str = "preemptible"
    max_retries: str = "max_retries"
    mem: str = ""          # command_mem token; defaults to "mem_<name>"
    mem_derive: str = ""   # if set, emit `Int <mem> = <mem_derive>` before the block
    overhead: str = "mem_machine_overhead"
    time: str = ""         # runtime_minutes RHS; defaults to "time_startup + time_<name>"
    disk: str = "disk"

    def command_mem(self):
        return self.mem or f"mem_{self.name}"

    def runtime_minutes(self):
        return self.time or f"time_startup + time_{self.name}"


# In collection-map order (also the struct field order).
RUNTIMES = [
    # --- Preprocessing ---
    RT("parse_input", "gatk_docker"),
    RT("index_feature_file", "gatk_docker", jar=True),
    RT("get_tumor_sample_names", "ubuntu_docker", jar=True),
    RT("get_sample_name", "gatk_docker", jar=True),
    RT("annotate_intervals", "gatk_docker", jar=True),
    RT("preprocess_intervals", "gatk_docker", jar=True),
    RT("split_intervals", "gatk_docker", jar=True),
    RT("reorder_sam", "gatk_docker", jar=True),
    # --- CNV ---
    RT("collect_callable_loci", "gatk_docker", jar=True, mem="mem_callable_loci",
       time="time_startup + time_callable_loci"),
    RT("collect_read_counts", "gatk_docker", jar=True),
    RT("denoise_read_counts", "gatk_docker", jar=True),
    RT("vcf_to_pileup_variants", "bcftools_docker"),
    RT("get_pileup_summaries", "gatk_docker", jar=True,
       time="time_startup + ceil(time_get_pileup_summaries / scatter_count_for_pileups)"),
    RT("gather_pileup_summaries", "gatk_docker", jar=True),
    RT("select_pileup_summaries", "gatk_docker"),
    RT("pileup_to_allelic_counts", "python_docker"),
    RT("harmonize_copy_ratios", "python_docker", cpu="cpu_harmonize_copy_ratios",
       mem_derive="mem_harmonize_copy_ratios_base + num_bams * mem_harmonize_copy_ratios_additional_per_sample"),
    RT("merge_allelic_counts", "python_docker"),
    RT("calculate_contamination", "gatk_docker", jar=True),
    RT("genotype_variants", "python_docker"),
    RT("model_segments", "gatk_docker", jar=True,
       mem_derive="mem_model_segments_base + num_bams * mem_model_segments_additional_per_sample"),
    RT("call_copy_ratio_segments", "gatk_docker", jar=True),
    RT("plot_modeled_segments", "gatk_docker", jar=True),
    RT("filter_copy_ratios", "python_docker"),
    RT("recount_markers", "python_docker"),
    # --- Clonal analysis ---
    RT("model_segments_to_acs_conversion", "python_docker"),
    RT("process_maf_for_absolute", "python_docker"),
    RT("absolute", "absolute_docker"),
    RT("absolute_extract", "absolute_docker"),
    RT("absolute_extract_postprocess", "python_docker"),
    RT("phylogicndt_task", "phylogicndt_docker"),
    # --- SNV ---
    RT("subset_bam_to_shard", "gatk_docker", jar=True),
    RT("mutect1", "mutect1_docker", cpu="cpu_mutect1", preemptible="preemptible_mutect1",
       max_retries="max_retries_mutect1", mem="mem_mutect1_base", overhead="mem_mutect1_overhead",
       time="time_startup + ceil(time_mutect1_total / scatter_count_for_variant_calling)"),
    RT("merge_mutect1_forcecall_vcfs", "gatk_docker"),
    RT("mutect2", "gatk_docker", jar=True, cpu="cpu_mutect2", preemptible="preemptible_mutect2",
       max_retries="max_retries_mutect2",
       mem_derive="mem_mutect2_base + num_bams * mem_mutect2_additional_per_sample",
       overhead="mem_mutect2_overhead",
       time="time_startup + ceil(time_mutect2_total / scatter_count_for_variant_calling)"),
    RT("learn_read_orientation_model", "gatk_docker", jar=True,
       mem_derive="mem_learn_read_orientation_model_base + num_bams * mem_learn_read_orientation_model_additional_per_sample"),
    RT("gather_vcfs", "gatk_docker", jar=True, mem="mem_merge_vcfs",
       time="time_startup + time_merge_vcfs", disk="1 + disk"),
    RT("merge_mafs", "ubuntu_docker"),
    RT("merge_mutect_stats", "gatk_docker", jar=True),
    RT("print_reads", "gatk_docker", jar=True),
    RT("filter_variant_calls", "gatk_docker", jar=True),
    RT("filter_alignment_artifacts", "gatk_docker", jar=True, cpu="cpu_filter_alignment_artifacts",
       mem_derive="mem_filter_alignment_artifacts_base + num_bams * mem_filter_alignment_artifacts_additional_per_sample",
       time="time_startup + ceil(time_filter_alignment_artifacts_total / scatter_count_for_variant_calling)"),
    RT("select_variants", "gatk_docker", jar=True),
    RT("funcotate", "gatk_docker", jar=True),
    RT("create_empty_annotation", "ubuntu_docker", mem="mem_create_empty_annotations",
       time="time_startup + time_create_empty_annotations", disk="disk_sizeGB"),
    # --- Assorted ---
    RT("create_cnv_panel", "gatk_docker", jar=True, disk="disk + disk_create_cnv_panel"),
    RT("create_mutect2_panel", "gatk_docker", jar=True, overhead="mem_create_mutect2_panel_overhead",
       disk="disk + disk_create_mutect2_panel"),
    RT("select_af_only_from_vcf", "gatk_docker"),
    # --- Ancestry ---
    RT("ancestry", "ancestry_docker"),
]
