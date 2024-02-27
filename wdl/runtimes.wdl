version development

# A NOTE ON MULTI-THREADING:
# Asking for a machine on gcp with specified cpu and memory will make gcp provide
# a machine with AT LEAST that much cpu and memory. E.g. a request for a machine
# with 4 cpus and 64GB mem can result in a delivered machine with 10 cpus. Thus,
# tasks should not necessarily request a specific number of threads, but use
# $(nproc) instead. That said, tasks with more threads might need higher memory.


struct Runtime {
    String docker
    File? jar_override
    Int preemptible
    Int max_retries
    Int cpu
    Int machine_mem
    Int command_mem
    Int runtime_minutes
    Int disk
    Int boot_disk_size
}


workflow UpdateRuntimeParameters {
    input {
        Runtime runtime_params
        String? docker
        File? jar_override
        Int? preemptible
        Int? max_retries
        Int? cpu
        Int? machine_mem
        Int? command_mem
        Int? runtime_minutes
        Int? disk
        Int? boot_disk_size
    }

    Runtime updated_runtime = {
        "docker": select_first([docker, runtime_params.docker]),
        # cannot use select_first for optional fields:
        "jar_override": if defined(jar_override) then jar_override else runtime_params.jar_override,
        "preemptible": select_first([preemptible, runtime_params.preemptible]),
        "max_retries": select_first([max_retries, runtime_params.max_retries]),
        "cpu": select_first([cpu, runtime_params.cpu]),
        "machine_mem": select_first([machine_mem, runtime_params.machine_mem]),
        "command_mem": select_first([command_mem, runtime_params.command_mem]),
        "runtime_minutes": select_first([runtime_minutes, runtime_params.runtime_minutes]),
        "disk": select_first([disk, runtime_params.disk]),
        "boot_disk_size": select_first([boot_disk_size, runtime_params.boot_disk_size])
    }

    output {
        Runtime params = updated_runtime
    }
}