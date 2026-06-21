#!/usr/bin/env python3
"""Generate WDL struct definitions and bulk single-pass constructors from spec.py.

Sketch / proof-of-concept. Usage:
    python generate.py structs           # emit the SequencingRun/Sample/Patient struct blocks
    python generate.py update_samples     # emit the full UpdateSamples workflow

The derivation rules below ("File?" -> overlay input "Array[File]?", default "[]",
length test "length(...)"; "Array[SequencingRun]" -> "Array[Array[SequencingRun]]?",
default "[[]]", length test "length(flatten(...))") are encoded once here. Adding or
renaming a field means editing spec.py only.
"""

import sys
import spec

I = "    "  # one indent level


def strip_opt(t):
    return t[:-1] if t.endswith("?") else t


def is_array(elem):
    return elem.startswith("Array[")


def overlay_input_type(f):
    return f"Array[{strip_opt(f.wdl_type)}]?"


def arr_decl(f):
    """`Array[<elem>] <name>_arr = select_first([<name>, <default>])`"""
    elem = strip_opt(f.wdl_type)
    default = "[[]]" if is_array(elem) else "[]"
    return f"Array[{elem}] {f.name}_arr = select_first([{f.name}, {default}])"


def length_test(f):
    var = f"{f.name}_arr"
    return f"length(flatten({var}))" if is_array(strip_opt(f.wdl_type)) else f"length({var})"


# ---------------------------------------------------------------- struct blocks

def emit_struct(name):
    fields = spec.STRUCTS[name]
    width = max(len(f.wdl_type) for f in fields)
    lines = [f"struct {name} {{"]
    for f in fields:
        line = f"{I}{f.wdl_type:<{width}} {f.name}"
        if f.origin:
            pad = " " * max(2, 60 - len(line))  # always >= 2 spaces before the comment
            line = f"{line}{pad}# {f.origin}"
        lines.append(line)
    lines.append("}")
    return "\n".join(lines)


# ---------------------------------------------------------------- UpdateSamples

def emit_update_samples():
    sample_over = [f for f in spec.SAMPLE if f.overlaid]
    out = []
    w = out.append

    w("version development\n")
    w('import "sample.wdl" as s')
    w('import "patient.wdl" as p\n\n')

    w("workflow UpdateSamples {")
    w(f"{I}input {{")
    w(f"{I}{I}Patient patient")
    for f in sample_over:
        w(f"{I}{I}{overlay_input_type(f)} {f.name}")
    w(f"{I}}}\n")

    w(f"{I}# Build each updated Sample in a single pass. Each provided overlay array carries")
    w(f"{I}# one entry per sample (aligned by index); an absent or empty array leaves the")
    w(f"{I}# field unchanged. INVARIANT: any provided overlay array has length == length(patient.samples).")
    for f in sample_over:
        w(f"{I}{arr_decl(f)}")
    w("")

    w(f"{I}scatter (i in range(length(patient.samples))) {{")
    w(f"{I}{I}Sample sample_i = patient.samples[i]")
    w(f"{I}{I}Sample updated_sample_i = object {{")
    body = []
    for f in spec.SAMPLE:
        if f.overlaid:
            body.append(f"{I}{I}{I}{f.name}: if {length_test(f)} > 0 then {f.name}_arr[i] else sample_i.{f.name}")
        else:
            body.append(f"{I}{I}{I}{f.name}: sample_i.{f.name}")
    w(",\n".join(body))
    w(f"{I}{I}}}")
    w(f"{I}}}")
    w(f"{I}Array[Sample] samples = updated_sample_i\n")

    w(TUMOR_NORMAL_SELECTION)

    w(f"{I}Patient updated_patient_obj = object {{")
    pbody = []
    for f in spec.PATIENT:
        if f.name == "normal_samples":
            rhs = "select_first([normal_samples, patient.normal_samples])"
        elif f.name == "matched_normal_sample":
            rhs = "if defined(matched_normal_sample) then matched_normal_sample else patient.matched_normal_sample"
        elif f.name in spec.PATIENT_COMPUTED:      # samples, tumor_samples
            rhs = f.name
        else:
            rhs = f"patient.{f.name}"
        pbody.append(f"{I}{I}{f.name}: {rhs}")
    w(",\n".join(pbody))
    w(f"{I}}}\n")

    w(f"{I}output {{")
    w(f"{I}{I}Patient updated_patient = updated_patient_obj")
    w(f"{I}}}")
    w("}")
    return "\n".join(out)


# Fixed (not field-driven) boilerplate: re-select tumor/normal/matched-normal samples
# from the freshly built `samples` by name.
TUMOR_NORMAL_SELECTION = """\
    scatter (tumor_sample in patient.tumor_samples) {
        scatter (sample in samples) {
            if (sample.name == tumor_sample.name) {
                Sample selected_tumor_sample = sample
            }
        }
        Sample tumor_samples = select_all(selected_tumor_sample)[0]
    }

    if (patient.has_normal) {
        scatter (normal_sample in patient.normal_samples) {
            scatter (sample in samples) {
                if (sample.name == normal_sample.name) {
                    Sample selected_normal_sample = sample
                }
            }
            Sample normal_samples = select_all(selected_normal_sample)[0]
        }
        if (defined(patient.matched_normal_sample)) {
            Sample previous_matched_normal_sample = select_first([patient.matched_normal_sample])
            scatter (sample in normal_samples) {
                if (sample.name == previous_matched_normal_sample.name) {
                    Sample matched_normal_samples = sample
                }
            }
            Sample matched_normal_sample = select_all(matched_normal_samples)[0]
        }
    }
"""


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "update_samples"
    if what == "structs":
        print("\n\n".join(emit_struct(n) for n in spec.STRUCTS))
    elif what == "update_samples":
        print(emit_update_samples())
    else:
        sys.exit(f"unknown target: {what}")
