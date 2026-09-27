#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Prepare an external binset so refinement and gap filling can resume.

The historical ``-r`` route called outlier removal and ignored ``--module``.
Later stages still expect the checkpoint and list files produced by a normal
autobinning run. This module creates those files from data-feeding outputs
without rerunning candidate generation.
"""

import os
import shutil

from basalt_runtime import append_checkpoint, read_checkpoint_step


OUTLIER_ONLY_NOTICE = (
    "External binset route: --module was not set, so only outlier screening "
    "ran. Contig retrieval, secondary dereplication, rOLC, and reassembly "
    "are available for both CheckM2 and legacy CheckM. From this directory, "
    "rerun the same -r/-a/-c/-s command with --module refinement, "
    "--module reassembly, or --module all. --module reassembly first "
    "finishes any missing refinement stages, then runs gap filling. "
    "BASALT-Air is a separate repository and is not changed by this command."
)

DEREPLICATION_NOTICE = (
    "The -b route performs cross-assembly dereplication only. It does not "
    "run rOLC or reassembly. After BestBinset/ and "
    "BestBinset_comparison_files/ exist, continue with "
    "-r BestBinset --module all (or --module reassembly)."
)


class ExternalBinsetError(ValueError):
    """Raised when an external binset cannot enter refinement or gap filling."""


def external_module_plan(module_explicit, functional_module):
    """Decide which external-binset stages an explicit ``--module`` requests."""
    if functional_module == "autobinning":
        raise ValueError(
            "--module autobinning cannot be combined with --refinement-binset; "
            "the bins already exist. Use --module refinement, reassembly, or all"
        )
    if not module_explicit:
        return "outlier"
    if functional_module == "refinement":
        return "refinement"
    if functional_module in ("reassembly", "all"):
        return "gap-filling"
    raise ValueError("unsupported --module value: {}".format(functional_module))


def connection_stem(path):
    """Return the assembly key encoded in a condensed-connection filename."""
    name = os.path.basename(path)
    prefix = "condense_connections_"
    if name.startswith(prefix) and name.endswith(".txt"):
        return name[len(prefix):-4]
    return ""


def _leading_index(name):
    token = os.path.basename(name).split("_", 1)[0]
    stem = os.path.splitext(token)[0]
    return stem if stem.isdigit() else ""


def pair_connections(assemblies, connection_paths):
    """Align connection files to assemblies.

    Data feeding names connections after the imported binset folder
    (``condense_connections_500_binner_A.txt``) and assemblies as
    ``500_assembly.fa``. Pair by stem, then by the shared numeric prefix.
    Positional pairing is used only when no name-based pair exists.
    """
    unused = list(connection_paths)
    paired = []
    for assembly in assemblies:
        stem = os.path.splitext(os.path.basename(assembly))[0]
        index = _leading_index(assembly)
        match = None
        for path in unused:
            candidate = connection_stem(path)
            if candidate and (candidate == stem or stem.startswith(candidate) or candidate.startswith(stem)):
                match = path
                break
        if match is None and index:
            indexed = [
                path for path in unused
                if connection_stem(path).startswith(index + "_")
                or _leading_index(connection_stem(path)) == index
            ]
            if len(indexed) == 1:
                match = indexed[0]
        paired.append(match)
        if match is not None:
            unused.remove(match)
    if all(item is None for item in paired) and len(assemblies) == len(connection_paths):
        return list(connection_paths)
    return paired


def _local_name(path, workdir, label):
    candidate = path if os.path.isabs(path) else os.path.join(workdir, path)
    if not os.path.isfile(candidate):
        raise ExternalBinsetError("Missing {} file: {}".format(label, path))
    if os.path.dirname(os.path.abspath(candidate)) != os.path.abspath(workdir):
        raise ExternalBinsetError(
            "{} must be in the working directory ({}); got {}. "
            "Run this route from the data-feeding output directory or link "
            "the file there. Absolute paths are not reliable in the legacy "
            "shell stages.".format(label, workdir, path)
        )
    return os.path.basename(candidate)


def _write_lines(path, lines, overwrite):
    if (not overwrite) and os.path.isfile(path) and os.path.getsize(path) > 0:
        return
    with open(path, "w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(str(line).rstrip("\n") + "\n")


def ensure_bestbinset(source, workdir):
    """Make the orchestrator's hard-coded ``BestBinset`` name available.

    Copying a folder onto itself is avoided. An existing ``BestBinset`` is
    kept so ``--mode continue`` does not replace a resumed binset.
    """
    target = os.path.join(workdir, "BestBinset")
    source_path = source if os.path.isabs(source) else os.path.join(workdir, source)
    if os.path.abspath(source_path) == os.path.abspath(target):
        if not os.path.isdir(target):
            raise ExternalBinsetError("Refinement binset does not exist: " + source)
        return target
    if os.path.isdir(target):
        return target
    if not os.path.isdir(source_path):
        raise ExternalBinsetError("Refinement binset does not exist: " + source)
    shutil.copytree(source_path, target)
    return target


def prepare_external_binset_route(refinement_binset, assemblies, coverages,
                                  functional_module, running_mode, workdir=None):
    """Write the state files used by steps 4-12 and return the orchestrator module.

    ``--module reassembly`` is mapped to the full post-binning orchestrator.
    Unfinished refinement stages are prerequisites for OLC and reassembly;
    a checkpoint at step 7 or later skips them.
    """
    workdir = workdir or os.getcwd()
    plan = external_module_plan(True, functional_module)
    if plan == "outlier":
        raise ExternalBinsetError("outlier-only runs do not need pipeline state")

    assembly_names = [_local_name(path, workdir, "assembly") for path in assemblies]
    coverage_names = [_local_name(path, workdir, "coverage") for path in coverages]
    if len(assembly_names) != len(coverage_names):
        raise ExternalBinsetError(
            "--assemblies and --coverage-list must contain the same number of entries"
        )
    ensure_bestbinset(refinement_binset, workdir)

    connection_paths = sorted(
        os.path.join(workdir, name)
        for name in os.listdir(workdir)
        if connection_stem(name)
    )
    paired = pair_connections(assembly_names, connection_paths)
    missing_connections = [
        assembly for assembly, connection in zip(assembly_names, paired)
        if not connection or not os.path.isfile(connection)
    ]
    if missing_connections:
        raise ExternalBinsetError(
            "Contig retrieval needs condense_connections_*.txt from data feeding "
            "for: {}. Run -d first and execute -r from that output directory. "
            "Found connection files: {}.".format(
                ", ".join(missing_connections),
                ", ".join(os.path.basename(path) for path in connection_paths) or "none",
            )
        )

    comparison_dir = os.path.join(workdir, "BestBinset_comparison_files")
    if plan == "gap-filling" and len(assembly_names) > 1 and not os.path.isdir(comparison_dir):
        raise ExternalBinsetError(
            "Gap filling needs BestBinset_comparison_files from the -b "
            "dereplication step. Run -b on the data-fed binsets first, then "
            "rerun -r BestBinset --module all from that directory."
        )

    overwrite = running_mode == "new"
    _write_lines(
        os.path.join(workdir, "Coverage_matrix_list.txt"),
        coverage_names,
        overwrite,
    )
    _write_lines(
        os.path.join(workdir, "Assembly_mo_list.txt"),
        assembly_names,
        overwrite,
    )
    _write_lines(
        os.path.join(workdir, "Bestbinset_list.txt"),
        ["BestBinset"] * len(assembly_names),
        overwrite,
    )
    _write_lines(
        os.path.join(workdir, "Assembly_MoDict.txt"),
        [
            "{}\t{}".format(connection_stem(connection), assembly)
            for assembly, connection in zip(assembly_names, paired)
        ],
        overwrite,
    )

    checkpoint = os.path.join(workdir, "Basalt_checkpoint.txt")
    if running_mode == "new":
        with open(checkpoint, "w", encoding="utf-8") as handle:
            handle.write("3rd external binset supplied; autobinning skipped\n")
        step = 3
    else:
        step = read_checkpoint_step(checkpoint)
        if step < 3:
            append_checkpoint(
                checkpoint,
                "3rd external binset supplied; autobinning skipped",
            )
            step = 3
        outlier_dir = os.path.join(workdir, "BestBinset_outlier_refined")
        if step < 4 and os.path.isdir(outlier_dir):
            append_checkpoint(checkpoint, "4th outlier removal done!")
            step = 4

    return {
        "plan": plan,
        "checkpoint_step": step,
        "connections": [os.path.basename(path) for path in paired],
        "orchestrator_module": "all" if plan == "gap-filling" else "refinement",
        "comparison_dir": comparison_dir,
    }
