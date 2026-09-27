#!/usr/bin/env python3

"""Small runtime checks shared by BASALT command modules."""

from __future__ import annotations

import fcntl
import hashlib
import os
import re
import shutil
import subprocess
from pathlib import Path


EXPECTED_MODEL_ENSEMBLES = 5
_CHECKPOINT_STEP = re.compile(r"^(\d+)")
_CHECKM2_METRIC_KEYS = ("N50", "Completeness", "Genome size", "Contamination")


def require_model_directory() -> Path:
    """Return a validated BASALT model directory or raise an actionable error."""
    configured = os.environ.get("BASALT_WEIGHT")
    if not configured:
        raise RuntimeError(
            "BASALT_WEIGHT is not set. Download the five BASALT model ensembles "
            "with BASALT_models_download.py, then export BASALT_WEIGHT to that "
            "absolute directory."
        )

    model_dir = Path(configured).expanduser().resolve()
    if not model_dir.is_dir():
        raise RuntimeError(
            "BASALT_WEIGHT does not name a directory: {}".format(model_dir)
        )

    descriptors = sorted(model_dir.glob("*_ensemble.csv"))
    if len(descriptors) != EXPECTED_MODEL_ENSEMBLES:
        raise RuntimeError(
            "BASALT_WEIGHT must contain {} top-level *_ensemble.csv files; "
            "found {} in {}.".format(
                EXPECTED_MODEL_ENSEMBLES, len(descriptors), model_dir
            )
        )

    missing_directories = [
        descriptor.with_suffix("").name
        for descriptor in descriptors
        if not descriptor.with_suffix("").is_dir()
    ]
    if missing_directories:
        raise RuntimeError(
            "BASALT_WEIGHT is missing checkpoint directories: {}".format(
                ", ".join(missing_directories)
            )
        )

    missing_checkpoints = []
    unsafe_checkpoints = []
    for descriptor in descriptors:
        checkpoint_dir = descriptor.with_suffix("").resolve()
        entries = [
            line.strip()
            for line in descriptor.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        if not entries:
            missing_checkpoints.append("{} (empty descriptor)".format(descriptor.name))
            continue
        for entry in entries:
            checkpoint = (checkpoint_dir / entry).resolve()
            if checkpoint_dir not in checkpoint.parents:
                unsafe_checkpoints.append("{}/{}".format(checkpoint_dir.name, entry))
            elif not checkpoint.is_file():
                missing_checkpoints.append("{}/{}".format(checkpoint_dir.name, entry))
    if unsafe_checkpoints:
        raise RuntimeError(
            "BASALT_WEIGHT descriptors contain unsafe checkpoint paths: {}".format(
                ", ".join(unsafe_checkpoints[:10])
            )
        )
    if missing_checkpoints:
        raise RuntimeError(
            "BASALT_WEIGHT has missing checkpoint files: {}".format(
                ", ".join(missing_checkpoints[:10])
                + (" ..." if len(missing_checkpoints) > 10 else "")
            )
        )

    return model_dir


def read_checkpoint_step(path="Basalt_checkpoint.txt"):
    """Return the highest numbered checkpoint step without modifying the file.

    Non-numeric lines, including the historical ``Outlier removal done!``
    message, are ignored. A missing or unreadable file means step 0. This
    deliberately does not truncate the checkpoint: a failed parse used to
    reset a resumable run to the beginning.
    """
    try:
        lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    except (OSError, ValueError):
        return 0
    step = 0
    for line in lines:
        match = _CHECKPOINT_STEP.match(line.strip())
        if match:
            step = max(step, int(match.group(1)))
    return step


def append_checkpoint(path, line):
    """Append one checkpoint marker, creating the file if needed."""
    checkpoint = Path(path)
    existing = ""
    if checkpoint.exists():
        existing = checkpoint.read_text(encoding="utf-8", errors="replace")
    prefix = ""
    if existing and not existing.endswith("\n"):
        prefix = "\n"
    elif existing:
        prefix = "\n"
    with checkpoint.open("a", encoding="utf-8") as handle:
        handle.write(prefix + str(line).rstrip("\n") + "\n")


def depth_sample_count(header_line):
    """Count samples in a ``jgi_summarize_bam_contig_depths`` header.

    MetaBAT 2.15+ can omit ``.bam`` from sample column names, so counting
    the literal ``.bam-var`` substring reports zero samples and later divides
    by zero. Variance columns still end in ``-var``.
    """
    fields = str(header_line).strip().split("\t")
    variance_columns = [field for field in fields if field.endswith("-var")]
    if variance_columns:
        return len(variance_columns)
    extra = len(fields) - 3
    if extra > 0 and extra % 2 == 0:
        return extra // 2
    raise ValueError(
        "Could not determine the sample count from a depth header. "
        "Expected MetaBAT columns ending in -var, got: {}".format(
            "\t".join(fields[:8])
        )
    )


def coerce_count(value):
    """Parse a count that may be an integer or a float-formatted string."""
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError) as exc:
        raise ValueError("expected a numeric count, got {!r}".format(value)) from exc


def positive_thread_count(value, default=1):
    """Return a positive thread count, falling back when the value is unusable."""
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed > 0 else default


def worker_thread_budget(total_threads, worker_count):
    """Split a thread allocation across concurrent OLC/BLAST workers.

    Giving every worker the full ``-t`` value while a process pool is also
    sized to that value oversubscribes the node. One thread per worker is
    the other extreme and leaves a serial bin-pair loop idle.
    """
    total = positive_thread_count(total_threads, 1)
    workers = max(1, int(worker_count))
    return max(1, total // workers)


def bundled_executable(script_name, caller_file):
    """Resolve a Perl helper next to the caller, then on ``PATH``.

    ``perl Cytoscapeviz.pl`` searches the working directory, not ``PATH``,
    so an installed copy is invisible unless the command is given its path.
    """
    local = Path(caller_file).resolve().parent / script_name
    if local.is_file():
        return str(local)
    found = shutil.which(script_name)
    if found:
        return found
    raise FileNotFoundError(
        "{} was not found next to {} or on PATH. Re-run install.sh from "
        "this BASALT checkout.".format(script_name, caller_file)
    )


def pe_read_names(mates):
    """Return idempotent ``PE_r1_/PE_r2_`` read names.

    Internal stages historically build ``PE_r1_`` + the original filename.
    The external-binset route receives the already-prefixed data-feeding
    files, which would double the prefix and break reassembly mapping.
    """
    r1, r2 = str(mates[0]), str(mates[1])
    if not r1.startswith("PE_r1_"):
        r1 = "PE_r1_" + r1
    if not r2.startswith("PE_r2_"):
        r2 = "PE_r2_" + r2
    return [r1, r2]


def fasta_has_sequence(path):
    """Return whether a FASTA file contains at least one non-empty sequence."""
    try:
        handle = open(path, "r", encoding="utf-8", errors="replace")
    except OSError:
        return False
    with handle:
        seen_header = False
        for line in handle:
            if line.startswith(">"):
                seen_header = True
            elif seen_header and line.strip():
                return True
    return False


def expand_alignment_group(group, alignment_dict):
    """Grow one BLAST connected component and always stop.

    The historical ``while num2 != num1`` expansion only adds identifiers,
    so it normally converges. A bound remains so a resumed or malformed
    alignment table cannot spin until the scheduler kills the job.
    """
    limit = max(2, len(alignment_dict) + 2)
    previous = -1
    for _ in range(limit):
        current = len(group)
        if current == previous:
            return group, True
        previous = current
        for alignment in alignment_dict:
            parts = alignment.split(" ")
            if len(parts) < 2:
                continue
            query_id, subject_id = parts[0], parts[1]
            rendered = str(group)
            if ("'" + query_id + "'") in rendered or ("'" + subject_id + "'") in rendered:
                group[query_id] = 1
                group[subject_id] = 1
    return group, len(group) == previous


def clean_stale_olc_intermediates(work_dir, target_bin):
    """Remove partial products for one bin before OLC elongation resumes.

    A wall-time kill can leave ``*_merged``, ``blast_*self_merged*``, and a
    per-bin ``*_checkm`` directory. Reusing them makes step 8 re-evaluate
    the same candidate set. Completed ``*_OLC`` result folders are not
    touched, and only the named bin is cleaned so parallel workers do not
    delete each other's files.
    """
    work = Path(work_dir)
    name = Path(str(target_bin)).name
    for suffix in ("_merged", "_checkm", "_split_blast_output"):
        path = work / (name + suffix)
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)
        elif path.exists():
            path.unlink()
    for pattern in (
        "blast_" + name + "_self_merged_*",
        "Merged_seqs_" + name + "_*",
        "Merged_*_" + name,
    ):
        for path in work.glob(pattern):
            if path.is_dir():
                shutil.rmtree(path, ignore_errors=True)
            elif path.is_file():
                path.unlink()


def _run(command, cwd=None):
    """Run an external command, raising RuntimeError on a non-zero exit."""
    rendered = [str(item) for item in command]
    completed = subprocess.run(
        rendered,
        cwd=None if cwd is None else str(cwd),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "command failed ({}): {}".format(completed.returncode, " ".join(rendered))
        )


def fasta_content_hash(path):
    """Return an MD5 over uppercase sequence content, ignoring headers.

    CheckM2 metrics depend only on sequence content, so bins that reappear
    unchanged across OLC convergence iterations share one cache entry even
    when their filenames or headers differ.
    """
    digest = hashlib.md5()
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            stripped = "".join(line.split()).upper()
            if stripped:
                digest.update(stripped.encode("ascii", "ignore"))
    return digest.hexdigest()


def _extract_version_token(text):
    """Return the CheckM2 version token from possibly banner-prefixed output."""
    text = (text or "").strip()
    match = re.search(r"version\s+v?([0-9][0-9A-Za-z.+-]*)", text, re.IGNORECASE)
    if match:
        return match.group(1)
    for line in reversed(text.splitlines()):
        candidate = line.strip()
        if re.fullmatch(r"v?[0-9][0-9A-Za-z.+-]*", candidate):
            return candidate.lstrip("v")
    return None


def _checkm2_signature():
    """Describe the CheckM2 build and database for cache invalidation.

    ``checkm2 --version`` emits a timestamped banner and may prepend warning
    lines, so only the version token is kept; taking the whole output would
    make the signature (and the cache) unstable.
    """
    parts = ["checkm2"]
    try:
        version = subprocess.run(
            ["checkm2", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=120,
        )
        token = _extract_version_token(version.stdout)
        parts.append(token if token else "noversion")
    except Exception:  # noqa: BLE001 - signature is best effort
        parts.append("unknown")
    database = os.environ.get("CHECKM2DB")
    if database and os.path.isfile(database):
        stat = os.stat(database)
        parts.append("{}:{}:{}".format(database, stat.st_size, int(stat.st_mtime)))
    return "|".join(parts)


def _load_checkm2_cache(path, signature):
    cache = {}
    target = Path(path)
    if not target.is_file():
        return cache
    try:
        lines = target.read_text(encoding="utf-8").splitlines()
    except OSError:
        return cache
    if not lines or not lines[0].startswith("#") or signature not in lines[0]:
        return cache
    for line in lines[1:]:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 5:
            continue
        digest, n50, completeness, genome_size, contamination = fields
        try:
            cache[digest] = {
                "N50": int(float(n50)),
                "Completeness": float(completeness),
                "Genome size": int(float(genome_size)),
                "Contamination": float(contamination),
            }
        except ValueError:
            continue
    return cache


def _write_checkm2_cache_header(path, signature):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("# {}\n".format(signature))


def _append_checkm2_cache(path, entries):
    """Append evaluation results under an exclusive lock.

    Concurrent OLC workers share one cache file. Duplicate entries are
    harmless because loading keeps the last occurrence per digest.
    """
    with open(path, "a", encoding="utf-8") as handle:
        locked = False
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            locked = True
        except (AttributeError, OSError):
            locked = False
        try:
            for digest, metrics in entries:
                handle.write(
                    "{}\t{}\t{}\t{}\t{}\n".format(
                        digest,
                        metrics["N50"],
                        metrics["Completeness"],
                        metrics["Genome size"],
                        metrics["Contamination"],
                    )
                )
            handle.flush()
        finally:
            if locked:
                try:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
                except (AttributeError, OSError):
                    pass


def _parse_quality_report(report_path):
    metrics = {}
    with open(report_path, "r", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            if index == 0:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                if len(fields) >= 5:
                    genome_size, completeness, contamination, n50 = (
                        fields[1], fields[2], fields[3], fields[4],
                    )
                else:
                    continue
            else:
                completeness, contamination, n50, genome_size = (
                    fields[1], fields[2], fields[6], fields[8],
                )
            try:
                metrics[fields[0].strip()] = {
                    "N50": int(float(n50)),
                    "Completeness": float(completeness),
                    "Genome size": coerce_count(genome_size),
                    "Contamination": float(contamination),
                }
            except ValueError:
                continue
    return metrics


_EVAL_COUNTER = [0]


def evaluate_bins_checkm2(folder, num_threads, cache_path="Checkm2_metrics_cache.tsv",
                          notify=print):
    """Evaluate candidate bin FASTAs with CheckM2, reusing identical content.

    OLC convergence loops re-evaluate mostly unchanged bins: the original
    target bin is copied into every ``*_merged`` batch and repeated merge
    candidates recur across iterations. This helper keys deterministic
    CheckM2 metrics by sequence-content MD5 so unchanged candidates cost
    nothing, and evaluates only new content in one batched multi-threaded
    ``checkm2 predict`` call.

    Parameters
    ----------
    folder : str
        Directory of candidate FASTA files (``.fa``/``.fasta``).
    num_threads : int
        Thread budget for the batched prediction.
    cache_path : str
        Append-only TSV cache shared by concurrent workers.
    notify : callable
        Progress reporter used for cache-hit summaries.

    Returns
    -------
    dict mapping ``bin_stem -> metrics`` with the same keys as the
    historical parsers, or ``None`` when the batched prediction failed so
    the caller can fall back to the direct per-folder invocation.
    """
    source = Path(folder)
    if not source.is_dir():
        return None
    candidates = sorted(
        path for path in source.iterdir()
        if path.is_file() and path.suffix in (".fa", ".fasta")
    )
    if not candidates:
        return {}

    signature = _checkm2_signature()
    cache = _load_checkm2_cache(cache_path, signature)
    if not Path(cache_path).exists():
        _write_checkm2_cache_header(cache_path, signature)

    results, fresh = {}, []
    for path in candidates:
        stem = path.name[: -len(path.suffix)]
        digest = fasta_content_hash(path)
        if digest in cache:
            results[stem] = dict(cache[digest])
        else:
            fresh.append((path, stem, digest))

    if not fresh:
        notify(
            "CheckM2 cache: {}/{} candidate(s) reused, 0 evaluated".format(
                len(results), len(candidates)
            )
        )
        return results

    _EVAL_COUNTER[0] += 1
    work = Path(
        "Checkm2_batch_eval_{}_{}".format(os.getpid(), _EVAL_COUNTER[0])
    )
    output = Path(str(work) + "_out")
    shutil.rmtree(work, ignore_errors=True)
    shutil.rmtree(output, ignore_errors=True)
    try:
        work.mkdir(parents=True)
        for path, _stem, _digest in fresh:
            target = work / path.name
            try:
                os.link(path, target)
            except OSError:
                shutil.copy2(path, target)
        _run(
            [
                "checkm2", "predict",
                "-t", positive_thread_count(num_threads),
                "-i", str(work),
                "-x", "fa",
                "-o", str(output),
                "--force",
            ]
        )
        report = output / "quality_report.tsv"
        evaluated = _parse_quality_report(report)
        appended = []
        for path, stem, digest in fresh:
            if stem in evaluated:
                results[stem] = evaluated[stem]
                appended.append((digest, evaluated[stem]))
        if appended:
            _append_checkm2_cache(cache_path, appended)
        missing = len(fresh) - len(appended)
        notify(
            "CheckM2 cache: {}/{} candidate(s) reused, {} evaluated".format(
                len(results) - len(appended), len(candidates), len(appended)
            )
            + (", {} without a report".format(missing) if missing else "")
        )
        return results
    except Exception as exc:  # noqa: BLE001 - fall back to the direct call
        notify("CheckM2 batch evaluation failed ({}); using direct prediction".format(exc))
        return None
    finally:
        shutil.rmtree(work, ignore_errors=True)
        shutil.rmtree(output, ignore_errors=True)
