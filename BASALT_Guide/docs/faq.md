# Troubleshooting and FAQ

Start from evidence. Preserve the complete command, stdout, stderr, `Basalt_log.txt`, `Basalt_checkpoint.txt`, software versions, database configuration, and input checksums before changing files or rerunning a stage.

## Fast diagnostic sequence

Activate the environment first, or add `micromamba run -n basalt` before environment commands:

```bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate basalt
# Conda alternative: conda activate basalt

pwd
command -v micromamba >/dev/null && micromamba env list
command -v conda >/dev/null && conda info --envs
command -v BASALT checkm2 metabat2 SemiBin2 bowtie2 samtools spades.py
test -n "$BASALT_WEIGHT" && find "$BASALT_WEIGHT" -maxdepth 1 -name '*_ensemble.csv'
tail -n 30 Basalt_checkpoint.txt
tail -n 100 Basalt_log.txt
grep -Ei 'error|failed|warning|traceback' *.log 2>/dev/null || true
```

## Installation and environment

### Environment solving or package downloads are slow

Use the China-mainland installer in dry-run mode to inspect the manager, channels, PyPI index, and model source before retrying:

```bash
python3 BASALT_setup_China_mainland.py \
  --manager micromamba \
  --bootstrap-micromamba \
  --mirror auto \
  --model-source none \
  --dry-run
```

Remove `--dry-run` when the selected route is correct. The installer prefers micromamba and keeps mirror settings local to its subprocesses. If one route is connected but incomplete, choose `--mirror tuna`, `bfsu`, `ustc`, or `upstream` explicitly and pass `--update` for an existing environment. Do not repeatedly append channels to the global `.condarc`; mixed channel priority is a common cause of slow or inconsistent solves.

### `BASALT: command not found`

Activate the environment used during installation and inspect the launcher:

```bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate basalt
# Conda alternative: conda activate basalt
command -v BASALT
ls -l "$CONDA_PREFIX/bin/BASALT" "$CONDA_PREFIX/bin/BASALT.py"
```

Without activation, verify the launcher with `micromamba run -n basalt BASALT --help`.

If the launcher is absent, return to the cloned repository and reinstall:

```bash
cd /path/to/BASALT
micromamba run -n basalt bash install.sh
# Conda alternative: conda run -n basalt bash install.sh
```

### `Permission denied`, or should I run `chmod 777`?

No. `chmod 777` makes a file or directory writable by every local user, and `chmod -R 777` can expose launchers, Python code, model weights, databases, and results to accidental or malicious replacement. `chmod +x FILE` only adds executable permission; it does not mean `777`.

The repository installer should be run as `bash install.sh`, so the source installer does not need an executable bit. It verifies that the active environment is writable and installs executable files with mode `0755`. Diagnose the target before changing any permission:

```bash
printf 'user: %s\nenvironment: %s\n' "$(id -un)" "$CONDA_PREFIX"
ls -ld "$CONDA_PREFIX" "$CONDA_PREFIX/bin"
ls -l \
  "$CONDA_PREFIX/bin/BASALT" \
  "$CONDA_PREFIX/bin/BASALT.py" \
  "$CONDA_PREFIX/bin/jgi_summarize_bam_contig_depths"
test -w "$CONDA_PREFIX/bin" && echo 'environment bin is writable'
```

If the environment is writable but an older installation has incorrect executable modes, return to the matching BASALT checkout and rerun the installer:

```bash
micromamba run -n basalt bash install.sh
micromamba run -n basalt BASALT --help
```

For Conda, replace `micromamba run` with `conda run`. If the environment belongs to another account or to root, create a new user-owned environment or ask the administrator to configure controlled group access. Do not repair ownership problems with `sudo chmod -R 777`. On clusters, shared software permissions and ACLs should follow the local administrator's policy.

### `LorBin: command not found`

The `-e l` adapter does not install LorBin automatically. It targets the BASALT-compatible [LorBin–BASALT Extra-binner fork](https://github.com/PKU-EMBL/LorBin-BASALT-Extrabinner), not an arbitrary upstream installation. From the cloned fork, install it into the same environment that launches BASALT and retain its source commit:

```bash
git -C /path/to/LorBin-BASALT-Extrabinner rev-parse HEAD
micromamba run -n basalt python -m pip install --no-deps \
  /path/to/LorBin-BASALT-Extrabinner
micromamba run -n basalt LorBin --help
```

Conda users can replace `micromamba run` with `conda run`. If the command exists only in a separate environment, BASALT will not discover it automatically. See [Extra binners](extra-binners.md#lorbin) for dependency checks, the standalone-environment fallback, and output validation.

### BASALT models are missing

The pipeline expects `BASALT_WEIGHT` to name a directory containing five ensemble descriptor files and the associated checkpoints.

```bash
export BASALT_WEIGHT=/absolute/path/to/BASALT_WEIGHT
find "$BASALT_WEIGHT" -maxdepth 1 -name '*_ensemble.csv' | wc -l
```

If the count is not `5`, rerun the model downloader and inspect extraction errors:

```bash
micromamba run -n basalt \
  BASALT_models_download.py --source auto --path "$BASALT_WEIGHT"
```

For a persistent correction, keep one validated absolute path in `~/.bashrc`:

```bash
grep -n 'BASALT_WEIGHT' ~/.bashrc
# Edit or remove stale entries, then retain one line such as:
export BASALT_WEIGHT="/absolute/persistent/path/BASALT_WEIGHT"
```

Run `source ~/.bashrc` in the current interactive shell. For a scheduler that does not read `~/.bashrc`, export the same path in the submission script.

Automatic mode tries the official Hugging Face repository before Figshare. For a manually obtained Baidu Netdisk ZIP or another trusted local copy:

```bash
micromamba run -n basalt BASALT_models_download.py \
  --source archive \
  --archive /absolute/path/BASALT.zip \
  --path "$BASALT_WEIGHT"
```

### `CheckM2 database not found`

```bash
micromamba run -n basalt checkm2 database --download
micromamba run -n basalt checkm2 database --current
export CHECKM2DB=/absolute/path/to/CheckM2_database/uniref100.KO.1.dmnd
```

Only export `CHECKM2DB` when automatic discovery fails. If it must persist, keep one validated absolute definition in `~/.bashrc` and repeat it explicitly in scheduler jobs that do not source the file.

Test the installed CheckM2 version with a small independent input before rerunning a long BASALT analysis.

### `libcrypto.so` or another shared library is missing

Do not create an ABI-changing symbolic link as a first response. That can make a program start while remaining binary-incompatible.

Inspect the executable and solve the environment consistently:

```bash
command -v samtools
ldd "$(command -v samtools)" | grep 'not found' || true
conda list samtools openssl
```

Reinstall compatible packages in a clean environment or use the project container.

## Inputs and paths

### `IndexError: list index out of range` near input parsing

Check the paired-end grammar:

```text
one sample:       R1.fastq,R2.fastq
two samples:      S1_R1.fastq,S1_R2.fastq/S2_R1.fastq,S2_R2.fastq
```

Also confirm that filenames contain no commas or `/` characters beyond the intended delimiters.

### Do absolute paths work?

Absolute paths are not handled consistently across all legacy shell-command paths. Use a dedicated working directory and simple local symlinks. BASALT-Air is the separate implementation intended for first-class absolute-path and separate work/output-directory support.

### Can compressed reads be used?

The pipeline contains handlers for `.gz`, `.tar.gz`, and `.zip`. Plain FASTA and FASTQ inputs remain easier to audit. Test the exact archive layout on a small run because automated extraction can introduce unexpected filenames or overwrite existing files.

### Can I combine single assemblies and co-assemblies?

Yes, provided the supplied read coverage is biologically compatible with each assembly. BASALT compares selected candidates across assemblies to reduce redundancy. Document which samples contributed to every assembly and coverage matrix.

## Runtime and checkpoints

### A scheduler reports success, but the output is incomplete

Some external tools are invoked through shell calls, so a failed subcommand may not always produce a non-zero final BASALT status. Audit:

1. the last checkpoint;
2. errors in captured stderr;
3. optional-binner warnings;
4. the final FASTA count and sizes;
5. the final quality report.

### How do I resume?

From the same working directory and unchanged environment:

```bash
BASALT --mode continue
```

Checkpoint parsing ignores non-numeric lines, so a trailing `Outlier removal done!` message from an external-binset screening run no longer resets progress to step 0.

Do not resume after changing inputs, code, model files, databases, or major dependency versions. Start a new run directory instead.

### Can gap filling (rOLC and reassembly) run on an external binset fed through `-d`/`-b`/`-r`?

Yes. The three-stage route is: `-d` for data feeding, `-b` for dereplication, then `-r BestBinset --module all` (or `--module reassembly`) from the same working directory. Without `--module`, `-r` still performs outlier screening only and prints the follow-up command. `--module reassembly` finishes any refinement stages the checkpoint has not reached before OLC and reassembly. The same dispatch is used for CheckM2 (`-q checkm2`, default) and legacy CheckM. BASALT-Air is a separate codebase; check its own repository for its current behaviour.

Requirements for the `-r --module` route:

- run from the data-feeding output directory with `condense_connections_*.txt`, coverage matrices, and `PE_r1_*`/`PE_r2_*` reads present;
- for multiple assemblies, complete `-b` first so `BestBinset_comparison_files/` exists;
- unchanged contig identifiers and matching list order across `-a`, `-b`, `-c`, and `-r`.

### A wall-time kill interrupted the OLC / bin-merging stage (step 8). Is the state safe to resume?

The OLC workers now discard partial per-bin products (`*_merged`, `blast_*self_merged*`, per-bin `*_checkm`, `Merged_seqs_*`) when a bin restarts, and an empty FASTA in a `*_OLC` result folder is no longer counted as finished work, so `--mode continue` rebuilds instead of re-evaluating the same candidate set indefinitely. If a legacy run still loops after resume, remove the stale step-8 intermediates by hand and resume again.

### The OLC stage uses only a fraction of my cores

Stage 8 and stage 10 now size their worker pool from `-t` and give each worker a share of the thread budget instead of issuing single-threaded `blastn` and one-thread CheckM2 calls. Report the worker/thread split printed at stage start (`OLC workers: N; threads per worker: M`) together with `-t` when comparing performance between releases.

CheckM2 evaluations inside the merge loops are additionally cached by sequence content in `Checkm2_metrics_cache.tsv`. Unchanged candidates — including the original target bin recopied into every merge batch and candidates recurring across convergence iterations — are reused without re-annotation, and only new content is evaluated in one batched multi-threaded `checkm2 predict` call (progress lines read `CheckM2 cache: X/Y candidate(s) reused, Z evaluated`). The cache resets on `--mode new` and invalidates itself when the CheckM2 version or `CHECKM2DB` database changes; after swapping the database manually, remove the file or start a new run. If the batched prediction fails, BASALT falls back to the historical direct per-folder prediction.

### The run is too slow or uses too much storage

- use `--sensitive quick --refinepara quick`;
- remove optional binners;
- run a pilot on one assembly;
- request adequate local scratch storage;
- monitor memory and I/O rather than assuming CPU is the bottleneck.

Changing the preset changes the candidate space. Report the change as a method setting, not only as a performance optimization.

## Outputs and quality

### Where are final MAGs?

For the default CheckM2 route, they are in the directory named by `-o`. Without `-o`, use `Final_binset/`.

### Why is `quality_report.tsv` missing?

Possible causes include:

- CheckM2 database or executable failure;
- no bins reaching a stage, including after low-coverage candidate generation;
- an external stage failing before quality assessment;
- output being written under a stage-specific directory;
- cleanup or manual movement.

Search logs before lowering scientific quality thresholds. Lowering `--min-cpn` changes which bins enter refinement and is not a general fix for a missing database or failed program.

Start by distinguishing a missing report from a low-coverage result:

```bash
find . -type f -name 'quality_report.tsv' -print
find . -maxdepth 2 -type f \( -name '*.fa' -o -name '*.fasta' \) -size +0c | wc -l
grep -Ei 'coverage|no bins|quality_report|diamond|error|failed' \
  Basalt_log.txt basalt.stderr.log 2>/dev/null
```

If the mapping and binning stages produced few or no candidates, inspect mapped-read fraction, depth distribution, assembly fragmentation, and read-to-assembly provenance. Additional sequencing or a biologically justified co-assembly may improve coverage; lowering a quality threshold cannot create missing read support. If candidates exist but the report does not, test CheckM2 independently on one non-empty bin and verify its database before rerunning BASALT.

### Does a folder named `refined` or `MAGs` guarantee acceptable genomes?

No. Directory names encode workflow state. Apply explicit completeness, contamination, taxonomic, size, and study-specific criteria to the final report.

### Can the same contig appear in several candidate bins?

Yes, especially before within- and cross-assembly selection. BASALT reduces candidate redundancy, but the final biological interpretation still requires a stated dereplication and taxonomic framework.

## Data types and acceleration

### Which read designs are supported?

BASALT accepts assembly FASTA files with paired-end short reads, ONT or PacBio CLR reads, PacBio HiFi reads, or compatible combinations. Reassembly and polishing stages depend on which read types are present and can be skipped when their required inputs are absent.

### Does BASALT use a GPU?

Some dependencies, including SemiBin 2, PyTorch, and optional VAMB configurations, may use compatible accelerators. Availability depends on the installed packages and drivers. Do not assume GPU use from the BASALT command alone. Record device utilization and software versions if acceleration affects runtime or reproducibility.

## Optional binners

### BASALT finished after an optional binner failed

This is expected failure handling in the current CheckM2 orchestration. Optional-adapter exceptions are logged and the pipeline continues. Inspect the log and verify which candidates actually entered selection.

See [Extra binners](extra-binners.md).

## BASALT versus BASALT-Air

BASALT-Air is a separate Pixi-based implementation with different path and output controls. The Conda BASALT commands on this site do not support `--workdir`, `--outdir`, `--version`, or `--check-deps`.

Do not assume two implementations produce byte-identical outputs without a version-matched regression test. Record which repository, release, dependency lock, models, and databases were used.

## Reporting a bug

Open a [GitHub issue](https://github.com/PKU-EMBL/BASALT/issues) with:

- BASALT commit or release;
- operating system and resource allocation;
- Conda environment export;
- complete command with sensitive paths redacted consistently;
- `Basalt_log.txt`, `Basalt_checkpoint.txt`, stdout, and stderr;
- model and database configuration;
- relevant external-tool versions;
- a minimal shareable reproducer when possible.

Do not upload human, clinical, or otherwise restricted sequencing data to a public issue.
