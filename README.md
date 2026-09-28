# BASALT: Binning Across a Series of Assemblies Toolkit

[![License: MIT](https://img.shields.io/badge/license-MIT-2f855a.svg)](LICENSE)
[![Release: v1.2.2](https://img.shields.io/badge/release-v1.2.2-2563eb.svg)](https://github.com/PKU-EMBL/BASALT/releases/tag/V1.2.2)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776ab.svg)](https://www.python.org/)
[![Nature Communications](https://img.shields.io/badge/Nature%20Communications-2024-b31b1b.svg)](https://doi.org/10.1038/s41467-024-46539-7)
[![Documentation](https://readthedocs.org/projects/basalt-guide/badge/?version=latest)](https://basalt-guide.readthedocs.io/en/latest/)
[![BASALT-Air](https://img.shields.io/badge/BASALT--Air-v1.0.0-55ad62.svg)](https://github.com/PKU-EMBL/BASALT-Air)

BASALT is a metagenomic workflow for recovering and refining metagenome-assembled genomes (MAGs) across one or more assemblies. It combines multi-binner candidate generation, within- and between-assembly bin selection, contig-level contamination screening, contig retrieval, and reassembly in one checkpointed pipeline.

> [!IMPORTANT]
> **New: [BASALT-Air v1.0.0](https://github.com/PKU-EMBL/BASALT-Air) is now available.** BASALT-Air is the newer Pixi-based implementation for users who need absolute input paths, separate working and output directories, reproducible run manifests, and built-in dependency checks. Its executable is lowercase `basalt`, and its CLI and checkpoints are not interchangeable with the Conda-based BASALT documented in this repository.

> [!NOTE]
> **Recent development improvements (ahead of the next tagged release).** The external-binset route (`-d`/`-b`/`-r`) now reaches full gap filling — run `BASALT -r BestBinset --module all` from the data-feeding output directory to continue through contig retrieval, secondary dereplication, restrained OLC, and reassembly on both CheckM2 and legacy CheckM. The OLC verification stages are parallelized with a `-t`-derived per-worker thread budget, CheckM2 evaluations are reused across convergence iterations through a content-keyed cache, and interrupted runs resume cleanly: stale per-bin intermediates are discarded automatically and non-numeric checkpoint lines no longer reset progress. See the [external-binset workflow](BASALT_Guide/docs/usage.md#external-binsets-staged-workflow), [resume semantics](BASALT_Guide/docs/usage.md#resume-semantics), and [release notes](BASALT_Guide/docs/release-notes.md#unreleased).

The published evaluation reports improved MAG recovery and downstream genome-resolved analyses relative to the tested workflows and datasets. These results define the evidence boundary: performance depends on community complexity, sequencing depth, assembly quality, read type, and the selected parameters. See the [Nature Communications article](https://doi.org/10.1038/s41467-024-46539-7) for the complete experimental design.

## Core advantage: coherent multi-assembly refinement

BASALT can jointly evaluate sample-specific assemblies and a biologically matched co-assembly. Individual assemblies preserve sample-restricted populations while avoiding some pooling complexity; co-assembly can increase effective depth for shared, low-abundance genomes. Across related samples, concordant coverage changes provide an abundance fingerprint that complements tetranucleotide composition and sequence overlap during bin selection, dereplication, contig screening, retrieval, and reassembly. This design is especially useful for longitudinal, reactor time-series, nearby-site, extraction-method, intervention, growth-stage, and environmental-gradient studies. Validate it in a pilot; it is not a reason to combine unrelated samples. See [Study design patterns](BASALT_Guide/docs/study-design.md).

<p align="center">
  <img src="fig/workflow.png" alt="BASALT workflow from assemblies and reads to refined MAGs" width="86%">
</p>

## Documentation

The [BASALT guide](https://basalt-guide.readthedocs.io/en/latest/) is the source of truth for installation, inputs, commands, outputs, troubleshooting, and reproducible reporting. Its editable sources are under [`BASALT_Guide/docs/`](BASALT_Guide/docs/index.md).

- [Installation](BASALT_Guide/docs/installation.md)
- [BASALT-Air: the new implementation](BASALT_Guide/docs/basalt-air.md)
- [Quick start](BASALT_Guide/docs/quickstart.md)
- [Input formats](BASALT_Guide/docs/inputs.md)
- [Study design patterns](BASALT_Guide/docs/study-design.md)
- [Command-line reference](BASALT_Guide/docs/usage.md)
- [Tutorial](BASALT_Guide/docs/tutorial.md)
- [Reproducibility and reporting](BASALT_Guide/docs/reproducibility.md)

## Choose a route

| Goal | Recommended route | Read first |
|---|---|---|
| Start a new deployment with modern path and environment handling | BASALT-Air v1.0.0 | [BASALT-Air overview](BASALT_Guide/docs/basalt-air.md) |
| Reproduce or resume an existing Conda BASALT analysis | This BASALT repository | [Installation](BASALT_Guide/docs/installation.md) |
| First successful run | One assembly and one matched paired-end sample | [Quick start](BASALT_Guide/docs/quickstart.md) |
| Compare individual and pooled assemblies | Run a biologically coherent multi-assembly pilot | [Study design patterns](BASALT_Guide/docs/study-design.md) |
| Recover candidates with lower compute demand | `--sensitive quick --refinepara quick` | [Command-line reference](BASALT_Guide/docs/usage.md) |
| Expand the candidate pool | `--sensitive sensitive` or `more-sensitive`; optionally add a validated binner | [Extra binners](BASALT_Guide/docs/extra-binners.md) |
| Import bins from another workflow | Data feeding (`-d`), dereplication (`-b`), then `-r BestBinset --module all` for refinement and gap filling | [External-binset workflow](BASALT_Guide/docs/usage.md#external-binsets-staged-workflow) |
| Screen an existing binset only | `BASALT -r <binset> ...` without `--module` (checkpointed; continue later with `--module`) | [Command-line reference](BASALT_Guide/docs/usage.md) |
| Resume an interrupted run | `BASALT --mode continue` in the unchanged run directory | [Resume semantics](BASALT_Guide/docs/usage.md#resume-semantics) |

### BASALT or BASALT-Air?

| Capability | BASALT in this repository | BASALT-Air v1.0.0 |
|---|---|---|
| Environment | Conda plus external installation steps | Pixi environment and lock file |
| Command | `BASALT` | `basalt` |
| Input paths | Dedicated run directory and local names recommended | First-class absolute paths |
| Work and output locations | Intermediates in the current directory | Separate `--workdir` and `--outdir` |
| Run provenance | Command, logs, and checkpoint files | Adds a structured run manifest and timestamped logs |
| Best suited to | Existing BASALT analyses and checkpoint reproduction | New deployments needing clearer path and environment isolation |

See the [BASALT-Air repository](https://github.com/PKU-EMBL/BASALT-Air) for its authoritative installation and command reference.

## What BASALT does

1. **Autobinning.** BASALT runs MetaBAT 2 and SemiBin 2 across parameter settings. The `sensitive` and `more-sensitive` presets add CONCOCT, and `more-sensitive` also adds MaxBin 2.0.
2. **Bin selection.** Candidate bins are compared within each assembly and then across assemblies to reduce redundancy.
3. **Refinement.** A five-model multilayer-perceptron ensemble screens contigs using sequence-composition and coverage-derived features. Read connectivity is then used to retrieve candidate contigs.
4. **Reassembly and gap filling.** Selected and redundant bins are reassembled and merged through restrained OLC when the required read types are available; this stage now also runs on external binsets through `-r ... --module all`, with parallel verification and content-cached CheckM2 evaluation.

Optional adapters support MetaBinner (`-e m`), VAMB (`-e v`), and LorBin (`-e l`). Their installation and input requirements are documented separately in the [extra-binner guide](BASALT_Guide/docs/extra-binners.md).

The LorBin adapter is maintained against [PKU-EMBL/LorBin-BASALT-Extrabinner](https://github.com/PKU-EMBL/LorBin-BASALT-Extrabinner), an in-house fork that adapts LorBin to the BASALT Python 3.12 environment and carries BASALT-specific runtime fixes. Record the fork's Git commit for every analysis that uses `-e l`; the package version alone does not uniquely identify these changes.

### Candidate-generation presets

| Preset | Candidate generators | Practical trade-off |
|---|---|---|
| `quick` | MetaBAT 2 and SemiBin 2 | Lower compute demand and fewer candidate sets |
| `sensitive` | MetaBAT 2, CONCOCT, and SemiBin 2 | Default balance of candidate diversity and runtime |
| `more-sensitive` | MaxBin 2.0, MetaBAT 2, CONCOCT, and SemiBin 2 | Higher compute and storage demand; not guaranteed to recover more acceptable MAGs |

## Supported input designs

Every run requires at least one assembly and at least one source of read coverage. BASALT accepts the following designs:

| Design | Assembly | Short reads | ONT or PacBio CLR | PacBio HiFi |
|---|---:|---:|---:|---:|
| Short-read only | required | required | no | no |
| Long-read only | required | no | required | no |
| HiFi only | required | no | no | required |
| Short reads + ONT or CLR | required | required | required | no |
| Short reads + HiFi | required | required | no | required |

Input paths are resolved from the working directory. Use filenames or relative paths without shell metacharacters; for the most reliable execution, place or symlink inputs into a dedicated run directory. Paired-end files within a sample are separated by a comma, and multiple short-read samples are separated by `/`.

For replicate or longitudinal studies, decide the biological strata before assembling or pooling. BASALT can compare sample-specific and matched pooled assemblies, but it cannot infer whether a pool erased a time, treatment, host, or habitat contrast. The [study-design guide](BASALT_Guide/docs/study-design.md) provides concrete layouts and a minimum design manifest.

## Requirements

- Linux x86-64
- Python 3.12
- 8 CPU threads and 128 GB RAM as a practical starting point
- 32 or more threads and 256 GB RAM for larger or more complex datasets
- Sufficient temporary storage for mapping, binning, and reassembly intermediates

Resource use is dataset-dependent. Run a pilot or request more memory for large assemblies, many samples, or `more-sensitive` mode.

## Installation

The Conda-compatible installation below is the maintained route for this repository. [Micromamba](https://mamba.readthedocs.io/en/stable/installation/micromamba-installation.html) is recommended for faster dependency solving and parallel package downloads; current Conda with the libmamba solver is also supported. See the [installation guide](BASALT_Guide/docs/installation.md) for databases, model weights, mainland-China mirrors, verification, and container execution notes.

```bash
git clone --branch V1.2.2 --depth 1 \
  https://github.com/PKU-EMBL/BASALT.git BASALT
cd BASALT

micromamba create -n basalt -f basalt_environment.yml \
  --strict-channel-priority --yes
micromamba run -n basalt bash install.sh
```

The release tag makes this installation reproducible instead of following a moving `master` branch. Clone `master` without `--branch V1.2.2` only for development.

Run the installer with `bash install.sh`; the source file does not need to be executable. The installer assigns mode `0755` only to installed Python and Perl scripts and the generated `BASALT` launcher. The platform-compatible `jgi_summarize_bam_contig_depths` executable is supplied by the maintained MetaBAT 2 package and is validated during installation; BASALT no longer ships a stale Linux binary. Do not use `sudo`, `chmod 777`, or `chmod -R 777` on the repository or environment. If installation reports a permission error, use a user-owned Conda/Micromamba environment and follow the [permission troubleshooting guide](BASALT_Guide/docs/faq.md#permission-denied-or-should-i-run-chmod-777).

Activate the environment before using the unprefixed commands shown later in this README:

```bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate basalt
# Conda alternative: conda activate basalt
```

If shell activation is undesirable, prefix each command with `micromamba run -n basalt`, for example `micromamba run -n basalt BASALT --help`. Conda users can create the environment with `CONDA_CHANNEL_PRIORITY=strict conda env create -n basalt -f basalt_environment.yml --solver libmamba --yes` and replace `micromamba run` with `conda run`. [`basalt_environment.yml`](basalt_environment.yml) is the sole environment definition. It uses the portable `conda-forge` and `bioconda` channel names and resolves the complete base stack in one Conda transaction; it does not hard-code a site mirror.

As an alternative to the manual environment commands above, mainland-China users can run the network-aware installer. It probes TUNA, BFSU, USTC, and upstream endpoints, prefers `micromamba`, applies matching Conda and PyPI mirrors only to the installation subprocess, and leaves `~/.condarc` and pip configuration unchanged:

```bash
python3 BASALT_setup_China_mainland.py \
  --manager micromamba \
  --bootstrap-micromamba \
  --mirror auto \
  --model-dir /absolute/persistent/path/BASALT_WEIGHT
```

Use `--dry-run` to inspect all commands first, `--mirror tuna|bfsu|ustc|upstream` to pin a route, or `--update` to synchronize an existing environment. When this installer completes the model step successfully, do not download the same weights again; continue with validation and the shell configuration below.

Download the trained BASALT models into a persistent directory and verify them before configuring the shell. The downloader validates all five ensemble descriptors and every checkpoint they reference:

```bash
micromamba run -n basalt BASALT_models_download.py \
  --source auto \
  --path /absolute/persistent/path/BASALT_WEIGHT

find /absolute/persistent/path/BASALT_WEIGHT \
  -maxdepth 1 -name '*_ensemble.csv' | wc -l
```

The expected count is `5`. Inspect `~/.bashrc` before adding the validated path for future login and batch shells:

```bash
grep -n '^export BASALT_WEIGHT=' ~/.bashrc || true
# If no current definition was printed, append one:
printf '%s\n' \
  'export BASALT_WEIGHT="/absolute/persistent/path/BASALT_WEIGHT"' \
  >> ~/.bashrc
source ~/.bashrc
test -d "$BASALT_WEIGHT"
```

If `grep` found an existing definition, edit that line instead of appending another one. Keep exactly one authoritative `BASALT_WEIGHT` export.

Automatic mode performs a bounded reachability check, downloads concurrently from the official [PKU-EMBL/BASALT_WEIGHT](https://huggingface.co/PKU-EMBL/BASALT_WEIGHT) repository first, and falls back to the legacy Figshare archive when Hugging Face is unreachable. Adjust the check with `--hf-timeout SECONDS`. A local ZIP, a custom URL, and the documented Baidu Netdisk route are available when either service is unsuitable.

The downloader pins the release model revision by default. If a site-approved Hugging Face-compatible endpoint is required, pass `--hf-endpoint URL`; record the selected revision and hash the downloaded files for the analysis record.

Set the CheckM2 database path if it is not already configured by your environment:

```bash
checkm2 database --download
checkm2 database --current
```

Only when automatic discovery is unavailable:

```bash
export CHECKM2DB=/absolute/path/to/CheckM2_database/uniref100.KO.1.dmnd
```

If the export is required across logins, keep one validated absolute definition in `~/.bashrc`, following the same duplicate-check and scheduler precautions used for `BASALT_WEIGHT`.

Verify the command and external tools before a production run:

```bash
BASALT --help
command -v \
  checkm2 metabat2 jgi_summarize_bam_contig_depths SemiBin2 \
  bowtie2 samtools spades.py
```

## Quick start

Create a clean working directory and link the inputs into it:

```bash
mkdir -p basalt_run
cd basalt_run
ln -s /absolute/path/assembly.fasta .
ln -s /absolute/path/sample_R1.fastq .
ln -s /absolute/path/sample_R2.fastq .
```

Run the full short-read workflow:

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -t 32 \
  -m 128 \
  --mode new \
  -o study_basalt
```

The final MAGs and final CheckM2 report are written under `study_basalt/`. Intermediate files, `BASALT_command.txt`, `Basalt_log.txt`, and `Basalt_checkpoint.txt` remain in the working directory.

To resume an interrupted run from that same directory:

```bash
BASALT --mode continue
```

Do not combine files from different runs in one working directory. `--mode new` initializes checkpoint state, but it does not provide complete cleanup or provenance isolation.

## Common commands

Short reads across two assemblies and two samples:

```bash
BASALT \
  -a assembly_1.fasta,assembly_2.fasta \
  -s sample_1_R1.fastq,sample_1_R2.fastq/sample_2_R1.fastq,sample_2_R2.fastq \
  -t 64 -m 256 --mode new -o study_basalt
```

Short reads plus ONT or PacBio CLR reads:

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -l nanopore.fastq \
  -t 64 -m 256 --mode new -o study_basalt
```

PacBio HiFi reads:

```bash
BASALT \
  -a assembly.fasta \
  -hf hifi.fastq \
  -t 32 -m 128 --mode new -o study_basalt
```

Faster candidate generation and refinement:

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  --sensitive quick --refinepara quick \
  -t 32 -m 128 --mode new -o study_basalt
```

See the [command-line reference](BASALT_Guide/docs/usage.md) before using standalone refinement, data feeding, dereplication, or module-specific execution. These routes require intermediate files from a compatible BASALT run.

## Output contract

For the default CheckM2 route, `-o NAME` produces the final directory `NAME/`; without `-o`, the default is `Final_binset/`. A scientifically complete result should contain non-empty bin FASTA files and a matching quality report, while the working directory should retain `BASALT_command.txt`, `Basalt_log.txt`, and `Basalt_checkpoint.txt`.

Directory names such as `BestBinset`, `refined`, or `MAGs` describe workflow state, not publication-ready quality. Apply explicit completeness, contamination, taxonomy, contiguity, read-support, and dereplication criteria. See [Outputs and quality control](BASALT_Guide/docs/output.md).

## Reproducible use

For each analysis, archive:

- the exact BASALT command and BASALT Git commit or release;
- the Conda environment export and versions of external executables;
- CheckM or CheckM2 database identity and quality thresholds;
- input assembly and read checksums;
- `BASALT_command.txt`, `Basalt_log.txt`, and `Basalt_checkpoint.txt`;
- the final quality report and all exclusions applied after BASALT.

The [reporting guide](BASALT_Guide/docs/reproducibility.md) provides a Methods template and a minimum reporting checklist.

## Citation

If BASALT contributes to an analysis, cite:

> Qiu, Z., Yuan, L., Lian, C.-A. *et al.* BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis. *Nature Communications* **15**, 2179 (2024). [https://doi.org/10.1038/s41467-024-46539-7](https://doi.org/10.1038/s41467-024-46539-7)

```bibtex
@article{qiu2024basalt,
  author  = {Qiu, Zhiguang and Yuan, Li and Lian, Chun-Ang and Lin, Bin and
             Chen, Jie and Mu, Rong and Qiao, Xuejiao and Zhang, Liyu and
             Xu, Zheng and Fan, Lu and Zhang, Yunzeng and Wang, Shanquan and
             Li, Junyi and Cao, Huiluo and Li, Bing and Chen, Baowei and
             Song, Chi and Liu, Yongxin and Shi, Lili and Tian, Yonghong and
             Ni, Jinren and Zhang, Tong and Zhou, Jizhong and Zhuang, Wei-Qin and
             Yu, Ke},
  title   = {BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis},
  journal = {Nature Communications},
  volume  = {15},
  pages   = {2179},
  year    = {2024},
  doi     = {10.1038/s41467-024-46539-7}
}
```

If an optional binner materially contributes to the reported analysis, cite that software and record its version and parameters.

## Support and licence

Report reproducible bugs through [GitHub Issues](https://github.com/PKU-EMBL/BASALT/issues). Include the command, log, checkpoint, software versions, database configuration, and the smallest shareable input that reproduces the problem.

Scientific or collaboration enquiries:

- Ke Yu: [yuke.sz@pku.edu.cn](mailto:yuke.sz@pku.edu.cn)
- Zhaorui (Elijah) Jiang: [zrjiang25@stu.pku.edu.cn](mailto:zrjiang25@stu.pku.edu.cn)

BASALT is released under the [MIT License](LICENSE).
