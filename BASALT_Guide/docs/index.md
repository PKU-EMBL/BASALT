# BASALT documentation

::::{container} basalt-hero
<p class="basalt-hero__eyebrow">Genome-resolved metagenomics · PKU-EMBL Laboratory</p>

<p class="basalt-hero__tagline">Recover and refine metagenome-assembled genomes across one or more assemblies in a checkpointed workflow.</p>

<p class="basalt-hero__summary">BASALT connects multi-binner candidate generation, within- and between-assembly bin selection, contig-level contamination screening, read-supported retrieval, and reassembly.</p>

<p class="basalt-hero__actions"><a class="basalt-button basalt-button--primary" href="quickstart.html">Get started</a><a class="basalt-button basalt-button--secondary" href="basalt-air.html">Discover BASALT-Air</a><a class="basalt-button basalt-button--quiet" href="study-design.html">Plan a study</a><a class="basalt-button basalt-button--quiet" href="pipeline.html">Explore the method</a></p>

<p class="basalt-hero__badges"><a href="release-notes.html">BASALT v1.2.2</a><a href="https://doi.org/10.1038/s41467-024-46539-7">Nature Communications 2024</a><a href="https://github.com/PKU-EMBL/BASALT-Air">New: BASALT-Air v1.0.0</a><span>Multi-assembly</span><span>MIT licensed</span></p>
::::

:::{admonition} Current release: BASALT v1.2.2
:class: tip

Version 1.2.2 improves environment reproducibility, model downloads, mainland-China installation, and documentation. Read the concise [release notes](release-notes.md#basalt-122) before updating an existing environment or resuming an older analysis.
:::

:::{admonition} Recent development improvements (ahead of the next tagged release)
:class: note

The external-binset route now supports full gap filling: after data feeding (`-d`) and dereplication (`-b`), run `BASALT -r BestBinset --module all` from the data-feeding output directory to continue through contig retrieval, secondary dereplication, restrained OLC, and reassembly. OLC verification is parallelized with a per-worker thread budget, CheckM2 evaluations are reused across convergence iterations via a content-keyed cache, and interrupted runs resume cleanly. See the [external-binset workflow](usage.md#external-binsets-staged-workflow) and the [unreleased notes](release-notes.md#unreleased).
:::

:::{admonition} New implementation: BASALT-Air v1.0.0
:class: note

[BASALT-Air](https://github.com/PKU-EMBL/BASALT-Air) is now available for new deployments. It uses Pixi, accepts absolute input paths, separates working and output directories, and records a structured run manifest. Read [Choosing BASALT-Air](basalt-air.md) before switching: the executable, sample delimiters, path options, and checkpoints are not interchangeable with this Conda edition.
:::

:::{admonition} Evidence boundary
:class: important

The published study reports improved MAG recovery and downstream genome-resolved analyses relative to the workflows and datasets tested. This is not a guarantee for every community or sequencing design. Dataset complexity, sequencing depth, assembly quality, read type, database choice, and parameter settings can change the outcome.
:::

## Core advantage: coherent multi-assembly refinement

BASALT can jointly evaluate sample-specific assemblies and a biologically matched co-assembly. The former can preserve sample-restricted populations while avoiding some pooling complexity; the latter can increase effective depth for genomes shared across samples. In related longitudinal, spatial, technical, intervention, stage, or environmental-gradient designs, concordant cross-sample coverage changes provide an abundance fingerprint that complements tetranucleotide composition and sequence overlap during candidate selection, dereplication, contig screening, retrieval, and reassembly.

This strategy is most informative when samples share substantial genomic content while retaining biologically meaningful abundance variation. It does not justify pooling unrelated samples or guarantee better MAGs. Use a prespecified pilot to compare individual-only, pooled-only, and combined designs under identical quality criteria. See [Study design patterns](study-design.md) for supported layouts and decision boundaries.

## Start here

Choose the route that matches the decision you need to make. Commands and defaults on this site are checked against the current repository.

Examples beginning with `BASALT` assume that the `basalt` environment is active. For activation-free Micromamba execution, replace `BASALT` with `micromamba run -n basalt BASALT`.

::::{grid} 1 2 3 3
:gutter: 3

:::{grid-item-card} Install and verify
:link: installation
:link-type: doc

Configure BASALT, model weights, quality-control databases, external executables, and containers.
:::

:::{grid-item-card} Run a minimal analysis
:link: quickstart
:link-type: doc

Launch an isolated short-read run with explicit resources, mode, and output name.
:::

:::{grid-item-card} Design assemblies and coverage
:link: study-design
:link-type: doc

Choose between sample-specific, pooled, longitudinal, and multi-assembly strategies.
:::

:::{grid-item-card} Understand the method
:link: pipeline
:link-type: doc

Follow each executable stage and distinguish software behaviour from interpretation.
:::

:::{grid-item-card} Audit the outputs
:link: output
:link-type: doc

Locate final MAGs, inspect quality reports, and diagnose incomplete stages.
:::

:::{grid-item-card} Report the analysis
:link: reproducibility
:link-type: doc

Archive provenance and use the minimum Methods and reporting checklist.
:::
::::

## BASALT framework

```{figure} img/workflow.png
:alt: BASALT framework from multi-assembly candidate generation through bin selection, refinement, retrieval, polishing, and reassembly
:class: basalt-framework
:align: center

The BASALT framework integrates multi-binner candidate generation, within- and cross-assembly selection, model-based contig screening, read-supported retrieval, polishing, and reassembly. See [Pipeline and methods](pipeline.md) for stage-level behaviour and interpretation boundaries.
```

## Method at a glance

| Stage | Function | Main evidence used | Principal boundary |
|---|---|---|---|
| Autobinning | Generate candidate bins under several algorithms or parameter settings | Sequence composition and read depth | Candidate diversity depends on the selected preset and optional binners |
| Bin selection | Compare bins within and across assemblies | Bin quality estimates, overlap, composition, and coverage | Dereplication does not establish biological strain identity |
| Refinement | Screen candidate contaminant contigs and retrieve supported contigs | Composition, coverage-derived features, and read connectivity | Predictions depend on the trained models and available coverage signals |
| Reassembly | Reassemble selected bins when compatible reads are present | Reads assigned to individual bins | Reassembly can be skipped or fail for bins with insufficient read support |

The [pipeline page](pipeline.md) separates what the software executes from how results should be interpreted.

## Supported analysis designs

Every run requires one or more assembly FASTA files and at least one read source that can provide coverage information.

| Design | Short reads | ONT or PacBio CLR | PacBio HiFi |
|---|---:|---:|---:|
| Short-read only | required | no | no |
| Long-read only | no | required | no |
| HiFi only | no | no | required |
| Short reads + ONT or CLR | required | required | no |
| Short reads + HiFi | required | no | required |

Read [Input formats](inputs.md) and [Study design patterns](study-design.md) before combining assemblies or samples. BASALT does not infer which reads belong to which biological sample.

## BASALT and BASALT-Air

This site primarily documents the Conda-based BASALT repository. [BASALT-Air v1.0.0](https://github.com/PKU-EMBL/BASALT-Air) is the newer Pixi-based implementation with absolute-path handling, separate work and output directories, dependency checks, and structured run manifests. See [Choosing BASALT-Air](basalt-air.md) for a concise comparison.

Use the documentation that belongs to the implementation and version you ran. Do not mix BASALT-Air separators or path options into a BASALT command.

## Evidence and citation

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

Also cite optional binners and quality-control software when they materially affect the reported result.

## Scope of this documentation

The guide distinguishes three kinds of statement:

- **Executable behaviour**, checked against the current repository.
- **Scientific interpretation**, limited to what the workflow and published evidence support.
- **Operational guidance**, identified as a recommendation when it is not a program requirement.

If documentation and observed behaviour differ, retain the full command, `Basalt_log.txt`, `Basalt_checkpoint.txt`, software versions, and a minimal reproducer. Report the discrepancy through [GitHub Issues](https://github.com/PKU-EMBL/BASALT/issues).

```{toctree}
:hidden:
:maxdepth: 2
:caption: Getting started

installation
basalt-air
quickstart
tutorial
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: User guide

inputs
study-design
usage
pipeline
output
extra-binners
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: Reproducibility

reproducibility
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: Help and development

faq
release-notes
developer
```
