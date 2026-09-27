#!/usr/bin/env python3

"""Bounded release-smoke tests for BASALT's dependency-light interfaces."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from zipfile import ZipFile
from pathlib import Path
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
MODULES = REPOSITORY / "BASALT"
sys.path.insert(0, str(MODULES))

from basalt_runtime import (  # noqa: E402
    append_checkpoint,
    clean_stale_olc_intermediates,
    coerce_count,
    depth_sample_count,
    evaluate_bins_checkm2,
    expand_alignment_group,
    fasta_content_hash,
    fasta_has_sequence,
    pe_read_names,
    read_checkpoint_step,
    require_model_directory,
    worker_thread_budget,
)
from external_workflow import (  # noqa: E402
    ExternalBinsetError,
    external_module_plan,
    pair_connections,
    prepare_external_binset_route,
)
from BASALT_models_download import unpack_model, validate_models  # noqa: E402
from S1e_extra_binners import (  # noqa: E402
    _metabinner_executable,
    _read_vamb_assignments,
    _write_assigned_bins,
    lorbin,
    vamb,
)


class ModelDirectoryTests(unittest.TestCase):
    def test_accepts_five_descriptor_and_checkpoint_pairs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model_dir = Path(temporary)
            for index in range(5):
                descriptor = model_dir / f"model_{index}_ensemble.csv"
                descriptor.write_text("checkpoint.pth\n", encoding="utf-8")
                checkpoint_dir = descriptor.with_suffix("")
                checkpoint_dir.mkdir()
                (checkpoint_dir / "checkpoint.pth").write_bytes(b"test checkpoint")

            with patch.dict(os.environ, {"BASALT_WEIGHT": str(model_dir)}):
                self.assertEqual(require_model_directory(), model_dir.resolve())
                self.assertEqual(len(validate_models(model_dir)), 5)

    def test_rejects_an_unset_model_directory(self) -> None:
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(RuntimeError, "BASALT_WEIGHT is not set"):
                require_model_directory()

    def test_rejects_a_descriptor_with_a_missing_checkpoint(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model_dir = Path(temporary)
            for index in range(5):
                descriptor = model_dir / f"model_{index}_ensemble.csv"
                descriptor.write_text("missing.pth\n", encoding="utf-8")
                descriptor.with_suffix("").mkdir()

            with self.assertRaisesRegex(RuntimeError, "Missing model checkpoint files"):
                validate_models(model_dir)

    def test_rejects_archive_path_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / "unsafe.zip"
            with ZipFile(archive, "w") as zipped:
                zipped.writestr("../outside.txt", "unsafe")

            with self.assertRaisesRegex(ValueError, "Unsafe path"):
                unpack_model(archive, root / "models")
            self.assertFalse((root / "outside.txt").exists())

    def test_rejects_checkpoint_path_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model_dir = Path(temporary)
            for index in range(5):
                descriptor = model_dir / f"model_{index}_ensemble.csv"
                descriptor.write_text("../outside.pth\n", encoding="utf-8")
                descriptor.with_suffix("").mkdir()
            (model_dir / "outside.pth").write_bytes(b"not a model")

            with self.assertRaisesRegex(RuntimeError, "Unsafe checkpoint paths"):
                validate_models(model_dir)


class OptionalBinnerTests(unittest.TestCase):
    def test_metabinner_home_requires_an_executable_script(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            home = Path(temporary)
            script = home / "run_metabinner.sh"
            script.write_text("#!/bin/sh\n", encoding="utf-8")
            with patch.dict(
                os.environ,
                {"METABINNER_HOME": str(home), "PATH": ""},
                clear=False,
            ):
                with self.assertRaisesRegex(RuntimeError, "executable"):
                    _metabinner_executable()
                script.chmod(0o755)
                self.assertEqual(_metabinner_executable(), script.resolve())

    def test_reads_vamb_5_unsplit_cluster_table(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            (output_dir / "vae_clusters_unsplit.tsv").write_text(
                "clustername\tcontigname\n1\tcontig_a\n2\tcontig_b\n",
                encoding="utf-8",
            )
            self.assertEqual(
                _read_vamb_assignments(output_dir),
                {"contig_a": "1", "contig_b": "2"},
            )

    def test_materializes_cluster_assignments_as_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assembly = root / "assembly.fa"
            assembly.write_text(
                ">contig_a\nAAAA\n>contig_b\nCCCCCC\n", encoding="utf-8"
            )
            output_dir = root / "bins"
            paths = _write_assigned_bins(
                assembly,
                {"contig_a": "cluster 1", "contig_b": "cluster 2"},
                output_dir,
                "candidate",
                minimum_size=4,
            )
            self.assertEqual([path.name for path in paths], [
                "candidate.cluster_1.fa",
                "candidate.cluster_2.fa",
            ])
            self.assertTrue(all(path.stat().st_size > 0 for path in paths))

    def test_vamb_adapter_uses_v5_command_and_sorted_bam(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "1_assembly.fa").write_text(
                ">contig_a\n" + "A" * 500_000 + "\n", encoding="utf-8"
            )
            (root / "1_DNA-1_sorted.bam").write_bytes(b"bam")
            commands = []

            def fake_run(command, cwd=None):
                commands.append([str(item) for item in command])
                output_dir = Path(command[command.index("--outdir") + 1])
                output_dir.mkdir()
                (output_dir / "vae_clusters_unsplit.tsv").write_text(
                    "clustername\tcontigname\n1\tcontig_a\n", encoding="utf-8"
                )

            with patch("S1e_extra_binners._run", side_effect=fake_run), patch(
                "S1e_extra_binners._run_quality_check"
            ):
                vamb(
                    "1_assembly.fa",
                    {"1": ["reads_R1.fastq", "reads_R2.fastq"]},
                    8,
                    str(root),
                    "checkm2",
                )

            self.assertEqual(commands[0][:3], ["vamb", "bin", "default"])
            self.assertIn("--bamfiles", commands[0])
            self.assertIn("--minfasta", commands[0])
            self.assertIn("-p", commands[0])
            self.assertIn(str((root / "1_DNA-1_sorted.bam").resolve()), commands[0])
            self.assertTrue((root / "1_assembly.fa_100_vamb_genomes").is_dir())

    def test_lorbin_adapter_uses_sorted_bams_and_multi_mode(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "1_assembly.fa").write_text(">contig_a\nAAAA\n", encoding="utf-8")
            for index in (1, 2):
                (root / f"1_DNA-{index}_sorted.bam").write_bytes(b"bam")
            commands = []

            def fake_run(command, cwd=None):
                commands.append([str(item) for item in command])
                output_dir = Path(command[command.index("-o") + 1])
                output_dir.mkdir()
                (output_dir / "lorbin.1.fa").write_text(
                    ">contig_a\nAAAA\n", encoding="utf-8"
                )

            with patch("S1e_extra_binners._run", side_effect=fake_run), patch(
                "S1e_extra_binners._run_quality_check"
            ):
                lorbin(
                    "1_assembly.fa",
                    {
                        "1": ["sample_1_R1.fastq", "sample_1_R2.fastq"],
                        "2": ["sample_2_R1.fastq", "sample_2_R2.fastq"],
                    },
                    8,
                    str(root),
                    "checkm2",
                )

            self.assertEqual(commands[0][:2], ["LorBin", "bin"])
            self.assertIn("--multi", commands[0])
            self.assertIn("--num_process", commands[0])
            self.assertIn(str((root / "1_DNA-1_sorted.bam").resolve()), commands[0])
            self.assertIn(str((root / "1_DNA-2_sorted.bam").resolve()), commands[0])
            self.assertTrue((root / "1_assembly.fa_100_lorbin_genomes").is_dir())


class CommandLineTests(unittest.TestCase):
    def run_cli(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(MODULES / "BASALT.py"), *arguments],
            cwd=REPOSITORY,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_help_is_available_without_pipeline_imports(self) -> None:
        result = self.run_cli("--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--assemblies", result.stdout)

    def test_rejects_nonpositive_threads(self) -> None:
        result = self.run_cli("--threads", "0")
        self.assertEqual(result.returncode, 2)
        self.assertIn("must be a positive integer", result.stderr)

    def test_rejects_malformed_paired_reads(self) -> None:
        result = self.run_cli("-a", "assembly.fa", "-s", "reads_R1.fastq")
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires R1,R2", result.stderr)

    def test_rejects_a_normal_run_without_an_assembly(self) -> None:
        result = self.run_cli()
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires at least one assembly", result.stderr)

    def test_external_binset_rejects_autobinning_module(self) -> None:
        result = self.run_cli(
            "-r", "existing_bins",
            "-a", "1_assembly.fa",
            "-c", "Coverage_matrix_for_binning_1_assembly.fa.txt",
            "-s", "PE_r1_x.fastq,PE_r2_x.fastq",
            "--module", "autobinning",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("cannot be combined with --refinement-binset", result.stderr)

    def test_external_binset_gap_filling_reports_missing_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = subprocess.run(
                [
                    sys.executable, str(MODULES / "BASALT.py"),
                    "-r", "existing_bins",
                    "-a", "1_assembly.fa",
                    "-c", "Coverage_matrix_for_binning_1_assembly.fa.txt",
                    "-s", "PE_r1_x.fastq,PE_r2_x.fastq",
                    "--module", "reassembly",
                ],
                cwd=temporary,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(result.returncode, 2)
            self.assertTrue(
                "Missing assembly file" in result.stderr
                or "refinement binset does not exist" in result.stderr,
                result.stderr,
            )


class CheckpointParsingTests(unittest.TestCase):
    def test_ignores_non_numeric_lines_and_keeps_the_file(self) -> None:
        """The plain "Outlier removal done!" line must not reset progress."""
        with tempfile.TemporaryDirectory() as temporary:
            checkpoint = Path(temporary) / "Basalt_checkpoint.txt"
            checkpoint.write_text(
                "1st autobinner done!\n"
                "2nd bin selection within group done!\n"
                "3rd bin selection within multiple groups done!\n"
                "4th outlier removal done!\n"
                "Outlier removal done!\n",
                encoding="utf-8",
            )
            self.assertEqual(read_checkpoint_step(str(checkpoint)), 4)
            self.assertIn("Outlier removal done!", checkpoint.read_text(encoding="utf-8"))

    def test_missing_file_means_step_zero(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            self.assertEqual(
                read_checkpoint_step(str(Path(temporary) / "absent.txt")), 0
            )

    def test_append_checkpoint_adds_separator(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            checkpoint = Path(temporary) / "Basalt_checkpoint.txt"
            checkpoint.write_text("4th outlier removal done!", encoding="utf-8")
            append_checkpoint(str(checkpoint), "5th contig retrieve done!")
            content = checkpoint.read_text(encoding="utf-8")
            self.assertEqual(
                content, "4th outlier removal done!\n5th contig retrieve done!\n"
            )
            self.assertEqual(read_checkpoint_step(str(checkpoint)), 5)


class DepthHeaderTests(unittest.TestCase):
    def test_reads_legacy_header_with_bam_suffix(self) -> None:
        header = (
            "contigName\tcontigLen\tcontigDepth\t"
            "sample.bam-totalAvgDepth\tsample.bam-var"
        )
        self.assertEqual(depth_sample_count(header), 1)

    def test_reads_new_header_without_bam_suffix(self) -> None:
        """jgi_summarize_bam_contig_depths v2.18+ drops .bam from names."""
        header = (
            "contigName\tcontigLen\tcontigDepth\t"
            "sample1-totalAvgDepth\tsample1-var\t"
            "sample2-totalAvgDepth\tsample2-var"
        )
        self.assertEqual(depth_sample_count(header), 2)

    def test_rejects_an_unrecognised_header(self) -> None:
        with self.assertRaisesRegex(ValueError, "sample count"):
            depth_sample_count("contigName\tcontigLen\tcontigDepth")


class CountParsingTests(unittest.TestCase):
    def test_accepts_float_formatted_genome_sizes(self) -> None:
        self.assertEqual(coerce_count("1860636.0"), 1860636)
        self.assertEqual(coerce_count(" 42 "), 42)
        self.assertEqual(coerce_count(7), 7)

    def test_rejects_non_numeric_values(self) -> None:
        with self.assertRaises(ValueError):
            coerce_count("n/a")


class OlcRuntimeTests(unittest.TestCase):
    def test_thread_budget_splits_allocation(self) -> None:
        self.assertEqual(worker_thread_budget(24, 6), 4)
        self.assertEqual(worker_thread_budget(1, 8), 1)
        self.assertEqual(worker_thread_budget("not-a-number", 4), 1)

    def test_pe_read_names_are_idempotent(self) -> None:
        """External-binset runs pass the data-feeding PE files directly."""
        self.assertEqual(
            pe_read_names(["r1.fq", "r2.fq"]),
            ["PE_r1_r1.fq", "PE_r2_r2.fq"],
        )
        self.assertEqual(
            pe_read_names(["PE_r1_r1.fq", "PE_r2_r2.fq"]),
            ["PE_r1_r1.fq", "PE_r2_r2.fq"],
        )

    def test_alignment_group_expansion_converges(self) -> None:
        group = {"a": 1}
        alignments = {"a b": 0, "b c": 0, "x y": 0}
        expanded, converged = expand_alignment_group(group, alignments)
        self.assertTrue(converged)
        self.assertEqual(sorted(expanded), ["a", "b", "c"])

    def test_fasta_has_sequence_detects_empty_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            empty = Path(temporary) / "empty.fa"
            empty.write_text(">bin.1\n", encoding="utf-8")
            filled = Path(temporary) / "filled.fa"
            filled.write_text(">bin.1\nACGT\n", encoding="utf-8")
            self.assertFalse(fasta_has_sequence(empty))
            self.assertTrue(fasta_has_sequence(filled))
            self.assertFalse(fasta_has_sequence(Path(temporary) / "absent.fa"))

    def test_stale_olc_cleanup_keeps_result_folders(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            stale = [
                work / "bin.1.fa_merged",
                work / "bin.1.fa_checkm",
                work / "blast_bin.1.fa_self_merged_1.txt",
            ]
            for index, path in enumerate(stale):
                if path.suffix == ".txt":
                    path.write_text("hit\n", encoding="utf-8")
                else:
                    path.mkdir()
                    (path / "part.fa").write_text(">x\nA\n", encoding="utf-8")
            results = work / "target_OLC"
            results.mkdir()
            (results / "bin.1.fa").write_text(">b\nACGT\n", encoding="utf-8")

            clean_stale_olc_intermediates(str(work), "bin.1.fa")

            for path in stale:
                self.assertFalse(path.exists(), path)
            self.assertTrue(results.is_dir())


class ExternalWorkflowTests(unittest.TestCase):
    def make_workdir(self, temporary, assemblies=("500_binner_A", "501_binner_B")):
        work = Path(temporary)
        for name in assemblies:
            (work / (name + "_assembly.fa")).write_text(">c\nACGT\n", encoding="utf-8")
            (work / ("condense_connections_" + name + ".txt")).write_text(
                "node1\tinteraction\tnode2\tconnections\n", encoding="utf-8"
            )
            (work / ("Coverage_matrix_for_binning_" + name + "_assembly.fa.txt")).write_text(
                "Name\tLength\n", encoding="utf-8"
            )
        binset = work / "BestBinset"
        binset.mkdir()
        (binset / "bin.1.fa").write_text(">c\nACGT\n", encoding="utf-8")
        comparison = work / "BestBinset_comparison_files"
        comparison.mkdir()
        (comparison / "Selected_bins_test.txt").write_text(
            "Bin1\tbin.1---bin.2\t{'Completeness': 98} {'Contamination': 0.4}\n",
            encoding="utf-8",
        )
        return work

    def test_module_plan(self) -> None:
        self.assertEqual(external_module_plan(False, "all"), "outlier")
        self.assertEqual(external_module_plan(True, "refinement"), "refinement")
        self.assertEqual(external_module_plan(True, "reassembly"), "gap-filling")
        self.assertEqual(external_module_plan(True, "all"), "gap-filling")
        with self.assertRaises(ValueError):
            external_module_plan(True, "autobinning")

    def test_connection_pairing_by_stem_and_index(self) -> None:
        assemblies = ["500_binner_A_assembly.fa", "501_binner_B_assembly.fa"]
        connections = [
            "condense_connections_501_binner_B.txt",
            "condense_connections_500_binner_A.txt",
        ]
        paired = pair_connections(assemblies, connections)
        self.assertEqual(
            paired,
            [
                "condense_connections_500_binner_A.txt",
                "condense_connections_501_binner_B.txt",
            ],
        )

    def test_prepare_route_writes_state_and_promotes_outlier_step(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            work = self.make_workdir(temporary)
            (work / "Basalt_checkpoint.txt").write_text(
                "3rd external binset supplied; autobinning skipped\n",
                encoding="utf-8",
            )
            (work / "BestBinset_outlier_refined").mkdir()

            state = prepare_external_binset_route(
                "BestBinset",
                ["500_binner_A_assembly.fa", "501_binner_B_assembly.fa"],
                [
                    "Coverage_matrix_for_binning_500_binner_A_assembly.fa.txt",
                    "Coverage_matrix_for_binning_501_binner_B_assembly.fa.txt",
                ],
                "reassembly",
                "continue",
                workdir=str(work),
            )

            self.assertEqual(state["orchestrator_module"], "all")
            self.assertEqual(state["checkpoint_step"], 4)
            self.assertEqual(
                (work / "Coverage_matrix_list.txt").read_text(encoding="utf-8").split(),
                [
                    "Coverage_matrix_for_binning_500_binner_A_assembly.fa.txt",
                    "Coverage_matrix_for_binning_501_binner_B_assembly.fa.txt",
                ],
            )
            self.assertEqual(
                (work / "Bestbinset_list.txt").read_text(encoding="utf-8").split(),
                ["BestBinset", "BestBinset"],
            )
            self.assertEqual(
                read_checkpoint_step(str(work / "Basalt_checkpoint.txt")), 4
            )

    def test_prepare_route_requires_connections(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            work = self.make_workdir(temporary)
            (work / "condense_connections_500_binner_A.txt").unlink()
            with self.assertRaisesRegex(ExternalBinsetError, "condense_connections"):
                prepare_external_binset_route(
                    "BestBinset",
                    ["500_binner_A_assembly.fa", "501_binner_B_assembly.fa"],
                    [
                        "Coverage_matrix_for_binning_500_binner_A_assembly.fa.txt",
                        "Coverage_matrix_for_binning_501_binner_B_assembly.fa.txt",
                    ],
                    "refinement",
                    "continue",
                    workdir=str(work),
                )

    def test_prepare_route_requires_comparison_files_for_multi_assembly_gap_filling(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            work = self.make_workdir(temporary)
            import shutil as _shutil

            _shutil.rmtree(work / "BestBinset_comparison_files")
            with self.assertRaisesRegex(ExternalBinsetError, "BestBinset_comparison_files"):
                prepare_external_binset_route(
                    "BestBinset",
                    ["500_binner_A_assembly.fa", "501_binner_B_assembly.fa"],
                    [
                        "Coverage_matrix_for_binning_500_binner_A_assembly.fa.txt",
                        "Coverage_matrix_for_binning_501_binner_B_assembly.fa.txt",
                    ],
                    "reassembly",
                    "continue",
                    workdir=str(work),
                )


class CheckM2CacheTests(unittest.TestCase):
    """Batched CheckM2 evaluation with content-keyed reuse (issue 99)."""

    def make_bins(self, temporary):
        root = Path(temporary)
        folder = root / "merged"
        folder.mkdir()
        (folder / "bin.1.fa").write_text(">a\n" + "ACGT" * 10 + "\n", encoding="utf-8")
        (folder / "bin.2.fa").write_text(">b\n" + "TTTT" * 10 + "\n", encoding="utf-8")
        return root, folder

    def fake_predict(self, calls):
        def _predict(command, cwd=None):
            calls.append([str(item) for item in command])
            output = Path(command[command.index("-o") + 1]) / "quality_report.tsv"
            output.parent.mkdir(parents=True, exist_ok=True)
            with open(output, "w", encoding="utf-8") as handle:
                handle.write(
                    "Name\tCompleteness\tContamination\tModel\tTable\tDensity\t"
                    "Contig_N50\tScaffolds\tGenome_Size\n"
                )
                for name, cpn, ctn, n50, size in (
                    ("bin.1", 98.5, 0.42, 40000, 1860636.0),
                    ("bin.2", 71.2, 3.10, 12000, 2555704),
                ):
                    handle.write(
                        "{}\t{}\t{}\tNeural Network\t11\t0.89\t{}\t3\t{}\n".format(
                            name, cpn, ctn, n50, size
                        )
                    )

        return _predict

    def test_first_call_evaluates_and_second_reuses(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, folder = self.make_bins(temporary)
            cache = root / "Checkm2_metrics_cache.tsv"
            calls = []
            messages = []
            with patch("basalt_runtime._run", side_effect=self.fake_predict(calls)), patch(
                "basalt_runtime._checkm2_signature", return_value="checkm2|test"
            ):
                first = evaluate_bins_checkm2(
                    str(folder), 8, cache_path=str(cache), notify=messages.append
                )
                second = evaluate_bins_checkm2(
                    str(folder), 8, cache_path=str(cache), notify=messages.append
                )

            self.assertEqual(len(calls), 1)
            self.assertIn("-t", calls[0])
            self.assertEqual(calls[0][calls[0].index("-t") + 1], "8")
            self.assertEqual(first["bin.1"]["Completeness"], 98.5)
            self.assertEqual(first["bin.1"]["Genome size"], 1860636)
            self.assertEqual(first["bin.2"]["N50"], 12000)
            self.assertEqual(second, first)
            self.assertTrue(any("reused, 0 evaluated" in m for m in messages))
            self.assertFalse(list(root.glob("Checkm2_batch_eval_*")))

    def test_header_changes_do_not_invalidate_the_cache(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, folder = self.make_bins(temporary)
            cache = root / "Checkm2_metrics_cache.tsv"
            (folder / "bin.2.fa").write_text(
                ">renamed_header\n" + "TTTT" * 10 + "\n", encoding="utf-8"
            )
            calls = []
            with patch("basalt_runtime._run", side_effect=self.fake_predict(calls)), patch(
                "basalt_runtime._checkm2_signature", return_value="checkm2|test"
            ):
                first = evaluate_bins_checkm2(str(folder), 4, cache_path=str(cache))
                second = evaluate_bins_checkm2(str(folder), 4, cache_path=str(cache))
            self.assertEqual(len(calls), 1)
            self.assertEqual(second, first)

    def test_version_change_resets_the_cache(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, folder = self.make_bins(temporary)
            cache = root / "Checkm2_metrics_cache.tsv"
            calls = []
            with patch("basalt_runtime._run", side_effect=self.fake_predict(calls)), patch(
                "basalt_runtime._checkm2_signature", return_value="checkm2|v1"
            ):
                evaluate_bins_checkm2(str(folder), 4, cache_path=str(cache))
            with patch("basalt_runtime._run", side_effect=self.fake_predict(calls)), patch(
                "basalt_runtime._checkm2_signature", return_value="checkm2|v2"
            ):
                evaluate_bins_checkm2(str(folder), 4, cache_path=str(cache))
            self.assertEqual(len(calls), 2)

    def test_predict_failure_returns_none_for_the_direct_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, folder = self.make_bins(temporary)
            with patch(
                "basalt_runtime._run", side_effect=RuntimeError("checkm2 missing")
            ), patch(
                "basalt_runtime._checkm2_signature", return_value="checkm2|test"
            ):
                result = evaluate_bins_checkm2(str(folder), 4, cache_path=str(root / "c.tsv"))
            self.assertIsNone(result)

    def test_content_hash_ignores_headers_and_case(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            a = root / "a.fa"
            b = root / "b.fa"
            a.write_text(">x\nACGT\nACgt\n", encoding="utf-8")
            b.write_text(">totally-different\nacgt acgt\n", encoding="utf-8")
            self.assertEqual(fasta_content_hash(a), fasta_content_hash(b))

    def test_version_token_extraction_tolerates_banners(self) -> None:
        from basalt_runtime import _extract_version_token

        self.assertEqual(_extract_version_token("1.1.0\n"), "1.1.0")
        self.assertEqual(
            _extract_version_token("WARNING: database location\n1.1.0\n"), "1.1.0"
        )
        self.assertEqual(
            _extract_version_token("Running CheckM2 version 1.1.0\n"), "1.1.0"
        )
        self.assertEqual(_extract_version_token("v1.1.0"), "1.1.0")
        self.assertIsNone(_extract_version_token(""))
        self.assertIsNone(_extract_version_token("some failure\n"))


if __name__ == "__main__":
    unittest.main()
