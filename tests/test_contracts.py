"""Contract tests for the parts of the pipeline that are pure logic.

These run without reference data, without conda environments and without a
cluster, so they can gate every commit. Each test corresponds to a defect that
reached a reviewer:

* the metadata schema, because the old parser guessed groups by stripping
  ``_vs_Input`` from file names and therefore could not represent a second
  condition;
* ``resolve_mode``, because a single-condition design silently produced
  "differential" numbers;
* the gene assignment, because distances were measured to a promoter peak
  instead of to a transcription start site, and because minus-strand TSSs were
  taken as the interval start;
* the annotation registry and the ``<assembly>-<source>`` path segment, because a
  wrong resolution writes one genome's results into another genome's directory;
* stage 1 linking and the Picard metrics parser, because ChIP Input libraries are
  not metadata rows and because blank lines before the metrics header aborted QC;
* ``workspace.marker`` and the dashboard's module reload, because a long-running
  panel kept the old signature and every results page load raised ``TypeError``;
* the signal-tab interval clipping and peak merge, because a fixed 1 Mb opening
  fell outside the mock's 60 kb contigs and stacked narrowPeak rows painted one
  locus several times over;
* the launcher's "all stages" target, because it meant "the last row in the
  table" and asked snakemake for 5_diffbind on a single-condition design, which
  the Snakefile never compiled in - the run died with MissingRuleException
  seconds after launch;
* the naming table, because a fixed-row editor drafted a row for the Input
  library that is already cited as input_id, and there was no way to delete it;
* the annotation page's metric block, because it described the annotation
  already written to config, so picking GRCh38 still showed the mouse TxDb and
  the peak annotation library never followed the selection;
* the annotation build panel, because a build halfway through ``hisat2-build``
  looked exactly like one that had finished - one status line and a log tail,
  no stage, no percentage, nothing that moved;
* the annotation file list, because it took the last component of
  ``rna_index`` (``genome``) as the index directory, so a complete GRCm38
  build was reported as five missing files under a ``genome/genome.N.ht2``
  path that never exists;
* the download panel, because between "started" and "done" it showed a log
  tail and nothing that moved, and its button waited for a FASTQ set that only
  a finished download could produce.

Run with: python -m unittest discover -s tests -v
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

SCRIPT_ROOT = os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts")
sys.path.insert(0, SCRIPT_ROOT)

from common import annotate, annotation, metadata, workspace
from common.metadata import comparisons, resolve_mode

SIX_COLUMN = """sample,mark,condition,replicate,input_id,fastq_prefix
A1,H3K27ac,Veh,1,IN_Veh,
A2,H3K27ac,Veh,2,IN_Veh,
B1,H3K27ac,TGFb,1,IN_TGFb,
B2,H3K27ac,TGFb,2,IN_TGFb,
"""

THREE_COLUMN = """IP sample,Input,IP_name
SRR1,SRR0,H3K27ac_rep1
SRR2,SRR0,H3K27ac_rep2
"""

SINGLE_CONDITION = """sample,mark,condition,replicate,input_id,fastq_prefix
A1,H3K27ac,single,1,IN,
A2,H3K27ac,single,2,IN,
"""

GROUP_SCHEMA = """sample,sample_name,condition,input_sample,input_sample_name
A1,H3K27ac_veh,rep1,IN1,IN_veh
A2,H3K27ac_veh,rep2,IN1,IN_veh
B1,H3K27ac_tgfb,rep1,IN2,IN_tgfb
B2,H3K27ac_tgfb,rep2,IN2,IN_tgfb
"""

ROLE_SCHEMA = """sample,role,sample_name,condition,input_sample
IN_veh_1,Input,input_veh,rep1,
IN_veh_2,Input,input_veh,rep2,
A1,IP,H3K27ac_veh,rep1,IN_veh_1
A2,IP,H3K27ac_veh,rep2,IN_veh_2
B1,IP,H3K27ac_tgfb,rep1,IN_veh_1
B2,IP,H3K27ac_tgfb,rep2,IN_veh_2
"""

COMPARISONS = """reference,treatment
H3K27ac_veh,H3K27ac_tgfb
"""


def write(text: str) -> str:
    handle = tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False)
    handle.write(text)
    handle.close()
    return handle.name


class MetadataSchemaTest(unittest.TestCase):
    def test_six_column_schema_is_parsed_verbatim(self):
        samples = metadata.read_chip_metadata(write(SIX_COLUMN))
        self.assertEqual(len(samples), 4)
        first = samples[0]
        self.assertEqual((first.sample, first.mark, first.condition), ("A1", "H3K27ac", "Veh"))
        self.assertEqual(first.group, "H3K27ac_Veh")
        self.assertEqual(first.input_id, "IN_Veh")
        self.assertEqual(first.replicate, 1)

    def test_legacy_schema_forces_single_condition(self):
        samples = metadata.read_chip_metadata(write(THREE_COLUMN))
        self.assertTrue(all(s.is_single_condition for s in samples))
        self.assertEqual({s.mark for s in samples}, {"H3K27ac"})
        self.assertEqual([s.replicate for s in samples], [1, 2])

    def test_legacy_schema_is_one_group(self):
        samples = metadata.read_chip_metadata(write(THREE_COLUMN))
        self.assertEqual(metadata.groups(samples), ["H3K27ac_single"])

    def test_duplicate_replicates_are_rejected(self):
        broken = SIX_COLUMN.replace("A2,H3K27ac,Veh,2", "A2,H3K27ac,Veh,1")
        with self.assertRaises(ValueError):
            metadata.read_chip_metadata(write(broken))

    def test_duplicate_sample_ids_are_rejected(self):
        broken = SIX_COLUMN.replace("B1,H3K27ac", "A1,H3K27ac")
        with self.assertRaises(ValueError):
            metadata.read_chip_metadata(write(broken))

    def test_missing_mark_column_is_rejected(self):
        with self.assertRaises(ValueError):
            metadata.read_chip_metadata(write("sample,condition,replicate,input_id\nX,Veh,1,IN\n"))


class ModeResolutionTest(unittest.TestCase):
    def test_two_conditions_in_one_mark_give_differential(self):
        samples = metadata.read_chip_metadata(write(SIX_COLUMN))
        self.assertEqual(resolve_mode({}, samples), "differential")

    def test_single_condition_gives_descriptive(self):
        samples = metadata.read_chip_metadata(write(SINGLE_CONDITION))
        self.assertEqual(resolve_mode({}, samples), "descriptive")

    def test_three_column_schema_gives_descriptive(self):
        samples = metadata.read_chip_metadata(write(THREE_COLUMN))
        self.assertEqual(resolve_mode({}, samples), "descriptive")

    def test_explicit_descriptive_wins_over_available_design(self):
        samples = metadata.read_chip_metadata(write(SIX_COLUMN))
        config = {"integrated": {"mode": "descriptive"}}
        self.assertEqual(resolve_mode(config, samples), "descriptive")

    def test_one_replicate_per_condition_gives_descriptive(self):
        samples = metadata.read_chip_metadata(
            write(
                "sample,mark,condition,replicate,input_id\n"
                "A1,H3K27ac,Veh,1,IN1\n"
                "B1,H3K27ac,TGFb,1,IN2\n"
            )
        )
        self.assertEqual(resolve_mode({}, samples), "descriptive")

    def test_group_names_need_a_comparison_file(self):
        samples = metadata.read_chip_metadata(write(GROUP_SCHEMA))
        self.assertEqual(resolve_mode({}, samples), "descriptive")
        path = write(COMPARISONS)
        try:
            self.assertEqual(resolve_mode({}, samples, path), "differential")
        finally:
            os.unlink(path)

    def test_comparisons_use_the_declared_reference(self):
        samples = metadata.read_chip_metadata(write(SIX_COLUMN))
        config = {"chip": {"reference_condition": "Veh", "diff": {"min_replicates": 2}}}
        self.assertEqual(comparisons(samples, config), [("H3K27ac", "Veh", "TGFb")])

    def test_comparisons_default_to_alphabetical_reference(self):
        samples = metadata.read_chip_metadata(write(SIX_COLUMN))
        self.assertEqual(comparisons(samples, {}), [("H3K27ac", "TGFb", "Veh")])

    def test_bad_reference_condition_is_an_error(self):
        samples = metadata.read_chip_metadata(write(SIX_COLUMN))
        config = {"chip": {"reference_condition": "Mock"}}
        with self.assertRaises(ValueError):
            comparisons(samples, config)

    def test_no_comparison_without_replicates(self):
        thin = SIX_COLUMN.replace("A2,H3K27ac,Veh,2,IN_Veh,\n", "")
        samples = metadata.read_chip_metadata(write(thin))
        self.assertEqual(comparisons(samples, {}), [])

    def test_comparison_name_is_readable(self):
        self.assertEqual(
            metadata.comparison_name("H3K27ac", "Veh", "TGFb"), "H3K27ac_TGFb_vs_Veh"
        )


class GroupSchemaTest(unittest.TestCase):
    """The five-column schema: sample_name is the group, condition is the repeat."""

    def test_group_names_are_read_verbatim(self):
        first = metadata.read_chip_metadata(write(GROUP_SCHEMA))[0]
        self.assertEqual(first.sample, "A1")
        self.assertEqual(first.sample_name, "H3K27ac_veh")
        self.assertEqual(first.group, "H3K27ac_veh")
        self.assertEqual(first.comparison_key, "H3K27ac_veh")
        self.assertEqual(first.condition, metadata.SINGLE_CONDITION)
        self.assertEqual(first.replicate, 1)
        self.assertEqual(first.input_id, "IN1")
        self.assertEqual(first.input_name, "IN_veh")

    def test_a_blank_input_name_falls_back_to_the_input_sample(self):
        schema = GROUP_SCHEMA.replace(
            "A1,H3K27ac_veh,rep1,IN1,IN_veh", "A1,H3K27ac_veh,rep1,IN1,"
        )
        self.assertEqual(metadata.read_chip_metadata(write(schema))[0].input_name, "IN1")

    def test_groups_are_the_sample_names(self):
        samples = metadata.read_chip_metadata(write(GROUP_SCHEMA))
        self.assertEqual(metadata.groups(samples), ["H3K27ac_tgfb", "H3K27ac_veh"])

    def test_a_repeated_rep_inside_one_group_is_rejected(self):
        broken = GROUP_SCHEMA.replace("A2,H3K27ac_veh,rep2", "A2,H3K27ac_veh,rep1")
        with self.assertRaises(ValueError):
            metadata.read_chip_metadata(write(broken))

    def test_a_missing_group_name_is_rejected(self):
        broken = GROUP_SCHEMA.replace("A1,H3K27ac_veh,rep1", "A1,,rep1")
        with self.assertRaises(ValueError):
            metadata.read_chip_metadata(write(broken))

    def test_a_condition_outside_rep1_to_rep3_is_rejected(self):
        broken = GROUP_SCHEMA.replace("A1,H3K27ac_veh,rep1", "A1,H3K27ac_veh,Veh")
        with self.assertRaises(ValueError):
            metadata.read_chip_metadata(write(broken))

    def test_group_names_are_recognised_as_such(self):
        self.assertTrue(
            metadata.uses_group_names(metadata.read_chip_metadata(write(GROUP_SCHEMA)))
        )
        self.assertFalse(
            metadata.uses_group_names(metadata.read_chip_metadata(write(SIX_COLUMN)))
        )

    def test_the_mark_is_the_part_before_the_last_underscore(self):
        self.assertEqual(metadata.display_mark("H3K27ac_veh"), "H3K27ac")
        self.assertEqual(metadata.display_mark("H3K27ac"), "H3K27ac")


class RoleSchemaTest(unittest.TestCase):
    """The role-aware schema: one row per library, Input included."""

    def test_input_rows_are_not_contrast_members(self):
        samples = metadata.read_chip_metadata(write(ROLE_SCHEMA))
        self.assertEqual([s.sample for s in samples], ["A1", "A2", "B1", "B2"])
        self.assertTrue(all(s.role == "IP" for s in samples))

    def test_input_rows_are_read_as_their_own_libraries(self):
        inputs = metadata.read_chip_inputs(write(ROLE_SCHEMA))
        self.assertEqual([s.sample for s in inputs], ["IN_veh_1", "IN_veh_2"])
        self.assertEqual({s.sample_name for s in inputs}, {"input_veh"})
        self.assertEqual([s.replicate for s in inputs], [1, 2])

    def test_the_input_name_comes_from_the_input_row(self):
        samples = metadata.read_chip_metadata(write(ROLE_SCHEMA))
        self.assertEqual({s.input_name for s in samples}, {"input_veh"})

    def test_two_input_libraries_can_pair_with_their_own_replicate(self):
        samples = metadata.read_chip_metadata(write(ROLE_SCHEMA))
        self.assertEqual(samples[0].input_id, "IN_veh_1")
        self.assertEqual(samples[1].input_id, "IN_veh_2")

    def test_a_blank_input_replicate_is_numbered(self):
        schema = ROLE_SCHEMA.replace("IN_veh_1,Input,input_veh,rep1", "IN_veh_1,Input,input_veh,")
        schema = schema.replace("IN_veh_2,Input,input_veh,rep2", "IN_veh_2,Input,input_veh,")
        inputs = metadata.read_chip_inputs(write(schema))
        self.assertEqual([s.replicate for s in inputs], [1, 2])

    def test_an_ip_citing_an_undeclared_input_is_rejected(self):
        broken = ROLE_SCHEMA.replace("rep1,IN_veh_1\nA2", "rep1,IN_missing\nA2")
        with self.assertRaises(ValueError):
            metadata.read_chip_metadata(write(broken))

    def test_two_inputs_of_one_group_cannot_share_a_replicate(self):
        broken = ROLE_SCHEMA.replace("IN_veh_2,Input,input_veh,rep2", "IN_veh_2,Input,input_veh,rep1")
        with self.assertRaises(ValueError):
            metadata.read_chip_metadata(write(broken))

    def test_an_input_row_cannot_cite_another_input(self):
        broken = ROLE_SCHEMA.replace("IN_veh_1,Input,input_veh,rep1,", "IN_veh_1,Input,input_veh,rep1,IN_veh_2")
        with self.assertRaises(ValueError):
            metadata.read_chip_metadata(write(broken))

    def test_an_older_file_declares_no_input_rows(self):
        self.assertIsNone(metadata.read_chip_inputs(write(GROUP_SCHEMA)))

    def test_an_unknown_role_is_rejected(self):
        broken = ROLE_SCHEMA.replace("A1,IP,", "A1,IPx,")
        with self.assertRaises(ValueError):
            metadata.read_chip_metadata(write(broken))


class ComparisonFileTest(unittest.TestCase):
    """comparisons.csv decides which group is contrasted against which."""

    def test_pairs_are_read_in_file_order(self):
        path = write(COMPARISONS)
        try:
            self.assertEqual(
                metadata.read_comparisons(path), [("H3K27ac_veh", "H3K27ac_tgfb")]
            )
        finally:
            os.unlink(path)

    def test_half_empty_self_and_repeated_pairs_are_dropped(self):
        path = write("reference,treatment\nA,\n,B\nA,A\nA,B\nA,B\n")
        try:
            self.assertEqual(metadata.read_comparisons(path), [("A", "B")])
        finally:
            os.unlink(path)

    def test_a_missing_file_declares_no_comparison(self):
        with tempfile.TemporaryDirectory() as root:
            self.assertEqual(
                metadata.read_comparisons(os.path.join(root, "comparisons.csv")), []
            )

    def test_the_treatment_names_the_output_namespace(self):
        samples = metadata.read_chip_metadata(write(GROUP_SCHEMA))
        path = write(COMPARISONS)
        try:
            self.assertEqual(
                comparisons(samples, {}, path),
                [("H3K27ac_tgfb", "H3K27ac_veh", "H3K27ac_tgfb")],
            )
        finally:
            os.unlink(path)

    def test_a_group_with_one_replicate_is_skipped(self):
        thin = GROUP_SCHEMA.replace("A2,H3K27ac_veh,rep2,IN1,IN_veh\n", "")
        samples = metadata.read_chip_metadata(write(thin))
        path = write(COMPARISONS)
        try:
            self.assertEqual(comparisons(samples, {}, path), [])
        finally:
            os.unlink(path)

    def test_both_groups_enter_one_contrast_reference_first(self):
        samples = metadata.read_chip_metadata(write(GROUP_SCHEMA))
        picked = metadata.comparison_samples(samples, "H3K27ac_veh", "H3K27ac_tgfb")
        self.assertEqual([s.sample for s in picked], ["A1", "A2", "B1", "B2"])

    def test_the_legacy_schema_still_filters_by_mark(self):
        samples = metadata.read_chip_metadata(write(SIX_COLUMN))
        picked = metadata.comparison_samples(samples, "Veh", "TGFb", "H3K27ac")
        self.assertEqual([s.sample for s in picked], ["A1", "A2", "B1", "B2"])


class GeneAssignmentTest(unittest.TestCase):
    """geneA TSS = 1000 (plus strand); geneB TSS = 50999 (minus strand end-1)."""

    def setUp(self):
        path = write(
            "chr1\t1000\t2000\tgeneA\t0\t+\tgeneA\n"
            "chr1\t50000\t51000\tgeneB\t0\t-\tgeneB\n"
        )
        self.index = annotate.load_genes(path)
        os.unlink(path)

    def test_promoter_window_uses_tss_not_peak_start(self):
        symbol, distance, region_type = annotate.assign_gene(
            self.index, "chr1", 50900, 51000, 2000, 1000000
        )
        self.assertEqual((symbol, region_type), ("geneB", "promoter"))
        self.assertEqual(distance, 0)

    def test_promoter_window_is_inclusive(self):
        symbol, _distance, region_type = annotate.assign_gene(
            self.index, "chr1", 3000, 3100, 2000, 1000000
        )
        self.assertEqual((symbol, region_type), ("geneA", "promoter"))

    def test_distal_region_reports_distance_to_tss(self):
        symbol, distance, region_type = annotate.assign_gene(
            self.index, "chr1", 11000, 11100, 2000, 1000000
        )
        self.assertEqual((symbol, region_type), ("geneA", "distal"))
        self.assertEqual(distance, 10000)

    def test_region_beyond_max_distance_is_unassigned(self):
        symbol, _distance, region_type = annotate.assign_gene(
            self.index, "chr1", 900000, 900100, 2000, 100000
        )
        self.assertEqual((symbol, region_type), ("", "unassigned"))

    def test_unknown_contig_is_unassigned(self):
        symbol, _distance, region_type = annotate.assign_gene(
            self.index, "chrZ", 10, 20, 2000, 1000000
        )
        self.assertEqual((symbol, region_type), ("", "unassigned"))

    def test_interval_distance_semantics(self):
        self.assertEqual(annotate.interval_distance(100, 200, 150, 250), 0)
        self.assertEqual(annotate.interval_distance(100, 200, 300, 400), 100)
        self.assertEqual(annotate.interval_distance(300, 400, 100, 200), 100)


class GeneTest(unittest.TestCase):
    def test_strand_aware_tss(self):
        plus = annotate.Gene("chr1", 100, 200, "a", "a", "+")
        minus = annotate.Gene("chr1", 100, 200, "a", "a", "-")
        self.assertEqual(plus.tss, 100)
        self.assertEqual(minus.tss, 199)


class MergePeaksFailFastTest(unittest.TestCase):
    """A missing replicate peak bed must abort, not shrink the group.

    The previous pipeline wrote a placeholder and exited 0, so a run with
    silently missing replicates looked successful. This test pins the exit
    status and the refusal to write an output.
    """

    METADATA = (
        "sample,mark,condition,replicate,input_id\n"
        "A1,H3K27ac,Veh,1,IN\n"
        "A2,H3K27ac,Veh,2,IN\n"
    )

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.metadata = os.path.join(self.tmp, "meta.csv")
        with open(self.metadata, "w") as handle:
            handle.write(self.METADATA)
        self.peaks = os.path.join(self.tmp, "peaks")
        os.makedirs(self.peaks)
        self.out = os.path.join(self.tmp, "out")
        bindir = os.path.join(self.tmp, "bin")
        os.makedirs(bindir)
        stub = os.path.join(bindir, "bedtools")
        with open(stub, "w") as handle:
            handle.write("#!/bin/sh\nprintf 'chr1\\t100\\t200\\tA1,A2\\t2\\n'\n")
        os.chmod(stub, 0o755)
        self.env = dict(os.environ, PATH=bindir + os.pathsep + os.environ["PATH"])
        self.script = os.path.join(SCRIPT_ROOT, "ChIPseq_pipline", "4_5_merge_peaks.py")
        self.merged = os.path.join(self.out, "H3K27ac_Veh", "H3K27ac_Veh_peaks_merged.bed")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self):
        return subprocess.run(
            [
                sys.executable, self.script,
                "--metadata", self.metadata,
                "--peakdir", self.peaks,
                "--outdir", self.out,
            ],
            capture_output=True, text=True, env=self.env,
        )

    def _write_peaks(self, samples):
        for sample in samples:
            directory = os.path.join(self.peaks, sample)
            os.makedirs(directory, exist_ok=True)
            with open(os.path.join(directory, f"{sample}_peaks.bed"), "w") as handle:
                handle.write("chr1\t100\t200\tpeak\t0\t.\n")

    def test_missing_replicate_bed_exits_nonzero(self):
        self._write_peaks(["A1"])
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing replicate peak beds", result.stderr)
        self.assertFalse(os.path.exists(self.merged), "a failed merge must not leave an output")

    def test_supporting_replicate_count_is_written(self):
        self._write_peaks(["A1", "A2"])
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        with open(self.merged) as handle:
            fields = handle.readline().rstrip("\n").split("\t")
        self.assertEqual(fields[:3], ["chr1", "100", "200"])
        self.assertEqual(fields[3], "H3K27ac_Veh_region_1")
        self.assertEqual(fields[4], "2")


class AnnotationRegistryTest(unittest.TestCase):
    """The annotation version and source decide every path and every TxDb.

    Results are stored per ``<assembly>-<source>``, so a wrong resolution does
    not fail loudly: it silently writes GRCh37 output into a GRCh38 directory.
    """

    def test_hg19_alias_resolves_to_grch37(self):
        info = annotation.single("hg19")
        self.assertEqual(info["assembly"], "GRCh37")
        self.assertEqual(info["id"], "GRCh37-gencode")
        self.assertEqual(info["txdb"], "TxDb.Hsapiens.UCSC.hg19.knownGene")
        self.assertEqual(info["orgdb"], "org.Hs.eg.db")

    def test_grch37_and_grch38_stay_apart(self):
        grch37 = annotation.single("GRCh37")
        grch38 = annotation.single("GRCh38")
        self.assertNotEqual(grch37["id"], grch38["id"])
        self.assertNotEqual(grch37["genome"], grch38["genome"])
        self.assertNotEqual(grch37["txdb"], grch38["txdb"])

    def test_unknown_assembly_is_rejected(self):
        with self.assertRaises(ValueError):
            annotation.single("hg18")

    def test_unknown_source_is_rejected(self):
        with self.assertRaises(ValueError):
            annotation.single("GRCh37", "refseq")

    def test_source_defaults_to_the_only_source(self):
        self.assertEqual(annotation.declared({"annotation": {"assembly": "hg19"}})[1], "gencode")


class ResultLayoutTest(unittest.TestCase):
    """One project annotation, two omics: the same segment under both roots."""

    def test_both_omics_use_the_same_annotation_segment(self):
        anno = annotation.single("GRCh37")["id"]
        root = os.path.join("tmp", "result")
        rna = workspace.stage_dir("RNAseq", "D", "4_DEseq", root, anno)
        chip = workspace.stage_dir("ChIPseq", "D", "4_peak_result", root, anno)
        self.assertEqual(
            str(rna), os.path.join(root, "RNAseq", "D", "GRCh37-gencode", "4_deseq")
        )
        self.assertEqual(
            str(chip), os.path.join(root, "ChIPseq", "D", "GRCh37-gencode", "4_peaks")
        )

    def test_versions_of_one_dataset_do_not_collide(self):
        root = os.path.join("tmp", "result")
        old = workspace.result_dir("RNAseq", "D", root, annotation.single("hg19")["id"])
        new = workspace.result_dir("RNAseq", "D", root, annotation.single("GRCh38")["id"])
        self.assertNotEqual(str(old), str(new))

    def test_joint_result_is_also_versioned(self):
        anno = annotation.single("GRCh37")["id"]
        joint = workspace.integrated_dir("R", "C", os.path.join("tmp", "result"), anno)
        self.assertEqual(
            str(joint),
            os.path.join("tmp", "result", "Integrated", "R__C", "GRCh37-gencode"),
        )

    def test_blank_annotation_keeps_the_legacy_layout(self):
        root = os.path.join("tmp", "result")
        self.assertEqual(
            str(workspace.result_dir("RNAseq", "D", root)), os.path.join(root, "RNAseq", "D")
        )


class DeclaredInputsTest(unittest.TestCase):
    """ChIP Input libraries are named by ``input_sample`` and are not metadata rows.

    Stage 1 links reads by metadata row, so the Inputs have to be listed
    explicitly or peak calling opens a BAM that was never linked.
    """

    ROWS = [
        {"sample": "A1", "input_sample": "IN_A"},
        {"sample": "A2", "input_sample": "IN_A"},
        {"sample": "B1", "input_sample": "IN_B"},
    ]

    def test_inputs_are_listed_once_each(self):
        self.assertEqual(workspace.declared_inputs("ChIPseq", self.ROWS), ["IN_A", "IN_B"])

    def test_rna_rows_declare_no_inputs(self):
        self.assertEqual(workspace.declared_inputs("RNAseq", self.ROWS), [])

    def test_rows_without_an_input_are_skipped(self):
        self.assertEqual(
            workspace.declared_inputs("ChIPseq", [{"sample": "A1", "input_sample": ""}]),
            [],
        )

    def test_the_old_input_id_column_is_still_understood(self):
        rows = [{"sample": "A1", "input_id": "IN_A"}]
        self.assertEqual(workspace.declared_inputs("ChIPseq", rows), ["IN_A"])


class ChIPContractTest(unittest.TestCase):
    """The dataset contract the preparation page writes and the pipeline reads."""

    def test_chip_columns_are_the_role_aware_schema(self):
        self.assertEqual(
            workspace.columns_for("ChIPseq"),
            ("sample", "role", "sample_name", "condition", "input_sample"),
        )

    def test_comparisons_name_a_reference_and_a_treatment(self):
        self.assertEqual(workspace.COMPARISON_COLUMNS, ("reference", "treatment"))

    def test_comparisons_sit_next_to_the_metadata(self):
        root = os.path.join("tmp", "data")
        self.assertEqual(
            str(workspace.comparisons_path("ChIPseq", "demo", root)),
            os.path.join(root, "ChIPseq", "demo", "comparisons.csv"),
        )

    def test_rna_groups_are_a_closed_set_with_control_first(self):
        self.assertEqual(
            workspace.RNA_GROUPS, ("control", "treatA", "treatB", "treatC", "treatD")
        )

    def test_replicates_are_rep1_to_rep3(self):
        self.assertEqual(workspace.CHIP_REPLICATES, ("rep1", "rep2", "rep3"))


def load_script(name: str, filename: str):
    """Import a stage script, whose filename starts with a digit."""
    import importlib.util

    path = os.path.join(SCRIPT_ROOT, name, filename)
    spec = importlib.util.spec_from_file_location(filename[:-3].replace("_", ""), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Macs2ArgsTest(unittest.TestCase):
    """Synthetic references cannot fit MACS2's fragment model.

    They state ``--nomodel --extsize`` instead; the switches must reach every
    callpeak invocation, and an empty setting must add nothing.
    """

    def setUp(self):
        self.module = load_script("ChIPseq_pipline", "4_1_peakcalling.py")

    def test_extra_args_are_appended_in_order(self):
        self.assertEqual(
            self.module.extra_macs2_args({"chip": {"macs_extra_args": ["--nomodel", "--extsize", 150]}}),
            ["--nomodel", "--extsize", "150"],
        )

    def test_string_setting_is_split(self):
        self.assertEqual(
            self.module.extra_macs2_args({"chip": {"macs_extra_args": "--nomodel --extsize 150"}}),
            ["--nomodel", "--extsize", "150"],
        )

    def test_absent_setting_adds_nothing(self):
        self.assertEqual(self.module.extra_macs2_args({}), [])

    def test_genome_size_override_wins_over_species(self):
        info = {"species": "homo"}
        self.assertEqual(self.module.effective_genome({"chip": {"genome_size": 180000}}, info), "180000")
        self.assertEqual(self.module.effective_genome({"chip": {"genome_size": "auto"}}, info), "hs")


class PeakModeTest(unittest.TestCase):
    """A ``chip.peak_mode`` key matches a bare mark and a ``<mark>_<condition>`` group.

    The five-column schema has no mark column, so the key has to be found at the
    front of the group name; a dataset that writes ``H3K27ac_veh`` must not
    silently fall back to the default mode.
    """

    def setUp(self):
        self.module = load_script("ChIPseq_pipline", "4_1_peakcalling.py")

    def test_a_bare_key_matches_its_mark(self):
        config = {"chip": {"peak_mode": {"H3K27ac": "broad"}}}
        self.assertEqual(self.module.peak_mode_for(config, "H3K27ac"), "broad")

    def test_a_group_name_matches_the_key_before_its_underscore(self):
        config = {"chip": {"peak_mode": {"H3K27ac": "broad", "default": "narrow"}}}
        self.assertEqual(self.module.peak_mode_for(config, "H3K27ac_veh"), "broad")
        self.assertEqual(self.module.peak_mode_for(config, "H3K4me3_dox"), "narrow")

    def test_the_longest_prefix_wins_when_two_keys_overlap(self):
        config = {"chip": {"peak_mode": {"H3K27": "broad", "H3K27ac": "narrow"}}}
        self.assertEqual(self.module.peak_mode_for(config, "H3K27ac_veh"), "narrow")

    def test_an_unknown_mode_is_an_error(self):
        config = {"chip": {"peak_mode": {"H3K27ac": "wide"}}}
        with self.assertRaises(ValueError):
            self.module.peak_mode_for(config, "H3K27ac")


class PicardMetricsTest(unittest.TestCase):
    """A Picard metrics file has comment lines and blank lines before the header.

    Reading values before the ``LIBRARY`` header raised UnboundLocalError, so
    every ChIP QC summary aborted on the first library.
    """

    METRICS = """## htsjdk.samtools.metrics.StringHeader
# picard.sam.MarkDuplicates

## METRICS CLASS\tpicard.sam.DuplicationMetrics
LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE
Unknown Library\t200000\t0\t0\t0\t134010\t0\t0\t0.67\t

"""

    def setUp(self):
        self.module = load_script("ChIPseq_pipline", "4_4_qc_summary.py")
        self.tmp = tempfile.mkdtemp()
        self.path = os.path.join(self.tmp, "lib.metricsFile")
        with open(self.path, "w") as handle:
            handle.write(self.METRICS)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_blank_lines_before_the_header_are_skipped(self):
        self.assertAlmostEqual(self.module.duplication_rate(self.path), 0.67)

    def test_missing_file_is_not_an_error(self):
        self.assertIsNone(self.module.duplication_rate(os.path.join(self.tmp, "absent")))

    def test_header_without_values_is_not_an_error(self):
        with open(self.path, "w") as handle:
            handle.write(self.METRICS.split("Unknown Library")[0])
        self.assertIsNone(self.module.duplication_rate(self.path))



class MarkerArgumentTest(unittest.TestCase):
    """The dashboard passes the annotation as the fifth positional argument.

    The ``<assembly>-<source>`` segment added an argument to ``workspace.marker``;
    a panel process that kept the old module reported ``takes from 3 to 4
    positional arguments but 5 were given`` on every results page load.
    """

    def test_annotation_is_the_fifth_positional_argument(self):
        marker = workspace.marker("RNAseq", "D", "4_DEseq", "/tmp/results", "GRCh37-gencode")
        expected = workspace.stage_dir(
            "RNAseq", "D", "4_DEseq", "/tmp/results", "GRCh37-gencode"
        ) / workspace.FINISH_NAME
        self.assertEqual(marker, expected)
        self.assertIn("GRCh37-gencode", marker.parts)


class ResultViewSelectionTest(unittest.TestCase):
    """The results page reads a chosen ``<dataset>/<annotation>``, not the config's.

    When the configured annotation has no directory the page falls back to a
    version that does exist, so the path helpers have to take that version.
    """

    CONFIG = {"annotation": {"assembly": "GRCm38", "source": "gencode"}}

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import datasets

        self.datasets = datasets

    def test_annotation_argument_wins_over_config(self):
        path = self.datasets.stage(
            "RNAseq", "MOCK_HG19", "4_DEseq", self.CONFIG, "GRCh37-gencode"
        )
        self.assertEqual(path.parts[-3:], ("MOCK_HG19", "GRCh37-gencode", "4_deseq"))

    def test_config_annotation_is_used_without_an_argument(self):
        path = self.datasets.stage("RNAseq", "D", "4_DEseq", self.CONFIG)
        self.assertEqual(path.parts[-2], "GRCm38-gencode")

    def test_joint_stage_uses_the_picked_pair(self):
        path = self.datasets.joint_stage(
            self.CONFIG,
            "0_stats",
            rna="MOCK_HG19",
            chip="MOCK_HG19",
            annotation="GRCh37-gencode",
        )
        self.assertIn("MOCK_HG19__MOCK_HG19", path.parts)
        self.assertIn("GRCh37-gencode", path.parts)
        self.assertNotIn("GRCm38-gencode", path.parts)


class LiveReloadTest(unittest.TestCase):
    """A long-running panel kept pipeline modules that changed on disk.

    Streamlit's own watcher drops the module entry but leaves the package
    attribute behind, and ``workflow/scripts/common`` sits outside the app
    directory it watches, so the helper has to drop both packages itself.
    """

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import live_reload

        self.live_reload = live_reload
        self.live_reload._STAMP.clear()
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.source = os.path.join(self.tmp, "mod.py")
        with open(self.source, "w") as handle:
            handle.write("X = 1\n")

    def cached(self) -> None:
        sys.modules["fake_pkg"] = object()
        sys.modules["fake_pkg.mod"] = object()

    def test_first_call_only_records_the_sources(self):
        self.cached()
        self.assertEqual(self.live_reload.refresh(self.tmp, packages=("fake_pkg",)), [])
        self.assertIn("fake_pkg.mod", sys.modules)

    def test_unchanged_sources_keep_cached_modules(self):
        self.live_reload.refresh(self.tmp, packages=("fake_pkg",))
        self.cached()
        self.assertEqual(self.live_reload.refresh(self.tmp, packages=("fake_pkg",)), [])
        self.assertIn("fake_pkg.mod", sys.modules)

    def test_changed_source_drops_cached_modules(self):
        self.live_reload.refresh(self.tmp, packages=("fake_pkg",))
        self.cached()
        with open(self.source, "w") as handle:
            handle.write("X = 2\n")
        os.utime(self.source, (1, 1))
        self.assertEqual(
            self.live_reload.refresh(self.tmp, packages=("fake_pkg",)), [self.source]
        )
        self.assertNotIn("fake_pkg.mod", sys.modules)
        self.assertNotIn("fake_pkg", sys.modules)


class SignalTrackTest(unittest.TestCase):
    """The signal tab draws one row per track inside the contig's bounds.

    Two defects sat behind this: an opening window of 1 Mb fell off the mock's
    60 kb contigs (pyBigWig rejects the interval outright), and the same locus
    appeared once per sample in the peak list, so a highlight or a jump target
    was drawn several times over.
    """

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import tracks

        self.tracks = tracks

    def test_region_is_clipped_to_the_contig(self):
        locus = self.tracks.parse_region("chr1:50,000-200,000", {"chr1": 60_000})
        self.assertEqual((locus["start"], locus["end"]), (50_000, 60_000))

    def test_unknown_contig_names_the_ones_that_exist(self):
        with self.assertRaises(ValueError) as caught:
            self.tracks.parse_region("chrX:0-100", {"chr1": 60_000})
        self.assertIn("chr1", str(caught.exception))

    def test_overlapping_peaks_merge_into_one_band(self):
        peaks = [
            {"sample": "A", "chrom": "chr1", "start": 100, "end": 300, "score": 10.0},
            {"sample": "B", "chrom": "chr1", "start": 250, "end": 400, "score": 40.0},
            {"sample": "A", "chrom": "chr1", "start": 900, "end": 950, "score": 5.0},
        ]
        merged = self.tracks.peaks_in(peaks, "chr1", 0, 1_000)
        self.assertEqual(len(merged), 2)
        self.assertEqual((merged[0]["start"], merged[0]["end"]), (100, 400))
        self.assertEqual(merged[0]["score"], 40.0)

    def test_smoothing_does_not_pull_the_ends_down(self):
        # A flat track must stay flat at the edges; a zero-padded convolution
        # would leave the first and last bin at a fraction of the plateau.
        values = np.full(9, 7.0)
        smoothed = self.tracks.smooth(values, bin_bp=1.0, smooth_bp=3.0)
        self.assertTrue(np.allclose(smoothed, 7.0), smoothed)

    def test_rows_share_one_y_domain(self):
        frame = pd.DataFrame(
            {
                "track": ["A", "A", "B", "B"],
                "position": [0, 100, 0, 100],
                "signal": [1.0, 2.0, 10.0, 20.0],
            }
        )
        chart = self.tracks.build_chart(
            frame, chrom="chr1", start=0, end=100, order=["A", "B"], unified=True
        )
        spec = chart.to_dict()
        self.assertEqual(len(spec["vconcat"]), 2)
        domains = [
            layer["encoding"]["y"]["scale"]["domain"]
            for row in spec["vconcat"]
            for layer in row["layer"]
            if "y" in layer.get("encoding", {})
        ]
        self.assertTrue(domains)
        self.assertTrue(all(domain == [0, 20.0] for domain in domains), domains)

    def test_each_row_keeps_its_own_domain_when_unified_is_off(self):
        frame = pd.DataFrame(
            {
                "track": ["A", "A", "B", "B"],
                "position": [0, 100, 0, 100],
                "signal": [1.0, 2.0, 10.0, 20.0],
            }
        )
        chart = self.tracks.build_chart(
            frame, chrom="chr1", start=0, end=100, order=["A", "B"], unified=False
        )
        tops = [
            row["layer"][-1]["encoding"]["y"]["scale"]["domain"][1]
            for row in chart.to_dict()["vconcat"]
        ]
        self.assertEqual(tops, [2.0, 20.0])


class RunScopeTest(unittest.TestCase):
    """The "all stages" shortcut must point at a stage the Snakefile included."""

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import pipeline

        self.pipeline = pipeline
        self.config = {
            "annotation": {"assembly": "GRCm38", "source": "gencode"},
            "RNAseq": {"dataset": "demo"},
            "chip": {"dataset": "demo"},
            "integrated": {"mode": "descriptive"},
        }

    def test_single_condition_design_drops_diffbind(self):
        stages = self.pipeline.active_stages("ChIPseq", self.config)
        self.assertNotIn("5_diffbind", stages)
        self.assertEqual(stages[-1], "4_peak_result")
        self.assertIn("比较表", self.pipeline.stage_unavailable("ChIPseq", "5_diffbind", self.config))
        command = self.pipeline.scope_command("ChIPseq", self.config, 4)
        self.assertTrue(command[-1].endswith("4_peaks/finished.txt"), command[-1])
        # An absolute target made snakemake 7 abort with MissingRuleException
        # against rules that spell their outputs relative to the repo root.
        self.assertFalse(command[-1].startswith("/"), command[-1])

    def test_auto_mode_reads_the_comparison_file_next_to_the_metadata(self):
        # The page resolved the mode from metadata alone, so a dataset that
        # declares its pairs in comparisons.csv still looked descriptive and
        # 5_diffbind was greyed out while the Snakefile would have run it.
        config = dict(self.config, integrated={"mode": "auto"})
        samples = Path(write(ROLE_SCHEMA))
        pairs = Path(write(COMPARISONS))
        with mock.patch.object(
            self.pipeline.datasets.ws, "metadata_path", return_value=samples
        ), mock.patch.object(
            self.pipeline.datasets.ws, "comparisons_path", return_value=pairs
        ):
            self.assertEqual(self.pipeline.resolved_mode(config), "differential")
            self.assertIn("5_diffbind", self.pipeline.active_stages("ChIPseq", config))

    def test_auto_mode_without_a_comparison_file_stays_descriptive(self):
        config = dict(self.config, integrated={"mode": "auto"})
        with mock.patch.object(
            self.pipeline.datasets.ws, "metadata_path", return_value=Path(write(ROLE_SCHEMA))
        ), mock.patch.object(
            self.pipeline.datasets.ws,
            "comparisons_path",
            return_value=Path("/nonexistent/comparisons.csv"),
        ):
            self.assertEqual(self.pipeline.resolved_mode(config), "descriptive")
            self.assertNotIn("5_diffbind", self.pipeline.active_stages("ChIPseq", config))

    def test_two_condition_design_keeps_diffbind(self):
        config = dict(self.config, integrated={"mode": "differential"})
        self.assertIn("5_diffbind", self.pipeline.active_stages("ChIPseq", config))
        command = self.pipeline.scope_command("ChIPseq", config, 4)
        self.assertTrue(command[-1].endswith("5_diffbind/finished.txt"), command[-1])

    def test_switching_a_stage_off_pulls_the_target_back(self):
        config = dict(self.config)
        config["RNAseq_modules"] = {
            "1_prepare": True,
            "2_QC": True,
            "3_RC": True,
            "4_DEseq": False,
        }
        command = self.pipeline.scope_command("RNAseq", config, 4)
        self.assertTrue(command[-1].endswith("3_align_quant/finished.txt"), command[-1])

    def test_no_enabled_stage_is_a_clear_error(self):
        config = dict(self.config)
        config["RNAseq_modules"] = {
            name: False for name in ("1_prepare", "2_QC", "3_RC", "4_DEseq")
        }
        with self.assertRaises(ValueError) as caught:
            self.pipeline.scope_command("RNAseq", config, 4)
        self.assertIn("没有启用的阶段", str(caught.exception))

    def test_progress_counts_only_runnable_stages(self):
        stages = self.pipeline.active_stages("ChIPseq", self.config)
        _done, total = self.pipeline.progress("ChIPseq", self.config, only=stages)
        self.assertEqual(total, len(stages))


class LogCleanupTest(unittest.TestCase):
    """清理所有日志 deletes what the drawers list, and nothing else."""

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import env_paths, runtime

        self.runtime = runtime
        self.root = Path(tempfile.mkdtemp())
        self.logs = self.root / "logs"
        # runs sits under logs, so a naive walk counts every run log twice.
        self.runs = self.logs / "runs"
        (self.logs / "RNAseq").mkdir(parents=True)
        self.runs.mkdir(parents=True)
        self.live = self.runs / "20261028-070000_RNAseq-all.log"
        self.old = self.runs / "20261027-193358_download-RNAseq-demo.log"
        self.module = self.logs / "RNAseq" / "4_DEseq.log"
        self.server = self.logs / "streamlit.log"
        for path in (self.live, self.old, self.module, self.server):
            path.write_text("line\n" * 8)
        self.patches = [
            mock.patch.object(env_paths, "RUN_ROOT", self.runs),
            mock.patch.object(env_paths, "LOG_ROOT", self.logs),
            # A running job holds its log open, so the name must stay.
            mock.patch.object(
                runtime, "running_jobs", lambda: [{"scope": "RNAseq", "log": str(self.live)}]
            ),
        ]
        for patch in self.patches:
            patch.start()

    def tearDown(self):
        for patch in self.patches:
            patch.stop()
        shutil.rmtree(self.root, ignore_errors=True)

    def test_the_live_run_and_the_panel_own_log_are_untouched(self):
        count, size = self.runtime.log_footprint()
        self.assertEqual(count, 2)
        removed, freed = self.runtime.clear_logs()
        self.assertEqual(removed, 2)
        self.assertEqual(freed, size)
        self.assertTrue(self.live.is_file())
        self.assertFalse(self.old.exists())
        self.assertFalse(self.module.exists())
        self.assertTrue(self.server.is_file())


class PrepareTableTest(unittest.TestCase):
    """The naming table holds one row per library and numbers the Inputs."""

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import prepare

        self.prepare = prepare
        self.rows = [
            {
                "sample": "IN1",
                "role": "Input",
                "sample_name": "IN_veh",
                "condition": "rep1",
                "input_sample": "",
            },
            {
                "sample": "A1",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "IN1",
            },
            {
                "sample": "A2",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep2",
                "input_sample": "IN1",
            },
        ]

    def test_an_input_row_is_kept_beside_the_ip_rows(self):
        normalized = self.prepare.normalize_rows("ChIPseq", self.rows)
        self.assertEqual([row["role"] for row in normalized], ["Input", "IP", "IP"])
        self.assertEqual(
            [row["sample_name"] for row in normalized],
            ["IN_veh", "H3K27ac_veh", "H3K27ac_veh"],
        )

    def test_a_row_cited_as_input_gets_the_input_role(self):
        rows = [
            {"sample": "IN1", "sample_name": "", "condition": "", "input_sample": ""},
            {
                "sample": "A1",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "IN1",
            },
        ]
        normalized = self.prepare.normalize_rows("ChIPseq", rows)
        self.assertEqual(normalized[0]["role"], "Input")
        self.assertEqual(normalized[0]["condition"], "rep1")

    def test_two_libraries_of_one_input_get_their_own_replicate(self):
        rows = [
            {
                "sample": "IN1",
                "role": "Input",
                "sample_name": "IN_veh",
                "condition": "",
                "input_sample": "",
            },
            {
                "sample": "IN2",
                "role": "Input",
                "sample_name": "IN_veh",
                "condition": "",
                "input_sample": "",
            },
        ]
        normalized = self.prepare.normalize_rows("ChIPseq", rows)
        self.assertEqual([row["condition"] for row in normalized], ["rep1", "rep2"])

    def test_a_blank_replicate_does_not_reuse_the_explicit_one(self):
        rows = [
            {
                "sample": "IN1",
                "role": "Input",
                "sample_name": "IN_veh",
                "condition": "rep1",
                "input_sample": "",
            },
            {
                "sample": "IN2",
                "role": "Input",
                "sample_name": "IN_veh",
                "condition": "",
                "input_sample": "",
            },
        ]
        normalized = self.prepare.normalize_rows("ChIPseq", rows)
        self.assertEqual([row["condition"] for row in normalized], ["rep1", "rep2"])

    def test_a_missing_input_row_is_added_where_step_one_listed_it(self):
        rows = [
            {
                "sample": "A1",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "IN1",
            }
        ]
        added = self.prepare.add_input_rows(rows, ["IN1", "A1"], {"IN1": "IN_veh"})
        self.assertEqual([row["sample"] for row in added], ["IN1", "A1"])
        self.assertEqual(added[0]["role"], "Input")
        self.assertEqual(added[0]["sample_name"], "IN_veh")
        self.assertEqual(added[0]["condition"], "rep1")

    def test_a_hand_typed_input_id_with_no_library_adds_no_row(self):
        rows = [{"sample": "A1", "input_sample": "IN_missing"}]
        self.assertEqual(
            [row["sample"] for row in self.prepare.add_input_rows(rows, ["A1"])],
            ["A1"],
        )

    def test_an_input_row_only_needs_a_group_name(self):
        rows = [
            {
                "sample": "IN1",
                "role": "Input",
                "sample_name": "IN_veh",
                "condition": "rep1",
                "input_sample": "",
            }
        ]
        self.assertEqual(self.prepare.missing_fields("ChIPseq", rows), [])

    def test_an_ip_row_still_needs_its_input(self):
        rows = [
            {
                "sample": "A1",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "",
            }
        ]
        hints = self.prepare.missing_fields("ChIPseq", rows)
        self.assertTrue(any("input_sample" in hint for hint in hints), hints)

    def test_an_input_group_is_not_offered_as_a_comparison(self):
        self.assertEqual(
            self.prepare.group_summary("ChIPseq", self.rows), [("H3K27ac_veh", 2)]
        )

    def test_two_libraries_of_one_input_cannot_share_a_replicate(self):
        rows = self.rows[:1] + [dict(self.rows[0], sample="IN2")]
        notes = self.prepare.group_notes("ChIPseq", rows)
        self.assertTrue(any("IN_veh / rep1" in note for note in notes), notes)

    def test_rna_rows_are_returned_untouched(self):
        rows = [{"sample": "S1", "sample_name": "S1_1", "group": "control"}]
        self.assertIs(self.prepare.normalize_rows("RNAseq", rows), rows)

    def test_an_input_sample_pointing_at_the_sample_itself_is_refused(self):
        rows = [
            {
                "sample": "A1",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "A1",
            }
        ]
        errors, _warnings = self.prepare.datasets.validate_metadata("ChIPseq", rows)
        self.assertTrue(any("input_sample" in error for error in errors), errors)

    def test_chip_candidates_are_shared_by_sample_and_input(self):
        rows = [{"sample": "A1", "input_sample": "IN1"}]
        with mock.patch.object(
            self.prepare.datasets, "fastq_table", return_value={"A1": {}, "IN1": {}}
        ):
            options = self.prepare.sample_options("ChIPseq", "demo", ["A2", "IN1"], rows)
        self.assertEqual(options, ["A1", "A2", "IN1"])

    def test_an_old_row_is_migrated_onto_the_role_columns(self):
        row = {
            "sample": "A1",
            "mark": "H3K27ac",
            "condition": "Veh",
            "replicate": "1",
            "input_id": "IN1",
        }
        self.assertEqual(
            self.prepare.chip_row(row),
            {
                "sample": "A1",
                "role": "IP",
                "sample_name": "H3K27ac_Veh",
                "condition": "rep1",
                "input_sample": "IN1",
            },
        )

    def test_a_single_condition_old_row_keeps_the_bare_mark(self):
        row = {
            "sample": "A1",
            "mark": "H3K27ac",
            "condition": "single",
            "replicate": "2",
            "input_id": "IN1",
        }
        migrated = self.prepare.chip_row(row)
        self.assertEqual(migrated["sample_name"], "H3K27ac")
        self.assertEqual(migrated["condition"], "rep2")

    def test_rna_rows_are_returned_untouched_by_the_migrator(self):
        rows = [{"sample": "S1", "group": "control"}]
        self.assertIs(self.prepare.normalize_rows("RNAseq", rows), rows)

    def test_a_row_without_a_group_name_is_counted_as_unnamed(self):
        rows = [
            {"sample": "A1", "sample_name": ""},
            {"sample": "A2", "sample_name": "H3K27ac_veh"},
            {"sample": "A3", "sample_name": "H3K27ac_veh"},
        ]
        self.assertEqual(
            self.prepare.group_summary("ChIPseq", rows),
            [("H3K27ac_veh", 2), ("未填", 1)],
        )
        # The unnamed row is already named by missing_fields, and the group it
        # would have joined is complete; nothing to add here.
        self.assertEqual(self.prepare.group_notes("ChIPseq", rows), [])

    def test_a_group_with_one_rep_is_reported_while_typing(self):
        rows = [
            {"sample": "A1", "sample_name": "H3K27ac_veh", "condition": "rep1"},
            {"sample": "B1", "sample_name": "H3K27ac_tgfb", "condition": "rep1"},
        ]
        notes = self.prepare.group_notes("ChIPseq", rows)
        self.assertTrue(any("H3K27ac_veh 只有 1 例" in note for note in notes), notes)
        self.assertTrue(any("H3K27ac_tgfb 只有 1 例" in note for note in notes), notes)

    def test_two_samples_claiming_the_same_rep_are_reported(self):
        rows = [
            {"sample": "A1", "sample_name": "H3K27ac_veh", "condition": "rep1"},
            {"sample": "A2", "sample_name": "H3K27ac_veh", "condition": "rep1"},
        ]
        notes = self.prepare.group_notes("ChIPseq", rows)
        self.assertTrue(
            any("H3K27ac_veh / rep1" in note and "同时占用" in note for note in notes),
            notes,
        )

    def test_a_condition_that_is_not_a_rep_is_called_out(self):
        rows = [{"sample": "A1", "sample_name": "H3K27ac_veh", "condition": "Veh"}]
        notes = self.prepare.group_notes("ChIPseq", rows)
        self.assertTrue(any("rep1 / rep2 / rep3" in note for note in notes), notes)

    def test_names_differing_only_by_case_are_called_out(self):
        rows = [
            {"sample": "S1", "sample_name": "S1", "group": "control"},
            {"sample": "S2", "sample_name": "S2", "group": "Control"},
        ]
        notes = self.prepare.group_notes("RNAseq", rows)
        self.assertTrue(any("大小写" in note for note in notes), notes)

    def test_chip_group_names_differing_only_by_case_are_called_out(self):
        rows = [
            {"sample": "A1", "sample_name": "H3K27ac_veh", "condition": "rep1"},
            {"sample": "A2", "sample_name": "H3K27ac_VEH", "condition": "rep1"},
        ]
        notes = self.prepare.group_notes("ChIPseq", rows)
        self.assertTrue(any("大小写" in note for note in notes), notes)

    def test_rows_with_every_cell_empty_are_dropped(self):
        rows = [
            {"sample": "A1", "mark": "H3K27ac"},
            {"sample": "", "mark": ""},
            {"sample": None, "mark": None},
        ]
        self.assertEqual(
            self.prepare.drop_blank_rows(rows), [{"sample": "A1", "mark": "H3K27ac"}]
        )

    def test_a_saved_input_row_is_drafted_again_when_it_is_missing(self):
        stored = [row for row in self.rows if row["sample"] != "IN1"]
        with mock.patch.object(self.prepare.datasets, "metadata_rows", return_value=stored):
            drafted = self.prepare.draft_rows("ChIPseq", "demo", ["IN1", "A1", "A2"])
        # The IP rows cite IN1, so the table cannot be saved without its row.
        self.assertEqual([row["sample"] for row in drafted], ["IN1", "A1", "A2"])
        self.assertEqual(drafted[0]["role"], "Input")

    def test_deleted_rows_stay_deleted_for_the_same_sample_set(self):
        stored = [row for row in self.rows if row["sample"] != "A2"]
        saved = {"table": {"chosen": ["IN1", "A1", "A2"], "rows": ["IN1", "A1"]}}
        with mock.patch.object(
            self.prepare.datasets, "metadata_rows", return_value=stored
        ), mock.patch.object(self.prepare.datasets, "descriptors", return_value=saved):
            drafted = self.prepare.draft_rows("ChIPseq", "demo", ["IN1", "A1", "A2"])
        self.assertEqual([row["sample"] for row in drafted], ["IN1", "A1"])

    def test_saving_records_which_rows_and_samples_the_table_had(self):
        written = {}
        fastq = {name: {"r1": "x.fastq.gz"} for name in ("IN1", "A1", "A2")}
        with mock.patch.object(
            self.prepare.datasets, "validate_metadata", return_value=([], [])
        ), mock.patch.object(
            self.prepare.datasets, "fastq_table", return_value=fastq
        ), mock.patch.object(
            self.prepare.datasets, "save_metadata", return_value="metadata.csv"
        ), mock.patch.object(
            self.prepare.datasets, "descriptors", return_value={}
        ), mock.patch.object(
            self.prepare.datasets,
            "save_descriptor",
            side_effect=lambda _omics, _name, payload: written.update(payload),
        ):
            errors, _warnings, _path = self.prepare.complete(
                "ChIPseq", "demo", self.rows, "geo", chosen=["IN1", "A1", "A2"]
            )
        self.assertEqual(errors, [])
        self.assertEqual(written["table"]["rows"], ["IN1", "A1", "A2"])
        self.assertEqual(written["table"]["chosen"], ["IN1", "A1", "A2"])

    def test_a_selected_run_still_to_download_is_saved_with_a_warning(self):
        fastq = {"A1": {"r1": "x.fastq.gz"}}
        with mock.patch.object(
            self.prepare.datasets, "validate_metadata", return_value=([], [])
        ), mock.patch.object(
            self.prepare.datasets, "fastq_table", return_value=fastq
        ), mock.patch.object(
            self.prepare.datasets, "save_metadata", return_value="metadata.csv"
        ), mock.patch.object(
            self.prepare.datasets, "srr_list", return_value=["IN1", "A2"]
        ), mock.patch.object(
            self.prepare.datasets, "descriptors", return_value={}
        ), mock.patch.object(
            self.prepare.datasets, "save_descriptor"
        ):
            errors, warnings, path = self.prepare.complete(
                "ChIPseq", "demo", self.rows, "geo"
            )
        # Step 1 lists these runs, so the gap is a download to finish rather
        # than a name to correct; the table is still worth saving.
        self.assertEqual(errors, [])
        self.assertIsNotNone(path)
        self.assertTrue(any("还没下载 FASTQ" in item for item in warnings), warnings)
        self.assertFalse(any("找不到 FASTQ" in item for item in warnings), warnings)

    def test_a_row_outside_the_step_one_selection_is_reported_as_missing(self):
        fastq = {"A1": {"r1": "x.fastq.gz"}}
        with mock.patch.object(
            self.prepare.datasets, "validate_metadata", return_value=([], [])
        ), mock.patch.object(
            self.prepare.datasets, "fastq_table", return_value=fastq
        ), mock.patch.object(
            self.prepare.datasets, "save_metadata", return_value="metadata.csv"
        ), mock.patch.object(
            self.prepare.datasets, "srr_list", return_value=[]
        ), mock.patch.object(
            self.prepare.datasets, "descriptors", return_value={}
        ), mock.patch.object(
            self.prepare.datasets, "save_descriptor"
        ):
            errors, warnings, _path = self.prepare.complete(
                "ChIPseq", "demo", self.rows, "geo"
            )
        self.assertEqual(errors, [])
        self.assertTrue(any("找不到 FASTQ" in item for item in warnings), warnings)

    def test_a_new_sample_selection_drafts_every_row_again(self):
        saved = {"table": {"chosen": ["IN1", "A1", "A2"], "rows": ["A1", "A2"]}}
        with mock.patch.object(
            self.prepare.datasets, "metadata_rows", return_value=[]
        ), mock.patch.object(self.prepare.datasets, "descriptors", return_value=saved):
            drafted = self.prepare.draft_rows(
                "ChIPseq", "demo", ["IN1", "A1", "A2", "A3"]
            )
        self.assertEqual(
            [row["sample"] for row in drafted], ["IN1", "A1", "A2", "A3"]
        )

    def test_a_draft_survives_a_new_session(self):
        with tempfile.TemporaryDirectory() as root:
            with mock.patch.object(
                self.prepare.datasets, "directory", return_value=Path(root) / "demo"
            ):
                self.prepare.write_draft(
                    "ChIPseq",
                    "demo",
                    {"source": "geo", "chosen": ["A1", "A2"], "rows": self.rows[1:]},
                )
                draft = self.prepare.read_draft("ChIPseq", "demo")
        self.assertEqual(draft["chosen"], ["A1", "A2"])
        self.assertEqual(draft["rows"], self.rows[1:])
        self.assertIn("updated_at", draft)

    def test_discarding_a_selection_removes_the_draft(self):
        with tempfile.TemporaryDirectory() as root:
            with mock.patch.object(
                self.prepare.datasets, "directory", return_value=Path(root) / "demo"
            ):
                self.prepare.write_draft("ChIPseq", "demo", {"chosen": ["A1"]})
                self.prepare.clear_draft("ChIPseq", "demo")
                draft = self.prepare.read_draft("ChIPseq", "demo")
        self.assertEqual(draft, {})

    def test_a_draft_with_no_file_reads_as_empty(self):
        with tempfile.TemporaryDirectory() as root:
            with mock.patch.object(
                self.prepare.datasets, "directory", return_value=Path(root) / "missing"
            ):
                self.assertEqual(self.prepare.read_draft("ChIPseq", "missing"), {})


class PreparationEditorTest(unittest.TestCase):
    """The two tables' column definitions, read off the page module.

    The page renders itself on import, so the warnings it prints in a bare
    interpreter are muted; what matters here is what the editor is told to
    offer.
    """

    @classmethod
    def setUpClass(cls):
        import importlib.util
        import logging

        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        path = os.path.join(app_root, "views", "3_组学准备.py")
        spec = importlib.util.spec_from_file_location("prepare_page", path)
        module = importlib.util.module_from_spec(spec)
        logging.disable(logging.CRITICAL)
        try:
            spec.loader.exec_module(module)
        finally:
            logging.disable(logging.NOTSET)
        cls.page = module

    def test_chip_editor_has_the_role_aware_columns(self):
        columns = self.page.editor_columns("ChIPseq", ["A1", "IN1"])
        self.assertEqual(
            list(columns),
            ["sample", "role", "sample_name", "condition", "input_sample"],
        )

    def test_the_role_column_offers_only_ip_and_input(self):
        columns = self.page.editor_columns("ChIPseq", ["A1", "IN1"])
        self.assertEqual(columns["role"]["type_config"]["options"], ["IP", "Input"])

    def test_chip_condition_offers_only_rep1_to_rep3(self):
        columns = self.page.editor_columns("ChIPseq", ["A1", "IN1"])
        self.assertEqual(
            columns["condition"]["type_config"]["options"], ["", "rep1", "rep2", "rep3"]
        )

    def test_the_input_dropdown_shares_the_sample_candidates(self):
        columns = self.page.editor_columns("ChIPseq", ["A1", "IN1"])
        self.assertEqual(
            columns["input_sample"]["type_config"]["options"],
            columns["sample"]["type_config"]["options"],
        )

    def test_rna_group_is_a_fixed_dropdown(self):
        columns = self.page.editor_columns("RNAseq", ["A1"])
        self.assertEqual(
            columns["group"]["type_config"]["options"],
            ["control", "treatA", "treatB", "treatC", "treatD"],
        )

    def test_rna_sample_stays_a_text_column(self):
        # The page disables it in the editor call; the config itself just keeps
        # it a plain text box so the disabled rule has something to apply to.
        columns = self.page.editor_columns("RNAseq", ["A1"])
        self.assertEqual(columns["sample"]["type_config"]["type"], "text")


class MetadataValidationTest(unittest.TestCase):
    """The two tables' save-time rules, which the page mirrors while typing."""

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import datasets

        self.datasets = datasets

    @staticmethod
    def rna(*groups):
        return [
            {"sample": f"S{index}", "sample_name": f"S{index}_name", "group": group}
            for index, group in enumerate(groups)
        ]

    @staticmethod
    def chip(*rows):
        return [
            {
                "sample": f"S{index}",
                "sample_name": name,
                "condition": condition,
                "input_sample": "IN1",
            }
            for index, (name, condition) in enumerate(rows)
        ]

    def test_a_complete_rna_design_passes(self):
        errors, warnings = self.datasets.validate_metadata(
            "RNAseq", self.rna("control", "control", "treatA", "treatA")
        )
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])

    def test_rna_needs_a_control_reference(self):
        errors, _warnings = self.datasets.validate_metadata(
            "RNAseq", self.rna("treatA", "treatA")
        )
        self.assertTrue(any("control" in error for error in errors), errors)

    def test_rna_needs_a_treat_group(self):
        errors, _warnings = self.datasets.validate_metadata(
            "RNAseq", self.rna("control", "control")
        )
        self.assertTrue(any("treat" in error for error in errors), errors)

    def test_a_group_outside_the_five_names_is_refused(self):
        errors, _warnings = self.datasets.validate_metadata(
            "RNAseq", self.rna("control", "control", "Veh", "Veh")
        )
        self.assertTrue(any("Veh" in error for error in errors), errors)

    def test_a_treat_with_one_replicate_is_an_error(self):
        errors, _warnings = self.datasets.validate_metadata(
            "RNAseq", self.rna("control", "control", "treatA")
        )
        self.assertTrue(any("treatA" in error for error in errors), errors)

    def test_a_complete_chip_design_passes(self):
        errors, warnings = self.datasets.validate_metadata(
            "ChIPseq", self.chip(("H3K27ac_veh", "rep1"), ("H3K27ac_veh", "rep2"))
        )
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])

    def test_a_repeated_rep_is_an_error(self):
        errors, _warnings = self.datasets.validate_metadata(
            "ChIPseq", self.chip(("H3K27ac_veh", "rep1"), ("H3K27ac_veh", "rep1"))
        )
        self.assertTrue(any("duplicate" in error for error in errors), errors)

    def test_a_missing_input_sample_is_an_error(self):
        rows = [
            {
                "sample": "S0",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "",
            }
        ]
        errors, _warnings = self.datasets.validate_metadata("ChIPseq", rows)
        self.assertTrue(any("input_sample" in error for error in errors), errors)

    def test_a_lonely_group_is_a_warning_not_an_error(self):
        errors, warnings = self.datasets.validate_metadata(
            "ChIPseq", self.chip(("H3K27ac_veh", "rep1"))
        )
        self.assertEqual(errors, [])
        self.assertTrue(any("H3K27ac_veh" in warning for warning in warnings), warnings)

    def test_a_role_aware_design_passes(self):
        rows = [
            {
                "sample": "IN1",
                "role": "Input",
                "sample_name": "input_veh",
                "condition": "rep1",
                "input_sample": "",
            },
            {
                "sample": "A1",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "IN1",
            },
            {
                "sample": "A2",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep2",
                "input_sample": "IN1",
            },
        ]
        errors, warnings = self.datasets.validate_metadata("ChIPseq", rows)
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])

    def test_an_ip_citing_a_row_without_the_input_role_is_an_error(self):
        rows = [
            {
                "sample": "A1",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "IN1",
            }
        ]
        errors, _warnings = self.datasets.validate_metadata("ChIPseq", rows)
        self.assertTrue(
            any("not a row with role Input" in error for error in errors), errors
        )

    def test_an_input_row_cannot_name_an_input_sample(self):
        rows = [
            {
                "sample": "IN1",
                "role": "Input",
                "sample_name": "input_veh",
                "condition": "rep1",
                "input_sample": "IN2",
            }
        ]
        errors, _warnings = self.datasets.validate_metadata("ChIPseq", rows)
        self.assertTrue(
            any("cannot name an input_sample" in error for error in errors), errors
        )

    def test_an_input_nobody_cites_is_a_warning(self):
        rows = [
            {
                "sample": "IN1",
                "role": "Input",
                "sample_name": "input_veh",
                "condition": "rep1",
                "input_sample": "",
            },
            {
                "sample": "IN2",
                "role": "Input",
                "sample_name": "input_veh",
                "condition": "rep2",
                "input_sample": "",
            },
            {
                "sample": "A1",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "IN1",
            },
            {
                "sample": "A2",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep2",
                "input_sample": "IN1",
            },
        ]
        _errors, warnings = self.datasets.validate_metadata("ChIPseq", rows)
        self.assertTrue(any("IN2" in warning for warning in warnings), warnings)

    def test_an_unknown_role_is_an_error(self):
        rows = [
            {
                "sample": "A1",
                "role": "chip",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "IN1",
            }
        ]
        errors, _warnings = self.datasets.validate_metadata("ChIPseq", rows)
        self.assertTrue(any("role must be IP or Input" in error for error in errors), errors)

    def test_two_inputs_without_a_rep_are_a_warning(self):
        rows = [
            {
                "sample": "IN1",
                "role": "Input",
                "sample_name": "input_veh",
                "condition": "",
                "input_sample": "",
            },
            {
                "sample": "IN2",
                "role": "Input",
                "sample_name": "input_veh",
                "condition": "",
                "input_sample": "",
            },
            {
                "sample": "A1",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep1",
                "input_sample": "IN1",
            },
            {
                "sample": "A2",
                "role": "IP",
                "sample_name": "H3K27ac_veh",
                "condition": "rep2",
                "input_sample": "IN2",
            },
        ]
        _errors, warnings = self.datasets.validate_metadata("ChIPseq", rows)
        self.assertTrue(any("没有 condition" in warning for warning in warnings), warnings)

    def test_comparisons_round_trip_through_the_dataset_directory(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "comparisons.csv"
            with mock.patch.object(self.datasets.ws, "comparisons_path", return_value=path):
                self.datasets.save_comparisons(
                    "ChIPseq",
                    "demo",
                    [{"reference": "H3K27ac_veh", "treatment": "H3K27ac_tgfb"}],
                )
                self.assertEqual(
                    self.datasets.comparison_rows("ChIPseq", "demo"),
                    [{"reference": "H3K27ac_veh", "treatment": "H3K27ac_tgfb"}],
                )


class AnnotationSelectionTest(unittest.TestCase):
    """The build page describes the assembly picked in the selectbox."""

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import runtime

        self.runtime = runtime

    def test_switching_to_human_shows_the_human_txdb(self):
        config = {"annotation": {"assembly": "GRCm38", "source": "gencode"}, "threads": 8}
        merged = self.runtime.deep_merge(
            config, {"annotation": {"assembly": "GRCh38", "source": "gencode"}}
        )
        info = annotation.resolve(merged)
        self.assertEqual(info["id"], "GRCh38-gencode")
        self.assertEqual(info["txdb"], "TxDb.Hsapiens.UCSC.hg38.knownGene")
        self.assertEqual(info["orgdb"], "org.Hs.eg.db")
        self.assertEqual(info["taxon"], "human")

    def test_merging_the_selection_keeps_the_other_sections(self):
        config = {"annotation": {"assembly": "GRCm38", "source": "gencode"}, "threads": 8}
        merged = self.runtime.deep_merge(
            config, {"annotation": {"assembly": "GRCh38"}}
        )
        self.assertEqual(merged["threads"], 8)
        self.assertEqual(merged["annotation"]["source"], "gencode")

    def test_mouse_keeps_the_mouse_txdb(self):
        info = annotation.resolve({"annotation": {"assembly": "GRCm38", "source": "gencode"}})
        self.assertEqual(info["txdb"], "TxDb.Mmusculus.UCSC.mm10.knownGene")
        self.assertEqual(info["orgdb"], "org.Mm.eg.db")


class AnnotationFileLayoutTest(unittest.TestCase):
    """The file list points at what the build script actually writes."""

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import annotation as ui_annotation

        self.ui = ui_annotation

    def test_rna_index_sits_in_the_assembly_named_directory(self):
        items = self.ui.expected_files("rna", "GRCm38", "gencode")
        paths = [path for _label, path in items if path.suffix == ".ht2"]
        self.assertEqual(len(paths), 8)
        self.assertEqual({path.parent.name for path in paths}, {"grcm38"})
        self.assertEqual(
            {path.name for path in paths},
            {f"genome.{number}.ht2" for number in range(1, 9)},
        )

    def test_splice_files_are_named_after_the_index_directory(self):
        names = [
            path.name
            for _label, path in self.ui.expected_files("rna", "GRCm38", "gencode")
        ]
        self.assertIn("splicesites_grcm38.txt", names)
        self.assertIn("exons_grcm38.txt", names)

    def test_chip_index_keeps_the_assembly_prefix(self):
        names = [
            path.name
            for _label, path in self.ui.expected_files("chip", "GRCm38", "gencode")
        ]
        self.assertIn("grcm38.1.bt2", names)
        self.assertIn("grcm38.rev.2.bt2", names)

    def test_index_prefix_points_at_the_aligner_stem(self):
        rna = self.ui.index_prefix("rna", "GRCm38")
        self.assertEqual(rna.name, "genome")
        self.assertEqual(rna.parent.name, "grcm38")
        self.assertEqual(self.ui.index_prefix("chip", "GRCm38").name, "grcm38")


class AnnotationProgressTest(unittest.TestCase):
    """The build bar follows the JSON the script writes, weights included."""

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import annotation as ui_annotation

        self.ui = ui_annotation

    def test_missing_weights_split_the_bar_evenly(self):
        view = self.ui.progress_view(
            {
                "steps": ["prepare", "genome", "gtf", "index"],
                "step_index": 2,
                "step_progress": 0.5,
            }
        )
        self.assertAlmostEqual(view["overall"], 0.625)

    def test_weights_put_the_bar_where_the_time_goes(self):
        view = self.ui.progress_view(
            {
                "steps": ["prepare", "index"],
                "weights": [1, 9],
                "step_index": 1,
                "step_progress": 0.5,
            }
        )
        self.assertAlmostEqual(view["overall"], 0.55)
        self.assertEqual(view["rows"][0][2], "完成")
        self.assertEqual(view["rows"][1][2], "进行中")

    def test_a_finished_build_is_a_full_bar(self):
        view = self.ui.progress_view(
            {
                "steps": ["prepare", "index"],
                "weights": [1, 9],
                "step_index": 0,
                "step_progress": 0.2,
                "done": True,
                "message": "索引已存在，未重建",
            }
        )
        self.assertEqual(view["overall"], 1.0)
        self.assertEqual([row[2] for row in view["rows"]], ["完成", "完成"])
        self.assertIn("索引已存在", view["text"])

    def test_every_step_the_script_writes_has_a_label(self):
        script = Path(__file__).resolve().parents[1] / "scripts" / "build-annotation.sh"
        found = re.findall(r"STEPS=\(([^)]*)\)", script.read_text(encoding="utf-8"))
        self.assertEqual(len(found), 2, "rna and chip each declare their steps")
        for body in found:
            for step in body.split():
                self.assertIn(step, self.ui.STEP_LABEL)

    def test_starting_a_build_drops_the_previous_run(self):
        with tempfile.TemporaryDirectory() as root:
            target = Path(root) / "annotation-rna.json"
            target.write_text('{"done": true}')
            with mock.patch.object(self.ui, "progress_file", return_value=target):
                returned = self.ui.begin_progress("rna")
            self.assertEqual(returned, target)
            self.assertFalse(target.exists())

    def test_each_index_build_owns_a_job_slot(self):
        # One slot for both kinds serialised HISAT2 behind bowtie2 for no
        # reason: the two builds share no file.
        self.assertEqual(
            {self.ui.job_scope(kind) for kind in ("rna", "chip")},
            {"annotation-rna", "annotation-chip"},
        )
        self.assertEqual(self.ui.progress_file("chip").stem, self.ui.job_scope("chip"))

    def test_a_build_without_a_progress_file_reads_as_empty(self):
        with tempfile.TemporaryDirectory() as root:
            with mock.patch.object(
                self.ui, "progress_file", return_value=Path(root) / "missing.json"
            ):
                self.assertEqual(self.ui.read_progress("chip"), {})
                self.assertEqual(self.ui.progress_view({})["rows"], [])



class DownloadProgressTest(unittest.TestCase):
    """The download bar follows the JSON fetch_sra writes, bytes included."""

    def setUp(self):
        app_root = os.path.join(os.path.dirname(__file__), "..", "streamlit_app")
        if app_root not in sys.path:
            sys.path.insert(0, app_root)
        from utils import prepare

        self.prepare = prepare

    def test_missing_weights_split_the_bar_evenly(self):
        view = self.prepare.download_view(
            {
                "steps": ["download", "convert", "compress"],
                "step_index": 1,
                "step_progress": 0.5,
            }
        )
        self.assertAlmostEqual(view["overall"], 0.5)
        self.assertEqual([row[2] for row in view["rows"]], ["完成", "进行中", "未开始"])

    def test_byte_progress_lands_inside_the_download_step(self):
        view = self.prepare.download_view(
            {
                "steps": ["download", "convert", "compress"],
                "weights": [6, 3, 1],
                "step_index": 0,
                "step_progress": 0.5,
                "detail": "1.0 GB / 2.0 GB",
            }
        )
        self.assertAlmostEqual(view["overall"], 0.3)
        self.assertIn("30%", view["text"])
        self.assertIn("1.0 GB / 2.0 GB", view["text"])

    def test_a_download_that_ended_with_failures_is_not_a_full_bar(self):
        view = self.prepare.download_view(
            {
                "steps": ["download", "convert", "compress"],
                "weights": [6, 3, 1],
                "step_index": 2,
                "step_progress": 1.0,
                "done": True,
                "finished": 4,
                "total": 6,
                "error": "SRR9: prefetch failed",
            }
        )
        self.assertAlmostEqual(view["overall"], 4 / 6)
        self.assertEqual(view["error"], "SRR9: prefetch failed")
        self.assertEqual([row[2] for row in view["rows"]], ["结束", "结束", "结束"])
        self.assertIn("有样本失败", view["text"])

    def test_a_finished_download_is_a_full_bar(self):
        view = self.prepare.download_view(
            {
                "steps": ["download", "convert", "compress"],
                "weights": [6, 3, 1],
                "step_index": 1,
                "step_progress": 0.5,
                "done": True,
                "message": "GSE1 · 2 个样本",
            }
        )
        self.assertEqual(view["overall"], 1.0)
        self.assertEqual([row[2] for row in view["rows"]], ["完成", "完成", "完成"])
        self.assertIn("GSE1", view["text"])

    def test_each_run_row_carries_its_stage_and_bytes(self):
        view = self.prepare.download_view(
            {
                "steps": ["download"],
                "runs": [
                    {"run": "SRR1", "stage": "compress", "bytes": 1500, "error": ""},
                    {"run": "SRR2", "stage": "fail", "bytes": 0, "error": "prefetch failed"},
                ],
            }
        )
        self.assertEqual([row["状态"] for row in view["runs"]], ["压缩中", "失败"])
        self.assertEqual(view["runs"][0]["已下载"], "1.5 KB")
        self.assertEqual(view["runs"][1]["错误"], "prefetch failed")

    def test_every_stage_the_script_writes_has_a_label(self):
        script = (
            Path(__file__).resolve().parents[1]
            / "workflow"
            / "scripts"
            / "common"
            / "fetch_sra.py"
        )
        source = script.read_text(encoding="utf-8")
        steps = re.search(r"STEPS = \(([^)]*)\)", source)
        self.assertIsNotNone(steps, "the script declares its steps")
        for step in re.findall(r'"([a-z]+)"', steps.group(1)):
            self.assertIn(step, self.prepare.DOWNLOAD_LABEL)
        block = re.search(r"STAGE_INDEX = \{(.*?)\}", source, re.S)
        for stage in re.findall(r'"([a-z]+)":', block.group(1)):
            self.assertIn(stage, self.prepare.RUN_STAGE_LABEL)

    def test_the_progress_file_is_named_after_the_job_slot(self):
        self.assertEqual(self.prepare.download_scope("ChIPseq"), "download-ChIPseq")
        self.assertEqual(
            self.prepare.download_progress_file("ChIPseq").stem,
            self.prepare.download_scope("ChIPseq"),
        )

    def test_starting_a_download_drops_the_previous_run(self):
        with tempfile.TemporaryDirectory() as root:
            target = Path(root) / "download-RNAseq.json"
            target.write_text('{"done": true}', encoding="utf-8")
            with mock.patch.object(self.prepare, "download_progress_file", return_value=target):
                returned = self.prepare.begin_download("RNAseq")
            self.assertEqual(returned, target)
            self.assertFalse(target.exists())

    def test_the_command_points_at_the_script_and_its_progress_file(self):
        progress = Path("/tmp/download-RNAseq.json")
        command = self.prepare.download_command("RNAseq", "GSE1", ["SRR1", "SRR2"], progress)
        self.assertEqual(
            command[1], str(self.prepare.env_paths.SCRIPT_ROOT / "common" / "fetch_sra.py")
        )
        self.assertEqual(command[command.index("--srr") + 1], "SRR1")
        self.assertEqual(command[command.index("--progress-file") + 1], str(progress))

    def test_metadata_alone_unlocks_the_download_button(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "RNAseq" / "GSE1"
            path.mkdir(parents=True)
            (path / "metadata.csv").write_text(
                "sample,sample_name,group\n"
                "S1,S1_1,control\n"
                "S2,S2_1,control\n"
                "S3,S3_1,treatA\n"
                "S4,S4_1,treatA\n",
                encoding="utf-8",
            )
            with mock.patch.object(self.prepare.datasets, "DATA_ROOT", Path(root)):
                self.assertTrue(self.prepare.prepared("RNAseq", "GSE1"))
                self.assertEqual(self.prepare.datasets.fastq_table("RNAseq", "GSE1"), {})

    def test_an_incomplete_table_keeps_the_download_button_locked(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "RNAseq" / "GSE1"
            path.mkdir(parents=True)
            (path / "metadata.csv").write_text(
                "sample,sample_name,group\nS1,S1_1,\nS2,S2_1,\n", encoding="utf-8"
            )
            with mock.patch.object(self.prepare.datasets, "DATA_ROOT", Path(root)):
                self.assertFalse(self.prepare.prepared("RNAseq", "GSE1"))

    def test_names_are_not_blamed_for_fastq_that_is_not_downloaded_yet(self):
        with mock.patch.object(self.prepare.datasets, "fastq_table", return_value={}):
            self.assertEqual(
                self.prepare.unmatched_fastq(
                    "ChIPseq", "GSE1", [{"sample": "SRR1", "input_sample": "SRR0"}]
                ),
                [],
            )

    def test_the_run_list_is_read_in_the_order_it_was_written(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "GSE1"
            (path / "geo").mkdir(parents=True)
            (path / "geo" / "SRR.txt").write_text("SRR9\nSRR3\n\n", encoding="utf-8")
            with mock.patch.object(self.prepare.datasets, "directory", return_value=path):
                self.assertEqual(
                    self.prepare.datasets.srr_list("RNAseq", "GSE1"), ["SRR9", "SRR3"]
                )

    def test_a_run_started_by_an_older_page_counts_the_fastq_on_disk(self):
        with mock.patch.object(
            self.prepare, "downloaded", return_value=(1, ["SRR2"])
        ), mock.patch.object(self.prepare, "expected_bytes", return_value=1000):
            view = self.prepare.file_view("RNAseq", "GSE1", ["SRR1", "SRR2"])
        self.assertAlmostEqual(view["overall"], 0.5)
        self.assertFalse(view["done"])
        self.assertEqual(view["rows"][0][2], "进行中")
        self.assertEqual([row["状态"] for row in view["runs"]], ["完成", "下载中"])


class SraFetchProgressTest(unittest.TestCase):
    """The script reports the slowest run's stage and the bytes on disk."""

    def setUp(self):
        from common import fetch_sra

        self.fetch = fetch_sra

    def test_the_slowest_run_sets_the_stage(self):
        with tempfile.TemporaryDirectory() as cache:
            progress = self.fetch.Progress(
                None,
                "RNAseq",
                "demo",
                ["SRR1", "SRR2"],
                {"SRR1": 100, "SRR2": 100},
                Path(cache),
            )
            progress.enter("SRR1", "download")
            progress.prefetched("SRR1", 100)
            progress.enter("SRR2", "download")
            state = progress.snapshot()
        self.assertEqual(state["step"], "download")
        self.assertAlmostEqual(state["step_progress"], 0.5)
        self.assertIn("100.0 B / 200.0 B", state["detail"])
        self.assertEqual(state["finished"], 0)

    def test_a_finished_run_counts_in_the_done_column(self):
        with tempfile.TemporaryDirectory() as cache:
            progress = self.fetch.Progress(
                None, "RNAseq", "demo", ["SRR1"], {"SRR1": 100}, Path(cache)
            )
            progress.finish("SRR1")
            state = progress.snapshot()
        self.assertEqual(state["step"], "compress")
        self.assertEqual(state["finished"], 1)
        self.assertEqual(state["step_progress"], 1.0)

    def test_a_failed_run_is_not_counted_as_finished(self):
        with tempfile.TemporaryDirectory() as cache:
            progress = self.fetch.Progress(
                None, "RNAseq", "demo", ["SRR1"], {}, Path(cache)
            )
            progress.fail("SRR1", "prefetch failed")
            state = progress.snapshot()
        self.assertEqual(state["finished"], 0)
        self.assertFalse(state["done"])

    def test_the_file_is_replaced_atomically_and_marks_the_end(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "progress" / "download-RNAseq.json"
            progress = self.fetch.Progress(
                path, "RNAseq", "demo", ["SRR1"], {}, Path(root)
            )
            progress.start()
            progress.finish("SRR1")
            progress.close(error="")
            state = json.loads(path.read_text(encoding="utf-8"))
            self.assertTrue(state["done"])
            self.assertEqual(state["finished"], 1)
            self.assertFalse(list(path.parent.glob("*.tmp")))

    def test_bytes_written_to_the_part_file_count_as_downloaded(self):
        with tempfile.TemporaryDirectory() as cache:
            folder = Path(cache) / "SRR1"
            folder.mkdir()
            (folder / "SRR1.sra.part").write_bytes(b"x" * 40)
            progress = self.fetch.Progress(
                None, "RNAseq", "demo", ["SRR1"], {"SRR1": 100}, Path(cache)
            )
            progress.enter("SRR1", "download")
            state = progress.snapshot()
        self.assertEqual(state["downloaded"], 40)
        self.assertAlmostEqual(state["step_progress"], 0.4)

    def test_the_direct_fetch_renames_the_part_file(self):
        def fake_run(command, check):
            Path(command[command.index("-o") + 1]).write_bytes(b"x" * 10)

        with tempfile.TemporaryDirectory() as cache, mock.patch.object(
            self.fetch.shutil, "which", return_value="/usr/bin/curl"
        ), mock.patch.object(self.fetch.subprocess, "run", side_effect=fake_run) as run:
            path = self.fetch.https_sra(Path(cache), "SRR1")
            self.assertEqual(path, Path(cache) / "SRR1" / "SRR1.sra")
            self.assertTrue(path.is_file())
            self.assertFalse((Path(cache) / "SRR1" / "SRR1.sra.part").exists())
        command = run.call_args[0][0]
        self.assertEqual(command[0], "curl")
        self.assertIn("https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1/SRR1", command)

    def test_a_failed_direct_fetch_leaves_the_prefetch_fallback(self):
        with tempfile.TemporaryDirectory() as cache, mock.patch.object(
            self.fetch.shutil, "which", return_value="/usr/bin/curl"
        ), mock.patch.object(
            self.fetch.subprocess,
            "run",
            side_effect=self.fetch.subprocess.CalledProcessError(22, ["curl"]),
        ):
            result = self.fetch.https_sra(Path(cache), "SRR1")
        self.assertIsNone(result)
        self.assertFalse((Path(cache) / "SRR1" / "SRR1.sra.part").exists())


if __name__ == "__main__":
    unittest.main()
