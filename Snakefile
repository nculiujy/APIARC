configfile: "config/config.yaml"

import os
import sys

sys.path.insert(0, os.path.join(workflow.basedir, "workflow", "scripts"))

from common import workspace
from common.annotation import manifest as annotation_manifest
from common.annotation import resolve as resolve_annotation
from common.metadata import (
    comparisons as chip_comparisons,
    groups as chip_groups,
    read_chip_metadata,
    resolve_mode,
)
from common.resources import (
    assembly as resource_assembly,
    chip_index,
    gene_bed,
    genome_fasta,
    gtf,
    rna_index,
)

# One annotation for the whole project: RNA-seq and ChIP-seq cannot be run
# against two different reference versions, and the version/source pair is
# spelled into every result path so a switch never overwrites earlier output.
ANNO = resolve_annotation(config)
ASSEMBLY = ANNO["assembly"]
SPECIES = ANNO["species"]
ANNO_ID = ANNO["id"]

RNA_DATASET = workspace.default_name(config, "RNAseq")
CHIP_DATASET = workspace.default_name(config, "ChIPseq")
JOINT_ID = workspace.integrated_id(RNA_DATASET, CHIP_DATASET)

CHIP_METADATA = workspace.metadata_from_config(config, "ChIPseq")
RNA_METADATA = workspace.metadata_from_config(config, "RNAseq")
# The comparison pairs live beside the metadata, so a ``chip.metadata``
# override keeps working without a second setting to keep in sync.
CHIP_COMPARISON_FILE = os.path.join(
    os.path.dirname(CHIP_METADATA), workspace.COMPARISON_NAME
)

CONFIG_FILES = list(workflow.configfiles) or ["config/config.yaml"]

# ``--configfile`` files are merged on top of config/config.yaml, so no single
# input file holds what the rules actually see. Rules that re-read the config
# (peak calling, annotation) take this parse-time dump instead.
RUN_CONFIG = os.path.join(
    ".apiarc-runtime", "run-config", f"{ANNO_ID}__{RNA_DATASET}__{CHIP_DATASET}.yaml"
)


def check_descriptor(omics: str, name: str) -> None:
    """Refuse to align reads prepared against another assembly."""
    declared = str(workspace.read_descriptor(omics, name).get("assembly") or "").strip()
    if not declared:
        return
    if resource_assembly(declared) != ASSEMBLY:
        raise ValueError(
            f"{omics} dataset '{name}' was prepared against {declared} but the project "
            f"annotation is {ANNO_ID}; re-import the reads or switch the annotation "
            "in the annotation build page, both omics must use the same version"
        )


check_descriptor("RNAseq", RNA_DATASET)
check_descriptor("ChIPseq", CHIP_DATASET)


def primary_config_file():
    import yaml

    os.makedirs(os.path.dirname(RUN_CONFIG), exist_ok=True)
    with open(RUN_CONFIG, "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=True)
    return RUN_CONFIG


RNA_ROOT = f"result/RNAseq/{RNA_DATASET}/{ANNO_ID}"
CHIP_ROOT = f"result/ChIPseq/{CHIP_DATASET}/{ANNO_ID}"
JOINT_ROOT = f"result/Integrated/{JOINT_ID}/{ANNO_ID}"
PROVENANCE_ROOT = f"result/provenance/{JOINT_ID}/{ANNO_ID}"

PEAK_ROOT = f"{CHIP_ROOT}/4_peaks"
PEAK_DIR = f"{PEAK_ROOT}/peakcalling"
ANNO_DIR = f"{PEAK_ROOT}/annotation"
DEEPT_DIR = f"{PEAK_ROOT}/deeptools"
DIFF_ROOT = f"{CHIP_ROOT}/5_diffbind"
STATS_ROOT = f"{JOINT_ROOT}/0_stats"

CHIP_MODULES = config.get("ChIPseq_modules", {})
RNA_MODULES = config.get("RNAseq_modules", {})
INTEGRATED_MODULES = config.get("Integrated_modules", {})

try:
    CHIP_SAMPLES = read_chip_metadata(CHIP_METADATA)
except (OSError, ValueError) as error:
    if any(CHIP_MODULES.get(name, True) for name in ("1_prepare", "2_QC", "3_CC", "4_peak_result")):
        raise
    CHIP_SAMPLES = []

CHIP_GROUPS = chip_groups(CHIP_SAMPLES)
CHIP_MARKS = sorted({s.mark for s in CHIP_SAMPLES})
CHIP_MODE = resolve_mode(config, CHIP_SAMPLES, CHIP_COMPARISON_FILE)
CHIP_COMPARISONS = (
    chip_comparisons(CHIP_SAMPLES, config, CHIP_COMPARISON_FILE)
    if CHIP_MODE == "differential"
    else []
)

# --- RNA contrast for every ChIP comparison --------------------------------
# A ChIP comparison names its time point in the group ("H3K27ac_TGFb12h");
# the RNA metadata names the same time point in its sample_name column. The
# joint stages used to read one fixed file, DEG_result.csv, which is a copy of
# the first contrast only - so 12 h binding was tested against 2 h expression.
RNA_DESEQ_DIR = f"result/RNAseq/{RNA_DATASET}/{ANNO_ID}/4_deseq"
RNA_DEG_FALLBACK = f"{RNA_DESEQ_DIR}/DEG_result.csv"


def rna_group_by_condition():
    import csv

    mapping = {}
    try:
        with open(RNA_METADATA, newline="") as handle:
            for row in csv.DictReader(handle):
                name = (row.get("sample_name") or "").strip()
                group = (row.get("group") or "").strip()
                if name and group:
                    mapping.setdefault(name, group)
    except OSError:
        pass
    return mapping


RNA_GROUP_BY_CONDITION = rna_group_by_condition()


def condition_of(group):
    """H3K27ac_TGFb12h -> TGFb12h, the RNA sample name the group encodes.

    The comparison tuples carry group names such as H3K27ac_TGFb12h, so the
    mark is not available as a separate token to strip. The RNA metadata does
    know the time points, and the longest suffix that matches one of them is
    the condition; taking the part after the first underscore would break on a
    condition that itself contains one.
    """
    matches = [
        name for name in RNA_GROUP_BY_CONDITION
        if group == name or group.endswith(f"_{name}")
    ]
    return max(matches, key=len) if matches else group


def rna_table_for(condition):
    group = RNA_GROUP_BY_CONDITION.get(condition)
    if group and group != "control":
        candidate = f"{RNA_DESEQ_DIR}/DEG_result_{group}_vs_control.csv"
        if os.path.exists(candidate):
            return candidate
    return RNA_DEG_FALLBACK


def comparison_rna_pairs():
    """(comparison label, DESeq2 table of the same time point).

    The label has to be written exactly as 5_diffbind names it -
    treatment_vs_reference - or the gene table cannot find the contrast and
    silently falls back to another one.
    """
    return [
        (
            f"{treatment}_vs_{reference}",
            rna_table_for(condition_of(treatment)),
        )
        for _namespace, reference, treatment in CHIP_COMPARISONS
    ]


def comparison_rna_tables():
    tables = [path for _, path in comparison_rna_pairs()]
    return tables or [RNA_DEG_FALLBACK]


def rna_args_labelled():
    """--rna label=path, for the stage that joins expression to a comparison."""
    pairs = comparison_rna_pairs()
    if not pairs:
        return f"--rna {RNA_DEG_FALLBACK}"
    return " ".join(f"--rna {label}={path}" for label, path in pairs)


def rna_args_plain():
    """--rna path repeated once per contrast, for the expression universes."""
    seen = []
    for path in comparison_rna_tables():
        if path not in seen:
            seen.append(path)
    return " ".join(f"--rna {path}" for path in seen)


_requested_mode = str(config.get("integrated", {}).get("mode", "auto")).lower()
if _requested_mode == "differential" and not CHIP_COMPARISONS:
    raise ValueError(
        "integrated.mode is 'differential' but no comparison pairs two groups with "
        "at least chip.diff.min_replicates replicates each. Declare the pairs in "
        f"{CHIP_COMPARISON_FILE} or set integrated.mode to 'auto'/'descriptive'."
    )

include: "workflow/rules/0_validate.smk"
include: "workflow/rules/provenance.smk"

if RNA_MODULES.get("1_prepare", True):
    include: "workflow/rules/RNAseq_pipline/1_prepare.smk"
if RNA_MODULES.get("2_QC", True):
    include: "workflow/rules/RNAseq_pipline/2_QC.smk"

if RNA_MODULES.get("3_RC", True):
    include: "workflow/rules/RNAseq_pipline/3_RC_pipline.smk"

if RNA_MODULES.get("4_DEseq", True):
    include: "workflow/rules/RNAseq_pipline/4_DEseq.smk"

if CHIP_MODULES.get("1_prepare", True):
    include: "workflow/rules/ChIPseq_pipline/1_prepare.smk"

if CHIP_MODULES.get("2_QC", True):
    include: "workflow/rules/ChIPseq_pipline/2_QC.smk"

if CHIP_MODULES.get("3_CC", True):
    include: "workflow/rules/ChIPseq_pipline/3_CC_pipline.smk"

if CHIP_MODULES.get("4_peak_result", True):
    include: "workflow/rules/ChIPseq_pipline/4_peak_result.smk"

DIFFBIND_ENABLED = (
    CHIP_MODULES.get("5_diffbind", True)
    and CHIP_MODULES.get("4_peak_result", True)
    and CHIP_MODE == "differential"
)
if DIFFBIND_ENABLED:
    include: "workflow/rules/ChIPseq_pipline/5_diffbind.smk"

STATS_ENABLED = INTEGRATED_MODULES.get("0_stats", True) and bool(CHIP_GROUPS) and bool(
    RNA_MODULES.get("4_DEseq", True)
)
if STATS_ENABLED:
    include: "workflow/rules/Integrated/0_stats.smk"

ENRICHMENT_ENABLED = INTEGRATED_MODULES.get("1_enrichment", True) and STATS_ENABLED
if ENRICHMENT_ENABLED:
    include: "workflow/rules/Integrated/1_enrichment.smk"

# The common-gene stage needs differential binding: in descriptive mode every
# region is "bound" and there is no diff_up/diff_down to intersect with the DEG
# list, so the rule would write an empty table. Declared before the motif include
# because the per-region TF-gene edges are restricted to the common genes and
# only exist when this stage runs.
COMMON_ENABLED = (
    INTEGRATED_MODULES.get("4_common", True)
    and STATS_ENABLED
    and DIFFBIND_ENABLED
)

MOTIF_ENABLED = INTEGRATED_MODULES.get("2_motif", True) and STATS_ENABLED
if MOTIF_ENABLED:
    include: "workflow/rules/Integrated/2_motif.smk"

ENHANCER_ENABLED = INTEGRATED_MODULES.get("3_enhancer", True) and STATS_ENABLED
if ENHANCER_ENABLED:
    include: "workflow/rules/Integrated/3_enhancer.smk"

if COMMON_ENABLED:
    include: "workflow/rules/Integrated/4_common.smk"


RNA_STAGES = ("1_prepare", "2_QC", "3_RC", "4_DEseq")
CHIP_STAGES = ("1_prepare", "2_QC", "3_CC", "4_peak_result", "5_diffbind")

RNA_TARGETS = {
    name: f"{RNA_ROOT}/{workspace.STAGE_DIRS['RNAseq'][name]}/finished.txt"
    for name in RNA_STAGES
}
CHIP_TARGETS = {
    name: f"{CHIP_ROOT}/{workspace.STAGE_DIRS['ChIPseq'][name]}/finished.txt"
    for name in CHIP_STAGES
}

CHIP_ACTIVE = {name: bool(CHIP_MODULES.get(name, True)) for name in CHIP_STAGES[:4]}
CHIP_ACTIVE["5_diffbind"] = DIFFBIND_ENABLED

INTEGRATED_TARGETS = {
    name: f"{JOINT_ROOT}/{workspace.INTEGRATED_STAGES[name]}/finished.txt"
    for name in workspace.INTEGRATED_STAGES
}
INTEGRATED_ACTIVE = {
    "0_stats": STATS_ENABLED,
    "1_enrichment": ENRICHMENT_ENABLED,
    "2_motif": MOTIF_ENABLED,
    "3_enhancer": ENHANCER_ENABLED,
    "4_common": COMMON_ENABLED,
}


def deepest_target(modules, order, targets):
    enabled = [name for name in order if modules.get(name, True)]
    if not enabled:
        return []
    return [targets[enabled[-1]]]


TARGETS = []
TARGETS += [PROVENANCE_MANIFEST]
TARGETS += deepest_target(RNA_MODULES, RNA_STAGES, RNA_TARGETS)
TARGETS += deepest_target(CHIP_ACTIVE, CHIP_STAGES, CHIP_TARGETS)
TARGETS += [
    target for name, target in INTEGRATED_TARGETS.items() if INTEGRATED_ACTIVE[name]
]


rule all:
    input:
        TARGETS
