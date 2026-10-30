
from common.metadata import comparison_samples, display_mark, uses_group_names

GROUP_SCHEMA = uses_group_names(CHIP_SAMPLES)
DIFF_CFG = config.get("chip", {}).get("diff", {})
DIFF_METHODS = ",".join(DIFF_CFG.get("methods", ["count", "presence", "stat"]))
DIFF_FDR = DIFF_CFG.get("fdr", 0.05)
DIFF_LFC = DIFF_CFG.get("lfc", 1.0)
DIFF_MIN_SUPPORT = int(DIFF_CFG.get("min_replicates", 2))
DIFF_LABEL_PERMUTATIONS = int(DIFF_CFG.get("label_permutations", 20))
DIFF_MAX_FALSE_POSITIVE = int(DIFF_CFG.get("max_false_positive", 5))


def comparison_parts(comparison):
    """``treatment_vs_reference`` -> (treatment, reference)."""
    treatment, _, reference = comparison.partition("_vs_")
    if not treatment or not reference:
        raise ValueError(f"comparison '{comparison}' is not in treatment_vs_reference form")
    return treatment, reference


def peaks_for_comparison(wildcards):
    """Per-replicate peaks of the two groups in this contrast."""
    treatment, reference = comparison_parts(wildcards.comparison)
    samples = comparison_samples(CHIP_SAMPLES, reference, treatment, wildcards.mark)
    return [f"{PEAK_DIR}/{s.sample}/{s.sample}_peaks.bed" for s in samples]


def comparison_label(wildcards):
    """Mark column written into ``differential_regions.tsv``.

    The five-column schema names the path namespace after the treatment group
    (``H3K27ac_veh``); reports want the mark part only.
    """
    return display_mark(wildcards.mark) if GROUP_SCHEMA else wildcards.mark


DIFF_COMPARISON_TARGETS = []
for _mark, _reference, _treatment in CHIP_COMPARISONS:
    _base = f"{DIFF_ROOT}/{_mark}/{_treatment}_vs_{_reference}"
    DIFF_COMPARISON_TARGETS.extend([
        f"{_base}/consensus/consensus_peaks.bed",
        f"{_base}/counts.tsv",
        f"{_base}/count_model.tsv",
        f"{_base}/presence.tsv",
        f"{_base}/stat_shift.tsv",
        f"{_base}/differential_regions.tsv",
        f"{_base}/method_overlap.tsv",
        f"{_base}/diff_summary.tsv",
        f"{_base}/negative_control.tsv",
    ])


rule diffbind_consensus:
    input:
        peaks=peaks_for_comparison
    output:
        bed=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/consensus/consensus_peaks.bed",
        support=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/consensus/replicate_support.tsv"
    params:
        script="workflow/scripts/ChIPseq_pipline/5_1_consensus_peaks.py",
        metadata=CHIP_METADATA,
        reference=lambda wildcards: comparison_parts(wildcards.comparison)[1],
        treatment=lambda wildcards: comparison_parts(wildcards.comparison)[0],
        outdir=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/consensus"
    wildcard_constraints:
        mark="[^/]+",
        comparison="[^/]+"
    threads: 2
    log:
        "logs/ChIPseq/5_1_consensus_{mark}_{comparison}.log"
    conda:
        "../../envs/ChIPseq_5_diffbind.yaml"
    shell:
        """
        python {params.script} \\
            --metadata {params.metadata} \\
            --peakdir {PEAK_DIR} \\
            --mark {wildcards.mark} \\
            --reference {params.reference} \\
            --treatment {params.treatment} \\
            --outdir {params.outdir} > {log} 2>&1
        """


rule diffbind_counts:
    input:
        consensus=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/consensus/consensus_peaks.bed",
        manifest=f"{PEAK_DIR}/bam_manifest.tsv"
    output:
        counts=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/counts.tsv"
    params:
        script="workflow/scripts/ChIPseq_pipline/5_2_count_matrix.py",
        metadata=CHIP_METADATA,
        reference=lambda wildcards: comparison_parts(wildcards.comparison)[1],
        treatment=lambda wildcards: comparison_parts(wildcards.comparison)[0]
    wildcard_constraints:
        mark="[^/]+",
        comparison="[^/]+"
    threads: config.get("threads", 8)
    log:
        "logs/ChIPseq/5_2_counts_{mark}_{comparison}.log"
    conda:
        "../../envs/ChIPseq_5_diffbind.yaml"
    shell:
        """
        python {params.script} \\
            --metadata {params.metadata} \\
            --manifest {input.manifest} \\
            --consensus {input.consensus} \\
            --mark {wildcards.mark} \\
            --reference {params.reference} \\
            --treatment {params.treatment} \\
            --out {output.counts} \\
            --threads {threads} > {log} 2>&1
        """


rule diffbind_count_model:
    input:
        counts=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/counts.tsv",
        support=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/consensus/replicate_support.tsv"
    output:
        result=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/count_model.tsv",
        normalised=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/normalized_counts.tsv"
    params:
        script="workflow/scripts/ChIPseq_pipline/5_3_diff_count.R",
        metadata=CHIP_METADATA,
        reference=lambda wildcards: comparison_parts(wildcards.comparison)[1],
        treatment=lambda wildcards: comparison_parts(wildcards.comparison)[0]
    wildcard_constraints:
        mark="[^/]+",
        comparison="[^/]+"
    threads: 2
    log:
        "logs/ChIPseq/5_3_count_{mark}_{comparison}.log"
    conda:
        "../../envs/ChIPseq_5_diffbind.yaml"
    shell:
        """
        Rscript {params.script} \\
            --counts {input.counts} \\
            --metadata {params.metadata} \\
            --support {input.support} \\
            --mark {wildcards.mark} \\
            --reference {params.reference} \\
            --treatment {params.treatment} \\
            --min-support {DIFF_MIN_SUPPORT} \\
            --fdr {DIFF_FDR} \\
            --lfc {DIFF_LFC} \\
            --out {output.result} \\
            --norm-out {output.normalised} > {log} 2>&1
        """


rule diffbind_presence:
    input:
        support=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/consensus/replicate_support.tsv"
    output:
        result=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/presence.tsv"
    params:
        script="workflow/scripts/ChIPseq_pipline/5_4_presence.py",
        metadata=CHIP_METADATA,
        reference=lambda wildcards: comparison_parts(wildcards.comparison)[1],
        treatment=lambda wildcards: comparison_parts(wildcards.comparison)[0]
    wildcard_constraints:
        mark="[^/]+",
        comparison="[^/]+"
    threads: 2
    log:
        "logs/ChIPseq/5_4_presence_{mark}_{comparison}.log"
    conda:
        "../../envs/ChIPseq_5_diffbind.yaml"
    shell:
        """
        python {params.script} \\
            --support {input.support} \\
            --metadata {params.metadata} \\
            --mark {wildcards.mark} \\
            --reference {params.reference} \\
            --treatment {params.treatment} \\
            --fdr {DIFF_FDR} \\
            --out {output.result} > {log} 2>&1
        """


rule diffbind_stat_shift:
    input:
        consensus=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/consensus/consensus_peaks.bed",
        peaks=peaks_for_comparison
    output:
        result=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/stat_shift.tsv"
    params:
        script="workflow/scripts/ChIPseq_pipline/5_5_stat_shift.py",
        metadata=CHIP_METADATA,
        reference=lambda wildcards: comparison_parts(wildcards.comparison)[1],
        treatment=lambda wildcards: comparison_parts(wildcards.comparison)[0]
    wildcard_constraints:
        mark="[^/]+",
        comparison="[^/]+"
    threads: 2
    log:
        "logs/ChIPseq/5_5_stat_{mark}_{comparison}.log"
    conda:
        "../../envs/ChIPseq_5_diffbind.yaml"
    shell:
        """
        python {params.script} \\
            --metadata {params.metadata} \\
            --peakdir {PEAK_DIR} \\
            --consensus {input.consensus} \\
            --mark {wildcards.mark} \\
            --reference {params.reference} \\
            --treatment {params.treatment} \\
            --fdr {DIFF_FDR} \\
            --out {output.result} > {log} 2>&1
        """


rule diffbind_merge:
    input:
        count_model=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/count_model.tsv",
        presence=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/presence.tsv",
        stat=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/stat_shift.tsv",
        consensus=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/consensus/consensus_peaks.bed"
    output:
        regions=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/differential_regions.tsv",
        overlap=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/method_overlap.tsv",
        summary=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/diff_summary.tsv"
    params:
        script="workflow/scripts/ChIPseq_pipline/5_6_merge_methods.py",
        methods=DIFF_METHODS,
        mark_label=comparison_label,
        outdir=f"{DIFF_ROOT}/{{mark}}/{{comparison}}"
    wildcard_constraints:
        mark="[^/]+",
        comparison="[^/]+"
    threads: 2
    log:
        "logs/ChIPseq/5_6_merge_{mark}_{comparison}.log"
    conda:
        "../../envs/ChIPseq_5_diffbind.yaml"
    shell:
        """
        python {params.script} \\
            --count {input.count_model} \\
            --presence {input.presence} \\
            --stat {input.stat} \\
            --consensus {input.consensus} \\
            --mark {params.mark_label} \\
            --comparison {wildcards.comparison} \\
            --methods {params.methods} \\
            --fdr {DIFF_FDR} \\
            --lfc {DIFF_LFC} \\
            --outdir {params.outdir} > {log} 2>&1
        """


rule diffbind_negative_control:
    input:
        counts=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/counts.tsv",
        support=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/consensus/replicate_support.tsv"
    output:
        result=f"{DIFF_ROOT}/{{mark}}/{{comparison}}/negative_control.tsv"
    params:
        script="workflow/scripts/ChIPseq_pipline/5_7_negative_control.R",
        metadata=CHIP_METADATA,
        reference=lambda wildcards: comparison_parts(wildcards.comparison)[1],
        treatment=lambda wildcards: comparison_parts(wildcards.comparison)[0]
    wildcard_constraints:
        mark="[^/]+",
        comparison="[^/]+"
    threads: 2
    log:
        "logs/ChIPseq/5_7_negative_{mark}_{comparison}.log"
    conda:
        "../../envs/ChIPseq_5_diffbind.yaml"
    shell:
        """
        Rscript {params.script} \\
            --counts {input.counts} \\
            --metadata {params.metadata} \\
            --support {input.support} \\
            --mark {wildcards.mark} \\
            --reference {params.reference} \\
            --treatment {params.treatment} \\
            --min-support {DIFF_MIN_SUPPORT} \\
            --permutations {DIFF_LABEL_PERMUTATIONS} \\
            --fdr {DIFF_FDR} \\
            --lfc {DIFF_LFC} \\
            --max-false-positive {DIFF_MAX_FALSE_POSITIVE} \\
            --out {output.result} > {log} 2>&1
        """


rule diffbind_finished:
    input:
        DIFF_COMPARISON_TARGETS
    output:
        flag=f"{DIFF_ROOT}/finished.txt"
    run:
        summaries = [path for path in input if path.endswith("diff_summary.tsv")]
        with open(output.flag, "w") as handle:
            handle.write("mark\tcomparison\tcategory\tn_regions\n")
            for path in summaries:
                with open(path) as source:
                    next(source, None)
                    for line in source:
                        handle.write(line)
