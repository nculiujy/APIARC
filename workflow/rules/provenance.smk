
PROVENANCE_MANIFEST = f"{PROVENANCE_ROOT}/run_manifest.yaml"

VERSION_TOOLS = {
    "snakemake": ["snakemake", "--version"],
    "bowtie2": ["bowtie2", "--version"],
    "samtools": ["samtools", "--version"],
    "macs2": ["macs2", "--version"],
    "bedtools": ["bedtools", "--version"],
    "deeptools": ["bamCoverage", "--version"],
    "subread_featureCounts": ["featureCounts", "-v"],
    "meme_ame": ["ame", "--version"],
    "hisat2": ["hisat2", "--version"],
    "stringtie": ["stringtie", "--version"],
    "trim_galore": ["trim_galore", "--version"],
    "fastp": ["fastp", "--version"],
    "R": ["R", "--version"],
    "python": ["python", "--version"],
}


rule provenance_versions:
    input:
        validation=VALIDATION_REPORT
    output:
        versions=f"{PROVENANCE_ROOT}/software_versions.tsv"
    params:
        provenance_root=PROVENANCE_ROOT
    run:
        import shutil
        import subprocess

        tools = ""
        for name, command in sorted(VERSION_TOOLS.items()):
            if shutil.which(command[0]) is None:
                tools += f"{name}\tunavailable\t\n"
                continue
            try:
                completed = subprocess.run(
                    command, capture_output=True, text=True, timeout=60
                )
                text = (completed.stdout or completed.stderr or "").strip()
                first = text.splitlines()[0] if text else ""
            except (OSError, subprocess.SubprocessError) as error:
                first = f"error: {error}"
            tools += f"{name}\tavailable\t{first}\n"
        with open(output.versions, "w") as handle:
            handle.write("tool\tstatus\tversion\n")
            handle.write(tools)


rule provenance_environment:
    input:
        versions=f"{PROVENANCE_ROOT}/software_versions.tsv",
        chip_metadata=CHIP_METADATA,
        rna_metadata=RNA_METADATA,
        chip_config=CONFIG_FILES
    output:
        manifest=PROVENANCE_MANIFEST
    params:
        root=PROVENANCE_ROOT
    threads: 1
    run:
        import shutil
        import subprocess
        import yaml

        root = params.root
        for source, target in (
            (input.chip_metadata, f"{root}/ChIPseq_metadata.snapshot.csv"),
            (input.rna_metadata, f"{root}/RNAseq_metadata.snapshot.csv"),
        ):
            shutil.copyfile(source, target)

        with open(f"{root}/resolved_config.yaml", "w") as handle:
            yaml.safe_dump(config, handle, sort_keys=True)

        if shutil.which("conda"):
            try:
                completed = subprocess.run(
                    ["conda", "list", "--explicit"], capture_output=True, text=True, timeout=300
                )
                with open(f"{root}/conda_explicit.txt", "w") as handle:
                    handle.write(completed.stdout)
            except (OSError, subprocess.SubprocessError):
                pass

        git_lines = []
        if shutil.which("git"):
            for label, command in (
                ("commit", ["git", "rev-parse", "HEAD"]),
                ("describe", ["git", "describe", "--tags", "--always", "--dirty"]),
                ("status", ["git", "status", "--porcelain"]),
            ):
                try:
                    completed = subprocess.run(
                        command, capture_output=True, text=True, timeout=60
                    )
                    git_lines.append(f"{label}\t{completed.stdout.strip() or 'clean'}")
                except (OSError, subprocess.SubprocessError):
                    git_lines.append(f"{label}\tunavailable")
        with open(f"{root}/git_state.tsv", "w") as handle:
            handle.write("field\tvalue\n")
            handle.write("\n".join(git_lines) + "\n")

        with open(output.manifest, "w") as handle:
            yaml.safe_dump(
                {
                    "config": config,
                    "annotation": annotation_manifest(config),
                    "files": sorted(str(path) for path in input),
                    "mode": CHIP_MODE,
                    "marks": CHIP_MARKS,
                    "comparisons": [
                        {"mark": mark, "reference": reference, "treatment": treatment}
                        for mark, reference, treatment in CHIP_COMPARISONS
                    ],
                    "snakemake_report": (
                        f"run `snakemake --report {PROVENANCE_ROOT}/report.html` "
                        "for the full provenance report"
                    ),
                },
                handle,
                sort_keys=True,
            )
