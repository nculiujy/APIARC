
rule rnaseq_3_rc_pipline:
    input:
        qc_flag=f"{RNA_ROOT}/2_cleandata/finished.txt"
    output:
        gene_matrix=f"{RNA_ROOT}/3_align_quant/gene_count_matrix.csv",
        transcript_matrix=f"{RNA_ROOT}/3_align_quant/transcript_count_matrix.csv",
        sample_list=f"{RNA_ROOT}/3_align_quant/sample_list.txt",
        finished_flag=f"{RNA_ROOT}/3_align_quant/finished.txt"
    params:
        out_dir=f"{RNA_ROOT}/3_align_quant",
        input_dir=f"{RNA_ROOT}/2_cleandata",
        align_script="workflow/scripts/RNAseq_pipline/3_1_Align.py",
        quant_script="workflow/scripts/RNAseq_pipline/3_2_Quant.py",
        prepde_script="workflow/scripts/RNAseq_pipline/3_3_prepDE.py",
        metadata_file=RNA_METADATA,
        hisat2_index=rna_index(ASSEMBLY),
        gtf_file=gtf(ASSEMBLY),
        max_parallel=config.get("max_parallel", 4)
    threads: config.get("threads", 8)
    log:
        "logs/RNAseq/3_RC_pipline.log"
    conda:
        "../../envs/RNAseq_3_RC_pipline.yaml"
    shell:
        """
        mkdir -p {params.out_dir}

        python {params.align_script} \\
            --inputdir {params.input_dir} \\
            --outputdir {params.out_dir} \\
            --threads {threads} \\
            --maxparallel {params.max_parallel} \\
            --hisat2_index {params.hisat2_index} >> {log} 2>&1

        python {params.quant_script} \\
            --inputdir {params.out_dir} \\
            --outputdir {params.out_dir} \\
            --threads {threads} \\
            --maxparallel {params.max_parallel} \\
            --gtf_file {params.gtf_file} >> {log} 2>&1

        SAMPLE_LIST="{output.sample_list}"
        : > $SAMPLE_LIST

        awk -F ',' 'NR>1 {{print $1" "$2}}' {params.metadata_file} | while read SRR NAME; do
            GTF_PATH=$(find {params.out_dir} -name "transcripts.gtf" | grep "/$SRR/" | head -n 1 || true)
            if [ -n "$GTF_PATH" ]; then
                # The matrix is keyed by the run ID, not the group label: two
                # replicates share a sample_name and duplicate column names
                # abort DESeq2 on duplicate row.names.
                echo "${{SRR}} $GTF_PATH" >> $SAMPLE_LIST
            else
                echo "[ERROR] no transcripts.gtf for sample $SRR" | tee -a {log}
            fi
        done

        cat $SAMPLE_LIST >> {log}

        if [ ! -s $SAMPLE_LIST ]; then
            echo "[ERROR] sample_list.txt is empty; quantification produced nothing" | tee -a {log}
            exit 1
        fi

        python {params.prepde_script} \\
            -i $SAMPLE_LIST \\
            -g {output.gene_matrix} \\
            -t {output.transcript_matrix} >> {log} 2>&1

        touch {output.finished_flag}
        """
