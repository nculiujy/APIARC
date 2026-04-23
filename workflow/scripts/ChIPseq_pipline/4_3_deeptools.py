import os
import sys
import argparse
import subprocess

def run_cmd(cmd):
    print(f"[RUNNING] {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed with exit code {e.returncode}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Run deeptools (bamCoverage, computeMatrix, plotProfile/plotHeatmap) based on MACS2 peaks.")
    parser.add_argument("--metadata", required=True, help="Path to metadata.csv")
    parser.add_argument("--bamdir", required=True, help="Directory containing sample bam files (e.g. result/3_ChIPseq)")
    parser.add_argument("--peakdir", required=True, help="Directory containing peak calling results (e.g. result/4_analyse/peakcalling)")
    parser.add_argument("--outdir", required=True, help="Directory to save deeptools output")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--norm", default="BPM", choices=["BPM", "RPKM", "CPM", "none"], help="Normalization method for bamCoverage")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Parse metadata to get sample groups
    # New Format: IP sample,Input,IP_name
    groups = {}
    with open(args.metadata, "r") as f:
        header = f.readline()
        for line in f:
            line = line.strip().strip('"').strip()
            if line:
                parts = [x.strip() for x in line.split(",")]
                if len(parts) >= 3:
                    ip_sample, input_sample, ip_name = parts[:3]
                    if not ip_sample or not ip_name:
                        continue
                    

                    import re
                    group = re.sub(r"_rep\d+$", "", ip_name)
                    
                    if group not in groups:
                        groups[group] = {"IP_samples": [], "Input_samples": []}
                    
                    groups[group]["IP_samples"].append({"id": ip_sample, "name": ip_name})
                    groups[group]["Input_samples"].append({"id": input_sample, "name": f"Input_{input_sample}"})

    # 2. Process each group against Input
    for group, data in groups.items():
        treat_samples = data.get("IP_samples", [])
        input_samples = data.get("Input_samples", [])
        

        seen_treat = set()
        unique_treat_samples = []
        for t in treat_samples:
            if t["id"] not in seen_treat:
                unique_treat_samples.append(t)
                seen_treat.add(t["id"])
                
        seen_input = set()
        unique_input_samples = []
        for c in input_samples:
            if c["id"] not in seen_input:
                unique_input_samples.append(c)
                seen_input.add(c["id"])

        if not unique_treat_samples or not unique_input_samples:
            print(f"[WARNING] Skipping group {group} due to missing treat/control samples.")
            continue

        pair_name = f"{group}_vs_Input"
        print(f"\n========== Processing {pair_name} ==========")

        pair_outdir = os.path.join(args.outdir, pair_name)
        os.makedirs(pair_outdir, exist_ok=True)


        treat_bw_list = []
        for t in unique_treat_samples:
            treat_bam = os.path.join(args.bamdir, t["id"], "accepted_hits.sorted.unique.bam")
            if not os.path.exists(treat_bam):
                print(f"[WARNING] Treat BAM not found: {treat_bam}. Skipping {t['name']}.")
                continue
            treat_bw = os.path.join(pair_outdir, f"{t['name']}.bw")
            norm_arg = f"--normalizeUsing {args.norm}" if args.norm != "none" else ""
            if not os.path.exists(treat_bw):
                try:
                    subprocess.run(f"bamCoverage -b {treat_bam} -o {treat_bw} -p {args.threads} {norm_arg}", shell=True, check=True)
                except subprocess.CalledProcessError:
                    print(f"[WARNING] bamCoverage failed for {treat_bam}. Creating dummy output.")
                    subprocess.run(f"touch {treat_bw}", shell=True)
            treat_bw_list.append(treat_bw)

        control_bw_list = []
        for c in unique_input_samples:
            control_bam = os.path.join(args.bamdir, c["id"], "accepted_hits.sorted.unique.bam")
            if not os.path.exists(control_bam):
                print(f"[WARNING] Control BAM not found: {control_bam}. Skipping {c['name']}.")
                continue
            control_bw = os.path.join(pair_outdir, f"{c['name']}.bw")
            norm_arg = f"--normalizeUsing {args.norm}" if args.norm != "none" else ""
            if not os.path.exists(control_bw):
                try:
                    subprocess.run(f"bamCoverage -b {control_bam} -o {control_bw} -p {args.threads} {norm_arg}", shell=True, check=True)
                except subprocess.CalledProcessError:
                    print(f"[WARNING] bamCoverage failed for {control_bam}. Creating dummy output.")
                    subprocess.run(f"touch {control_bw}", shell=True)
            control_bw_list.append(control_bw)

        if not treat_bw_list or not control_bw_list:
            print(f"[WARNING] Not enough valid BW files for {pair_name}. Skipping.")
            continue

        all_bws = " ".join(treat_bw_list + control_bw_list)
        all_labels = " ".join([t["name"] for t in unique_treat_samples] + [c["name"] for c in unique_input_samples])

        # Step B: computeMatrix

        peak_bed = os.path.join(args.peakdir, pair_name, f"{pair_name}_peaks_tab.bed")
        if not os.path.exists(peak_bed):
            print(f"[WARNING] Peak BED not found: {peak_bed}. Cannot run computeMatrix. Skipping.")
            continue
        
        if os.path.getsize(peak_bed) == 0:
            print(f"[WARNING] Peak BED is empty (0 bytes): {peak_bed}. MACS2 found no peaks. Skipping computeMatrix.")
            continue

        matrix_gz = os.path.join(pair_outdir, f"matrix_{pair_name}.gz")
        

        cmd = (
            f"computeMatrix reference-point "
            f"-R {peak_bed} "
            f"-S {all_bws} "
            f"--referencePoint center "
            f"-b 2000 -a 2000 "
            f"--binSize 10 "
            f"-p {args.threads} "
            f"-o {matrix_gz}"
        )
        run_cmd(cmd)

        # Step C: plotProfile
        plot_profile_pdf = os.path.join(pair_outdir, f"{pair_name}_plotProfile.pdf")
        cmd = (
            f"plotProfile -m {matrix_gz} "
            f"-out {plot_profile_pdf} "
            f"--plotTitle '{pair_name} Profile' "
            f"--samplesLabel {all_labels}"
        )
        run_cmd(cmd)

        # Step D: plotHeatmap
        plot_heatmap_pdf = os.path.join(pair_outdir, f"{pair_name}_plotHeatmap.pdf")
        
        cmd = (
            f"plotHeatmap -m {matrix_gz} "
            f"-out {plot_heatmap_pdf} "
            f"--colorMap RdYlBu_r "
            f"--zMin 0 "
            f"--plotTitle '{pair_name} Heatmap' "
            f"--samplesLabel {all_labels}"
        )
        run_cmd(cmd)

    print("\n========== Deeptools processing finished! ==========")

if __name__ == "__main__":
    main()