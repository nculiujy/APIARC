import os
import argparse
import subprocess
import glob
import re
from concurrent.futures import ThreadPoolExecutor

def parse_args():
    parser = argparse.ArgumentParser(description="Step 3.1: HISAT2 Alignment")
    parser.add_argument("--inputdir", required=True, help="Input directory containing .clean.fastq files")
    parser.add_argument("--outputdir", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=4, help="Threads per task")
    parser.add_argument("--maxparallel", type=int, default=4, help="Max parallel tasks")
    parser.add_argument("--hisat2_index", type=str, required=True, help="Path to HISAT2 index (e.g. .../grcm38/genome)")
    return parser.parse_args()

def run_hisat2(sample_id, fq1, fq2, args, rel_dir):
    index_dir = args.hisat2_index
    
    print(f"[{sample_id}] Starting HISAT2 alignment using index {index_dir}...")
    
    hisat2_dir = os.path.join(args.outputdir, rel_dir, "hisat2file", sample_id)
    os.makedirs(hisat2_dir, exist_ok=True)
    
    bam_file = os.path.join(hisat2_dir, f"{sample_id}.sorted.bam")
    dedup_bam = os.path.join(hisat2_dir, f"{sample_id}.dedup.bam")
    qc_log = os.path.join(hisat2_dir, "QC_results.log")

    if os.path.exists(dedup_bam):
        print(f"[{sample_id}] Dedup BAM already exists. Skipping.")
        return True

    sam_file = os.path.join(hisat2_dir, "accepted_hits.sam")
    
    # HISAT2 command
    hisat2_cmd = [
        "hisat2", "-x", index_dir, "-p", str(args.threads), "--dta",
        "--rg-id", sample_id, "--rg", f"SM:{sample_id}"
    ]
    if fq2:
        hisat2_cmd.extend(["-1", fq1, "-2", fq2])
    else:
        hisat2_cmd.extend(["-U", fq1])
    
    hisat2_cmd.extend(["-S", sam_file])
    
    try:
        with open(qc_log, "w") as log:
            subprocess.run(hisat2_cmd, check=True, stderr=log)
        
        subprocess.run(f"samtools view -bS {sam_file} | samtools sort -@ {args.threads} -o {bam_file}", shell=True, check=True)
        subprocess.run(f"samtools index {bam_file}", shell=True, check=True)
        
        if os.path.exists(sam_file):
            os.remove(sam_file)

        picard_jar = "workflow/resources/picard-2.18.2/picard.jar"
        dedup_cmd = [
            "java", "-Xmx15g", "-jar", picard_jar, "MarkDuplicates",
            f"I={bam_file}", f"O={dedup_bam}",
            f"METRICS_FILE={hisat2_dir}/{sample_id}.metrics",
            "REMOVE_DUPLICATES=true", "ASSUME_SORT_ORDER=coordinate"
        ]
        subprocess.run(dedup_cmd, check=True)
        subprocess.run(f"samtools index {dedup_bam}", shell=True, check=True)

        if os.path.exists(bam_file):
            os.remove(bam_file)
            os.remove(bam_file + ".bai")

    except subprocess.CalledProcessError as e:
        print(f"[{sample_id}] Alignment failed: {e}")
        return False
    
    return True

def main():
    args = parse_args()
            
    tasks = []
    processed_single_ends = set()
    
    for root, dirs, files in os.walk(args.inputdir):
        for f in files:
            if f.endswith((".clean.fastq", ".clean.fq", ".clean.fastq.gz", ".clean.fq.gz")):
                rel_dir = os.path.relpath(root, args.inputdir)
                if rel_dir == ".":
                    rel_dir = ""
                
                match_paired = re.search(r"(.+)_([12])\.clean\.(fastq|fq)(\.gz)?$", f)
                if match_paired:
                    sid = match_paired.group(1)
                    read_type = match_paired.group(2)
                    ext = match_paired.group(3)
                    gz = match_paired.group(4) or ""
                    if read_type == "1":
                        fq1 = os.path.join(root, f)
                        fq2_name = f.replace(f"_1.clean.{ext}{gz}", f"_2.clean.{ext}{gz}")
                        fq2 = os.path.join(root, fq2_name)
                        if not os.path.exists(fq2):
                            fq2 = None
                        tasks.append((sid, fq1, fq2, args, rel_dir))
                else:
                    match_single = re.search(r"(.+)\.clean\.(fastq|fq)(\.gz)?$", f)
                    if match_single:
                        sid = match_single.group(1)
                        if sid not in processed_single_ends:
                            fq1 = os.path.join(root, f)
                            fq2 = None
                            tasks.append((sid, fq1, fq2, args, rel_dir))
                            processed_single_ends.add(sid)

    print(f"Found {len(tasks)} samples to process.")
    
    def process_task(task):
        sid, fq1, fq2, args, rel_dir = task
        return run_hisat2(sid, fq1, fq2, args, rel_dir)

    with ThreadPoolExecutor(max_workers=args.maxparallel) as executor:
        futures = [executor.submit(process_task, t) for t in tasks]
        for f in futures:
            f.result()

    with open(os.path.join(args.outputdir, "Align_finished.txt"), "w") as f:
        f.write("Alignment finished.\n")

if __name__ == "__main__":
    main()
