import os
import argparse
import subprocess
import glob
from concurrent.futures import ThreadPoolExecutor

def parse_args():
    parser = argparse.ArgumentParser(description="Step 3.2: StringTie Quantification")
    parser.add_argument("--inputdir", required=True, help="Directory containing hisat2 output")
    parser.add_argument("--outputdir", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=4, help="Threads per task")
    parser.add_argument("--maxparallel", type=int, default=4, help="Max parallel tasks")
    parser.add_argument("--gtf_file", type=str, required=True, help="Specific GTF file to use")
    return parser.parse_args()

def run_stringtie(hisat2_path, args, rel_dir, sample_id):
    gtf_file = args.gtf_file
    
    print(f"[{sample_id}] Starting StringTie quantification using GTF: {gtf_file}...")
    dedup_bam = os.path.join(hisat2_path, f"{sample_id}.dedup.bam")

    if not os.path.exists(dedup_bam):
        print(f"[{sample_id}] Warning: dedup.bam not found at {dedup_bam}. Skipping.")
        return False

    # outputdir/rel_dir/mRNA/stringtie/sample_id
    out_dir = os.path.join(args.outputdir, rel_dir, "mRNA/stringtie", sample_id)
    os.makedirs(out_dir, exist_ok=True)
    out_gtf = os.path.join(out_dir, "transcripts.gtf")
    gene_abund = os.path.join(out_dir, "gene_abund.tab")
    
    stringtie_cmd = [
        "stringtie", "-p", str(args.threads), "-e", "-B",
        "-G", gtf_file,
        "-A", gene_abund,
        "-o", out_gtf,
        dedup_bam
    ]
    try:
        subprocess.run(stringtie_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[{sample_id}] StringTie failed: {e}")
        return False
        
    return True

def main():
    args = parse_args()
    
    tasks = []
    print(f"[Info] Scanning {args.inputdir} for hisat2file directories...")
    for root, dirs, files in os.walk(args.inputdir):
        if "hisat2file" in dirs:
            rel_dir = os.path.relpath(root, args.inputdir)
            if rel_dir == ".":
                rel_dir = ""
            hisat2_base = os.path.join(root, "hisat2file")
            print(f"[Info] Found hisat2file dir at: {hisat2_base}")
            for sid in os.listdir(hisat2_base):
                sid_path = os.path.join(hisat2_base, sid)
                if os.path.isdir(sid_path):
                    tasks.append((sid_path, args, rel_dir, sid))

    if not tasks:
        print(f"[Warning] No samples found in {args.inputdir} for StringTie quantification.")

    with ThreadPoolExecutor(max_workers=args.maxparallel) as executor:
        futures = [executor.submit(run_stringtie, *t) for t in tasks]
        for f in futures:
            f.result()

    with open(os.path.join(args.outputdir, "Quant_finished.txt"), "w") as f:
        f.write("Quantification finished.\n")

if __name__ == "__main__":
    main()
