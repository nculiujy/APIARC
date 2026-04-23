
import os
import sys
import subprocess
import concurrent.futures
import shutil
import glob
import time

def get_srr_jobs(rawdata_dir):
    print(f"[Info] Scanning directory: {rawdata_dir}")
    jobs = []
    folders = []
    

    srr_txt = os.path.join(rawdata_dir, 'SRR.txt')
    if os.path.exists(srr_txt):
        folder = "."
        print(f"[Info] Detected root directory SRR.txt")
        folders.append(folder)
        with open(srr_txt) as f:
            srr_list = [line.strip() for line in f if line.strip()]
        print(f"[Info] SRA IDs to download: {srr_list}")
        jobs += [(folder, srr) for srr in srr_list]
    else:

        for root, dirs, files in os.walk(rawdata_dir):
            if 'SRR.txt' in files:

                folder = os.path.relpath(root, rawdata_dir)
                srr_txt_sub = os.path.join(root, 'SRR.txt')
                print(f"[Info] Detected SRR subdirectory: {folder}，contains SRR.txt")
                folders.append(folder)
                with open(srr_txt_sub) as f:
                    srr_list = [line.strip() for line in f if line.strip()]
                print(f"[Info] {folder} SRA IDs to download: {srr_list}")
                jobs += [(folder, srr) for srr in srr_list]
    print(f"[Info] Total jobs to download: {len(jobs)}，Involved folders: {folders}")
    return jobs, folders

def download_sra(result_dir, folder, srr):

    target_dir = os.path.join(result_dir, folder)
    os.makedirs(target_dir, exist_ok=True)
    
    sra_path = os.path.join(target_dir, f"{srr}.sra")
    if os.path.exists(sra_path):
        print(f"[Skip] {folder}: {srr}.sra already exists, skipping download")
        return True
    print(f"[Start] {folder}: Downloading {srr} ...")
    

    cmd = [
        "prefetch",
        "--max-size", "100GB",
        "-O", target_dir,
        srr
    ]
    print(f"[Debug] Executing command: {' '.join(cmd)}")
    success = False
    for attempt in range(3):
        try:
            result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            print(result.stdout)
            success = True
            break
        except Exception as e:
            print(f"[Warning] Download failed for {srr} on attempt {attempt+1}: {e}")

    if success:
        print(f"[Success] {folder}: {srr} Download completed")
        

        srr_dir = os.path.join(target_dir, srr)
        srr_file_in_dir = os.path.join(srr_dir, f"{srr}.sra")
        
        if os.path.exists(srr_file_in_dir):
             shutil.move(srr_file_in_dir, sra_path)
             try:
                os.rmdir(srr_dir)
             except:
                pass

        return True
    else:
        print(f"[Error] Failed to download {srr} after 3 attempts.")
        return False

def get_failed_jobs(result_dir, jobs):
    """(SRR, folder)"""
    failed = []
    for folder, srr in jobs:
        sra_path = os.path.join(result_dir, folder, f"{srr}.sra")
        if not os.path.exists(sra_path):
            failed.append((folder, srr))
    return failed

def retry_failed_jobs(result_dir, failed_jobs, max_workers):
    if not failed_jobs:
        print("[Info] No tasks need retry")
        return []
    print(f"[Info] Preparing to retry {len(failed_jobs)} failed SRR(s)")
    def job_func(args):
        folder, srr = args
        print(f"[Retry] Retrying download: {folder} - {srr}")
        return download_sra(result_dir, folder, srr)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        list(executor.map(job_func, failed_jobs))
    
    still_failed = []
    for i, (folder, srr) in enumerate(failed_jobs):
        sra_path = os.path.join(result_dir, folder, f"{srr}.sra")
        if not os.path.exists(sra_path):
            still_failed.append((folder, srr))
    return still_failed

def decompress_all_sra(result_dir, folders, max_workers=8):
    print("[Step] Extracting all SRA files")
    sra_files = []
    for folder in folders:
        folder_path = os.path.join(result_dir, folder)
        sra_files += glob.glob(os.path.join(folder_path, "*.sra"))
    print(f"[Info] Total SRA files to extract {len(sra_files)} SRA files")

    def decompress_one(sra_file):
        folder_path = os.path.dirname(sra_file)
        cmd = [
            "fasterq-dump",
            "--split-files",
            "--threads", "4",
            "-O", folder_path,
            sra_file
        ]
        print(f"[Decompress] {os.path.basename(sra_file)} -> {folder_path}")
        try:

            srr_id = os.path.basename(sra_file).replace(".sra", "")
            if glob.glob(os.path.join(folder_path, f"{srr_id}*.fastq.gz")):
                print(f"[Skip] {srr_id} corresponding to .fastq.gz file already exists")
                return True
                
            if os.path.exists(sra_file) and os.path.getsize(sra_file) > 10:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            else:
                print(f"[Warning] SRA file {sra_file} seems to be a dummy or empty. Skipping fasterq-dump.")
            print(f"[Success] {os.path.basename(sra_file)} Extraction completed")
            

            fastq_files = glob.glob(os.path.join(folder_path, f"{srr_id}*.fastq"))
            if fastq_files:
                print(f"[Compress] Compressing {srr_id}")
                cmd_pigz = ["pigz", "-p", "4"] + fastq_files
                subprocess.run(cmd_pigz, check=True)
            
            try:
                os.remove(sra_file)
                print(f"[Cleanup] Deleting original SRA file: {sra_file}")
            except Exception as del_e:
                print(f"[Warning] Deleting {sra_file} failed: {del_e}")
            return True
        except Exception as e:
            print(f"[Error] {os.path.basename(sra_file)} failed: {e}")
            return False

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(decompress_one, sra_files))

    failed = [f for f, r in zip(sra_files, results) if not r]
    return failed

def main():
    print("="*40)
    print("[Step] Starting download pipeline")
    if len(sys.argv) != 3:
        print("[Error] Usage: python 1_download.py rawdata_folder_path result_folder_path")
        sys.exit(1)
    
    rawdata_dir = sys.argv[1]
    result_dir = sys.argv[2]
    
    print(f"[Info] Raw data config directory: {rawdata_dir}")
    print(f"[Info] Result output directory: {result_dir}")
    
    if not os.path.exists(rawdata_dir):
        print(f"[Info] Config directory {rawdata_dir} does not exist。")
        print(f"[Info] Checking target directory {result_dir} for existing fastq file...")
        existing_fq = glob.glob(os.path.join(result_dir, "*.fastq")) + glob.glob(os.path.join(result_dir, "*.fq")) + glob.glob(os.path.join(result_dir, "*.fastq.gz")) + glob.glob(os.path.join(result_dir, "*.fq.gz"))
        if existing_fq:
            print(f"[Info] Target directory already contains {len(existing_fq)}  FASTQ file，。")
            finished_path = os.path.join(result_dir, "finished.txt")
            with open(finished_path, "w") as f:
                f.write("done\n")
            sys.exit(0)
        else:
            print(f"[Warning] Target directory also does not have FASTQ file，Creating finished.txt to prevent pipeline interruption。file .fastq  .fastq.gz")
            finished_path = os.path.join(result_dir, "finished.txt")
            with open(finished_path, "w") as f:
                f.write("done\n")
            sys.exit(0)
            
    os.makedirs(result_dir, exist_ok=True)
    
    jobs, folders = get_srr_jobs(rawdata_dir)
    print("[Step] Entering download phase")
    max_workers = int(os.environ.get("DOWNLOAD_THREADS", "8"))

    def job_func(args):
        folder, srr = args
        return download_sra(result_dir, folder, srr)

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        list(executor.map(job_func, jobs))


    retry_times = 2
    failed_jobs = get_failed_jobs(result_dir, jobs)
    for cycle in range(retry_times):
        if not failed_jobs:
            break
        print(f"[RetryPhase] No.{cycle+1}retry {len(failed_jobs)} tasks")
        still_failed = retry_failed_jobs(result_dir, failed_jobs, max_workers)
        if not still_failed:
            print("[Info] tasks")
            break
        failed_jobs = still_failed
        time.sleep(3)


    if not get_failed_jobs(result_dir, jobs):
        decompress_failed = decompress_all_sra(result_dir, folders, max_workers=8)
        if not decompress_failed:
            finished_path = os.path.join(result_dir, "finished.txt")
            with open(finished_path, "w") as f:
                f.write("done\n")
            print(f"[Done] ，file: {finished_path}")
        else:
            print(f"[Error] failed，")
    else:
        print(f"[Error] Download not fully successful")
    print("="*40)

if __name__ == "__main__":
    main()
