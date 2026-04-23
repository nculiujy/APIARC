
import os
import sys
import glob
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

def run_command(cmd):
    """"""
    print(f"[Command] {cmd}")
    status = os.system(cmd)
    if status != 0:
        raise RuntimeError(f"failed: {cmd}")

def download_sra(srr_id, download_dir):
    """
     prefetch  SRA file，file。
    prefetch Creating SRR_ID ， .sra file。
     .sra file already exists，。
    """
    sra_dir = os.path.join(download_dir, srr_id)
    sra_file = os.path.join(sra_dir, f"{srr_id}.sra")
    

    if os.path.exists(sra_file):
        print(f"[Skip] SRA filealready exists, skipping download: {sra_file}")
    else:
        print(f"[Process] Starting download: {srr_id}")

        lock_file = f"{sra_file}.lock"
        if os.path.exists(lock_file):
            print(f"[Warning] file {lock_file}，Deleting")
            os.remove(lock_file)
            
        cmd_prefetch = f"prefetch {srr_id} --output-directory {download_dir} --max-size 100G"
        # Try a few times to handle network timeout
        success = False
        for attempt in range(3):
            print(f"[Command] {cmd_prefetch} (Attempt {attempt+1})")
            status = os.system(cmd_prefetch)
            if status == 0:
                success = True
                break
            else:
                print(f"[Warning] Download failed for {srr_id} on attempt {attempt+1}")
        
        if not success:
            print(f"[Error] SRA download failed after 3 attempts: {srr_id}")
            raise RuntimeError(f"Failed to download {srr_id}")
    
    return sra_file

def fastq_dump_sra(srr_id, download_dir, rawdata_dir):
    """
     fasterq-dump  .sra file FASTQ， pigz 。
    file。
    ： ( _1.fastq.gz  _2.fastq.gz)。
    ，file _1/_2 ，。
    """
    sra_dir = os.path.join(download_dir, srr_id)
    sra_file = os.path.join(sra_dir, f"{srr_id}.sra")
    
    if not os.path.exists(sra_file):
        print(f"[Error] corresponding to SRA file: {sra_file}, cannot proceed with fastq-dump。")
        return False


    fq1 = os.path.join(download_dir, f"{srr_id}_1.fastq.gz")
    fq2 = os.path.join(download_dir, f"{srr_id}_2.fastq.gz")
    fq_single = os.path.join(download_dir, f"{srr_id}.fastq.gz")
    
    if (os.path.exists(fq1) and os.path.exists(fq2)) or os.path.exists(fq_single):
        print(f"[Skip] FASTQ file already exists，: {srr_id}")
    else:
        print(f"[Process] Starting conversion: {srr_id}")

        cmd_dump = f"fasterq-dump --split-3 {sra_file} -O {download_dir} -e 8"
        if os.path.exists(sra_file):
            run_command(cmd_dump)
        else:
            print(f"[Error] SRA file {sra_file} not found.")
        

        print(f"[Process] Compressing: {srr_id}")


        import glob
        fastq_files = glob.glob(os.path.join(download_dir, f"{srr_id}*.fastq"))
        if fastq_files:
            cmd_pigz = f"pigz -p 8 " + " ".join(fastq_files)
            run_command(cmd_pigz)
        else:
            print(f"[Warning] Not found to compress fastq file: {srr_id}")
    
    return True

def process_sample(srr_id, download_dir, rawdata_dir):
    """
    Sample：
    1.  SRA file
    2.  SRA file FASTQ 
    """
    try:
        download_sra(srr_id, download_dir)
        fastq_dump_sra(srr_id, download_dir, rawdata_dir)
        


        

        sra_dir = os.path.join(download_dir, srr_id)
        if os.path.exists(sra_dir):
            import shutil
            print(f"[Cleanup] Deleting SRA file: {sra_dir}")
            shutil.rmtree(sra_dir, ignore_errors=True)
            
        return True
    except Exception as e:
        print(f"[Error] Sample {srr_id} failed: {e}")
        return False

def main():
    print("="*40)
    print("[Step] Starting download pipeline")
    
    if len(sys.argv) != 3:
        print("[Error] Usage: python 1_download.py rawdata_folder_path result_folder_path")
        sys.exit(1)
        
    rawdata_dir = sys.argv[1]
    result_dir = sys.argv[2]
    
    srr_txt_path = os.path.join(rawdata_dir, "SRR.txt")
    
    if not os.path.exists(rawdata_dir) or not os.path.exists(srr_txt_path):
        print(f"[Warning] Cannot find SRR file: {srr_txt_path}。")
        print(f"[Info] Checking target directory {result_dir} for existing fastq file...")
        import glob
        existing_fq = glob.glob(os.path.join(result_dir, "*.fastq")) + glob.glob(os.path.join(result_dir, "*.fq")) + glob.glob(os.path.join(result_dir, "*.fastq.gz")) + glob.glob(os.path.join(result_dir, "*.fq.gz"))
        if existing_fq:
            print(f"[Info] Target directory already contains {len(existing_fq)}  FASTQ file，。")
        else:
            print(f"[Warning] Target directory also does not have FASTQ file，Creating finished.txt to prevent pipeline interruption。file .fastq  .fastq.gz")
        
        os.makedirs(result_dir, exist_ok=True)
        finished_path = os.path.join(result_dir, "finished.txt")
        with open(finished_path, "w") as f:
            f.write("done\n")
        sys.exit(0)

    print(f"[Info] file: {srr_txt_path}")
    print(f"[Info] Output result directory is: {result_dir}")


    with open(srr_txt_path, 'r') as f:
        srr_ids = [line.strip() for line in f if line.strip()]
        
    if not srr_ids:
        print(f"[Error] file {srr_txt_path} not found in SRR ID。")
        sys.exit(0)

    print(f"[Info] Found a total of {len(srr_ids)} Sampleneed processing: {', '.join(srr_ids)}")
    
    download_dir = result_dir
    os.makedirs(download_dir, exist_ok=True)
    
    max_workers = 4
    print(f"[Info] Starting multiprocessing, max workers: {max_workers}")
    
    success_count = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_sample, srr_id, download_dir, rawdata_dir): srr_id for srr_id in srr_ids}
        
        for future in as_completed(futures):
            srr_id = futures[future]
            try:
                result = future.result()
                if result:
                    success_count += 1
            except Exception as e:
                print(f"[Error] Sample {srr_id} encountered exception: {e}")

    print("="*40)
    print(f"[Summary] Sample: {len(srr_ids)}, Successfully processed: {success_count}")
    if success_count == len(srr_ids):
        print("[Step] Sample！")
        sys.exit(0)
    else:
        print("[Step] Samplefailed，。")
        sys.exit(1)

if __name__ == "__main__":
    main()
