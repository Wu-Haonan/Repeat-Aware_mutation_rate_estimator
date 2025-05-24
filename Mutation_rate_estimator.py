# Mutation_rate_estimator.py
import argparse
import os
import csv
import subprocess
import tempfile
from collections import defaultdict
from sourmash import load_one_signature, load_file_as_signatures
from scipy.optimize import newton,brentq

middle_dir = "./middle_results/"
os.makedirs(middle_dir, exist_ok=True)

# ---------- KMC wrapper ----------
def run_kmc(input_fasta, db_prefix, k):
    db_path = os.path.join(middle_dir, db_prefix)
    tmp_dir = os.path.join(middle_dir, "tmp_kmc")
    os.makedirs(tmp_dir, exist_ok=True)
    
    cmd = f"kmc -k{k} -ci1 -t8 -fm {input_fasta} {db_path} {tmp_dir}"
    print(f"[INFO] Running: {cmd}")
    subprocess.run("ulimit -n 2048", shell=True, check=True)
    subprocess.run(cmd, shell=True, check=True)
    
    return db_path

def build_histogram_cpp(db_path):
    """avoid dump"""
    
    cpp_executable = "./cpp/kmc_histogram"
    
    if not os.path.exists(cpp_executable):
        print(f"[ERROR] C++ histogram tool not found: {cpp_executable}")
        print("[INFO] Please compile it first:")
        print("  cd cpp && make")
        return None, None
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_file:
        csv_path = tmp_file.name
    
    try:
        cmd = [cpp_executable, db_path, csv_path]
        print(f"[INFO] Running C++ histogram: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if result.stdout:
            print(result.stdout)
        
        hist = {}
        total_kmers = 0
        
        with open(csv_path, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # jump head
            
            for row in reader:
                if len(row) >= 2:
                    count = int(row[0])
                    num_kmers = int(row[1])
                    hist[count] = num_kmers
                    total_kmers += num_kmers
        
        print(f"[INFO] C++ histogram completed: {total_kmers:,} total k-mers")
        return hist, total_kmers
        
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] C++ histogram failed: {e.stderr}")
        return None, None
    
    except Exception as e:
        print(f"[ERROR] Failed to process histogram: {e}")
        return None, None
    
    finally:
        if os.path.exists(csv_path):
            os.unlink(csv_path)

def dump_kmc_to_txt(db_prefix, output_file="kmc_dump.txt"):
    """alternative, dump"""
    output_path = os.path.join(middle_dir, output_file)
    cmd = f"kmc_dump {db_prefix} {output_path}"
    print(f"[INFO] Running fallback dump: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    
    return output_path

# ---------- Histogram ----------
def build_histogram_from_kmc_dump(dump_file):
    """alternative, read dump"""
    hist = defaultdict(int)
    total_kmers = 0
    with open(dump_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            count = int(parts[1])
            total_kmers += 1
            hist[count] += 1
    return hist, total_kmers

def read_histogram_csv(csv_file):
    hist = {}
    total_kmers = 0
    with open(csv_file) as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) != 2:
                continue
            i, n = map(int, row)
            total_kmers += n
            hist[i] = n
    return hist, total_kmers

# ---------- Sourmash ----------
def sketch(input_fasta, k, scaled, out_sig):
    out_sig_path = os.path.join(middle_dir, out_sig)
    cmd = f"sourmash sketch dna -p k={k},scaled={scaled},seed=42 {input_fasta} -o {out_sig_path}"
    print(f"[INFO] Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    
    return out_sig_path 

def compute_intersection(sig1_file, sig2_file):
    try:
        s1 = next(load_file_as_signatures(sig1_file))
        s2 = next(load_file_as_signatures(sig2_file))
        intersection_size = s1.minhash.count_common(s2.minhash)
        print(f"[INFO] Sketch intersection : {intersection_size}")
        return intersection_size
    except StopIteration:
        print("[ERROR] One of the files contains no signatures")
        return None
    except Exception as e:
        print(f"[ERROR] {e}")
        return None

# ---------- Newton Solver ----------
def solve_histogram_equation(hist, I, total_kmers):
    if I ==0:
        return 1
    else:
        def poly(x):
            return (total_kmers - I) - sum(count * x**i for i, count in hist.items()) 
        try:
            #root = newton(poly, x0=0.3, tol=1e-15)
            root = brentq(poly, 0, 1, xtol=1e-15)
        except Exception as e:
            print(f"[ERROR] Newton solver failed: {e}")
            return None
        return root

# ---------- Cleanup function ----------
def cleanup_middle_files():
    import shutil
    if os.path.exists(middle_dir):
        shutil.rmtree(middle_dir)

# ---------- Main pipeline ----------
def main():
    parser = argparse.ArgumentParser(description="Estimate mutation rate r between two sequences or k-mer sets")
    parser.add_argument("--mode", required=True, choices=["sequence", "mixture", "kmer"], help="Input mode")
    parser.add_argument("--input1", required=True, help="First input file (FASTA or k-mer set)")
    parser.add_argument("--input2", required=True, help="Second input file (FASTA or k-mer set)")
    parser.add_argument("--dist", help="Optional: histogram CSV file (only for kmer mode)")
    parser.add_argument("--k", type=int, default=31, help="k-mer size")
    parser.add_argument("--theta", type=float, default=1, help="Sketching threshold (default=1)")
    parser.add_argument("--cleanup", action="store_true", help="Clean up middle files after completion")
    parser.add_argument("--use-dump", action="store_true", help="Force use traditional dump method")
    args = parser.parse_args()

    print(f"[INFO] Middle files will be stored in: {middle_dir}")
    scaled = int(1 / args.theta)

    if args.mode in ["sequence", "mixture"]:
        db_path = run_kmc(args.input1, "kmc_db", args.k)
        
        # call kmc_histgram.cpp
        if not args.use_dump:
            hist, total_kmers = build_histogram_cpp(db_path)
        else:
            hist, total_kmers = None, None
        
        # if c++ fail
        if hist is None:
            print("[INFO] Falling back to traditional dump method...")
            dump_file_path = dump_kmc_to_txt(db_path, "kmc_dump.txt")
            hist, total_kmers = build_histogram_from_kmc_dump(dump_file_path)
        
    elif args.mode == "kmer":
        if not args.dist:
            raise ValueError("--dist must be provided in kmer mode")
        hist, total_kmers = read_histogram_csv(args.dist)

    sig1_path = sketch(args.input1, args.k, scaled, "sig1.sig")
    sig2_path = sketch(args.input2, args.k, scaled, "sig2.sig")
    
    I = compute_intersection(sig1_path, sig2_path) / args.theta
    
    print(f"[INFO] Total number of kmers of string 1 L0: {total_kmers}")
    print(f"[INFO] Estimated intersection I: {I}")
    
    if (I >= total_kmers):
        r = 0
    else:
        q = solve_histogram_equation(hist, I, total_kmers)
        r = 1 - (1 - q) ** (1/args.k)
    
    print(f"[RESULT] Estimated r : {r}")
    
    # option: clean middle data
    if args.cleanup:
        cleanup_middle_files()

if __name__ == "__main__":
    main()