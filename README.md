# Repeat-Aware_mutation_rate_estimator

This tool estimates the mutation rate between two sets of sequences, such as two reference or assembly genomes, or an assembly genome and raw sequencing reads. The method is designed to be robust to input sequences with high repetitiveness.

This repository builds upon the ideas of [this work](https://github.com/medvedevgroup/Repeat-Aware_Substitution_Rate_Estimator/), with substantial optimizations for scalability and speed. By leveraging KMC for k-mer counting and sourmash for sketching, our implementation enables fast mutation rate estimation between two ~3GB genomes in approximately 10 minutes of wall-clock time (with $k=31$ and `theta`=$1\times 10^{-2}$).

---

## Installation

We recommend using Conda to set up the environment:


Install all the packages based on `requirements.txt`

```bash
conda create -n <env_name> python=3.10 --file requirements.txt -c conda-forge -c bioconda -y
conda activate <env_name>
```

---

## Example Usage

We provide a small example under `example_data/`.

### Run from sequence mode (two FASTA files):

```bash
python Mutation_rate_estimator.py \
  --mode sequence \
  --input1 example_data/origin_seq.fasta \
  --input2 example_data/mutated_seq.fasta \
  --k 31 \
  --theta 0.01
```

### Run from mixture mode (one FASTA file + k-mer FASTA):

```bash
python Mutation_rate_estimator.py \
  --mode mixture \
  --input1 example_data/origin_seq.fasta \
  --input2 example_data/mutated_kmers.fasta \
  --k 31 \
  --theta 0.01
```

### Run from k-mer mode (two k-mer FASTA + histogram):

```bash
python Mutation_rate_estimator.py \
  --mode kmer \
  --input1 example_data/origin_kmers.fasta \
  --input2 example_data/mutated_kmers.fasta \
  --dist example_data/histogram.csv \
  --k 31 \
  --theta 0.01
```

---

## Usage

### Inputs

#### Sequence mode

Provide two complete sequences (could include multiple sequences) in separate FASTA files:

```bash
python Mutation_rate_estimator.py \
  --mode sequence \
  --input1 seq1.fasta \
  --input2 seq2.fasta \
  --k 31
```

#### Mixture mode

Provide sequences in FASTA format (could include multiple sequences) and one set of k-mers in FASTA format:

```bash
python Mutation_rate_estimator.py \
  --mode mixture \
  --input1 seq1.fasta \
  --input2 set2.fasta \
  --k 31
```

#### K-mer mode

Provide two sets of k-mers in FASTA format and a histogram CSV file. The histogram file should contain:

```csv
count,num_kmers
1,12345
2,6789
...
```

Then run:

```bash
python Mutation_rate_estimator.py \
  --mode kmer \
  --input1 set1.fasta \
  --input2 set2.fasta \
  --dist occ.csv \
  --k 31
```

### Using sketching

You can use the `--theta` parameter to apply FracMinHash sketching to approximate the intersection size between two k-mer sets.

- `theta` controls the sampling rate (default = 0.01). 
- The value must be in the range `(0, 1]`.
- For example, `--theta 0.01` retains approximately 1% of the hashed k-mers.

Example:

```bash
python Mutation_rate_estimator.py \
  --mode sequence \
  --input1 genome1.fasta \
  --input2 genome2.fasta \
  --k 31 \
  --theta 0.01
```



Note: Sketching could significantly accelerate large-scale computation. In this repository, we aim to estimate mutation rates between large-scale sequence datasets, such as reference genomes or whole-genome sequencing reads. To enable efficient computation, we use [sourmash](https://sourmash.readthedocs.io/en/latest/index.html) to implement FracMinHash sketching.

Note that sourmash requires a `--scaled` parameter of at least 100, which corresponds to a `theta` value of at most $0.01$ here. If you provide a larger `theta`, sourmash will automatically enforce it to be 0.01.

If you prefer exact intersection computation or are working with smaller datasets, we recommend using the original version [Repeat-Aware_Substitution_Rate_Estimator](https://github.com/medvedevgroup/Repeat-Aware_Substitution_Rate_Estimator). However, for large-scale data, we highly recommend using this sketch-based tool for its speed and scalability.



## Other Parameters

1. `--k`: The k-mer size, which must be an integer between 1 and 256 due to limitations of [KMC](https://github.com/refresh-bio/KMC). In practice, larger values of `k` are not useful, so this range is expected to be sufficient for most use cases.

2. `--cleanup`: If specified, all intermediate files generated during execution (such as KMC databases and sketch files) will be removed automatically after the program finishes.

3. `--use-dump`: By default, we call the C++ program `./cpp/kmc_histogram.cpp`, which uses the KMC API to read the kmc database and generate the abundance histogram directly in csv format—avoiding the overhead of dumping the k-mer count file to disk. If you specify `--use-dump`, the program will instead call `kmc_dump` to generate the full k-mer count file and then read it from disk. This may increase both runtime and disk space usage. If you want to retain the dumped k-mer count file, avoid using `--cleanup`.

