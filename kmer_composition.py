#!/usr/bin/env python
import argparse

from Bio import SeqIO
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

forward = "ACGT"
reverse = "TGCA"
trans_table = str.maketrans(forward, reverse)


def kmer_counter(args):
    """Must be a top-level function to be picklable by multiprocessing."""
    seq_id, seq, k = args
    forward = "ACGT"
    reverse = "TGCA"
    trans_table = str.maketrans(forward, reverse)

    l = len(seq)
    if l < k:
        return seq_id, {}, set()

    counter = {}
    local_kmers = set()

    for i in range(l - k + 1):
        kmer = seq[i:(i + k)]
        if "N" in kmer:
            continue
        kmer_rev = kmer.translate(trans_table)[::-1]
        if kmer > kmer_rev:
            kmer = kmer_rev
        local_kmers.add(kmer)
        counter[kmer] = counter.get(kmer, 0) + 1

    return seq_id, counter, local_kmers


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate normalized canonical k-mer composition per sequence from a FASTA file."
    )
    parser.add_argument("-i", "--input", required=True, help="Input assembly FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file path")
    parser.add_argument("-k", "--kmer-size", dest="kmer_size", type=int, required=True, help="k-mer length")
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of worker processes (default: 4)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    fasta_file = args.input
    kmer_length = args.kmer_size
    output_file = args.output
    n_threads = args.threads

    print("Input file:", fasta_file)
    print("Workers:", n_threads)

    print("Loading sequences...")
    records = [(str(r.id), str(r.seq).upper()) for r in SeqIO.parse(fasta_file, "fasta")]
    n_sequences = len(records)
    print(n_sequences, "sequences")
    print("Counting kmers of length", kmer_length)

    # Build args list for workers
    args = [(seq_id, seq, kmer_length) for seq_id, seq in records]

    kmer_counts = {}
    all_kmers = set()

    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        futures = {executor.submit(kmer_counter, arg): arg[0] for arg in args}
        for future in tqdm(as_completed(futures), total=n_sequences):
            seq_id, counter, local_kmers = future.result()
            if counter:
                kmer_counts[seq_id] = counter
                all_kmers.update(local_kmers)

    all_kmers_sorted = sorted(all_kmers)
    print("Writing results to:", output_file)

    seq_ids = []
    X = []

    with open(output_file, "w") as fo:
        fo.write("sequence_id\t" + "\t".join(all_kmers_sorted) + "\n")
        for seq_id, _ in records:  # preserve original order
            if seq_id not in kmer_counts:
                continue
            seq_ids.append(seq_id)
            total_kmers = sum(kmer_counts[seq_id].values())
            line = [
                kmer_counts[seq_id][kmer] / total_kmers if kmer in kmer_counts[seq_id] else 0
                for kmer in all_kmers_sorted
            ]
            fo.write(seq_id + "\t" + "\t".join(map(str, line)) + "\n")
            X.append(line)

if __name__ == "__main__":
    main()
