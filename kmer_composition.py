#!/usr/bin/env python

import sys

if len(sys.argv) < 4:
    print("Usage: python kmer_composition.py [assembly.fasta] [kmer_length] [output_file]")
    sys.exit(0)

from Bio import SeqIO
from tqdm import tqdm

fasta_file = sys.argv[1]
kmer_length = int(sys.argv[2])
output_file = sys.argv[3]

forward = "ACGT"
reverse = "TGCA"
trans_table = str.maketrans(forward, reverse)

# dict to store kmer counts per sequence
kmer_counts = {}

# set to store all observed kmers in any sequence
all_kmers = set()

def main():
    print("Input file:", fasta_file)
    
    # count sequences for progress bar
    n_sequences = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        n_sequences += 1

    print(n_sequences, "sequences")
    print("Counting kmers of length", kmer_length)

    progress = tqdm(total = n_sequences)

    for record in SeqIO.parse(fasta_file, "fasta"):
        kmer_counts[str(record.id)] = kmer_counter(str(record.seq).upper(), kmer_length)
        progress.update(1)

    progress.close()

    all_kmers_sorted = sorted(all_kmers)

    print("Writing results to:", output_file)

    fo = open(output_file, "w")

    line = "sequence_id\t" + "\t".join(all_kmers_sorted)
    fo.write(line + "\n")

    # print frequency table
    for seq in kmer_counts:
        line = []
        
        # count total number of kmers for sequence
        total_kmers = sum(kmer_counts[seq].values())

        for kmer in all_kmers_sorted:
            if kmer in kmer_counts[seq]:
                line.append(str(kmer_counts[seq][kmer] / total_kmers))
            else:
                line.append("0")
       
        fo.write(seq + "\t" +  "\t".join(line) + "\n")
    
def kmer_counter(seq, k):
    l = len(seq)
    counter = {}

    if l < k:
        return

    for i in range(l - k + 1):
        kmer = seq[i:(i+k)]

        if "N" in kmer:
            continue
        
        kmer_rev = kmer.translate(trans_table)[::-1]

        if kmer > kmer_rev:
            kmer = kmer_rev
        
        all_kmers.add(kmer)

        if kmer in counter:
            counter[kmer] += 1
        else:
            counter[kmer] = 1

    return counter

if __name__ == "__main__":
    main()
