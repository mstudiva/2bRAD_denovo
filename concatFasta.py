#!/usr/bin/env python3
"""
concatFasta_length_balanced.py

Convert multiple FASTA files (each representing a genome) into fake chromosomes.
Each input FASTA becomes K "chromosomes" by concatenating its contigs,
balancing lengths across the chromosomes.

Outputs:
  1) FASTA file of fake chromosomes
  2) TSV summary file listing stats per fake chromosome
  3) TSV mapping file of contig -> fake chromosome

Usage:
  python concatFasta_length_balanced.py \
    -o fakechrs.fa \
    -s summary.tsv \
    -m contig_map.tsv \
    genome1.fa genome2.fa ...
"""

import sys
import argparse
import os

def fasta_iter(path):
    """Yield (header, seq) from a FASTA file."""
    header = None
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].split()[0]  # take only first token as name
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        if header is not None:
            yield header, "".join(seq_parts)

def assign_by_length(contigs, k):
    """
    Assign contigs to k bins to balance total lengths.
    Greedy algorithm: largest contigs placed in the currently shortest bin.
    """
    bins = [[] for _ in range(k)]
    lengths = [0] * k
    sorted_contigs = sorted(contigs, key=lambda x: len(x[1]), reverse=True)
    for contig in sorted_contigs:
        i = min(range(k), key=lambda j: lengths[j])
        bins[i].append(contig)
        lengths[i] += len(contig[1])
    return bins

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fastas", nargs="+", help="Input FASTA files")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA")
    ap.add_argument("-s", "--summary", required=True, help="Output summary TSV")
    ap.add_argument("-m", "--mapping", required=True, help="Output contig-to-fakechrom mapping TSV")
    ap.add_argument("-k", "--k", type=int, default=5, help="Number of fake chromosomes per genome (default=5)")
    ap.add_argument("-g", "--gap", type=int, default=100, help="Number of Ns between contigs (default=100)")
    args = ap.parse_args()

    gap_seq = "N" * args.gap

    with open(args.output, "w") as out_fa, \
         open(args.summary, "w") as out_sum, \
         open(args.mapping, "w") as out_map:

        out_sum.write("genome_id\tfake_chr\tnum_contigs\ttotal_length\n")
        out_map.write("genome_id\tcontig\tfake_chr\n")

        for fasta in args.fastas:
            genome_id = os.path.splitext(os.path.basename(fasta))[0]
            contigs = list(fasta_iter(fasta))
            if not contigs:
                continue

            # distribute contigs across k bins, balancing by total length
            bins = assign_by_length(contigs, args.k)

            for i, bin_contigs in enumerate(bins, start=1):
                if not bin_contigs:
                    continue
                fakechr = f"{genome_id}_fakechr{i}"
                joined = gap_seq.join(seq for _, seq in bin_contigs)
                total_len = sum(len(seq) for _, seq in bin_contigs)

                # write FASTA
                header = f">{fakechr} original_contigs={len(bin_contigs)} total_length={total_len}"
                out_fa.write(header + "\n")
                for j in range(0, len(joined), 60):
                    out_fa.write(joined[j:j+60] + "\n")

                # write summary
                out_sum.write(f"{genome_id}\t{fakechr}\t{len(bin_contigs)}\t{total_len}\n")

                # write mapping
                for cname, _ in bin_contigs:
                    out_map.write(f"{genome_id}\t{cname}\t{fakechr}\n")

if __name__ == "__main__":
    main()
