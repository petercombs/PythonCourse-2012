from Bio import SeqIO
import sys

is_one = sys.argv[2]


all_seqs = []

for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    seq = [(1 if base in is_one else 0) for base in rec.seq]
    all_seqs.append(seq)
