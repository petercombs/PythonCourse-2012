from Bio import SeqIO
import sys

is_one = sys.argv[2:]


def find_all(string, substring):
    loc = string.find(substring)
    locs = []
    while loc != -1:
        locs.append(loc)
        loc = string.find(substring, loc + 1)
    return locs


def find_all_from_list(string, list):
    all_locs = []
    for substring in list:
        all_locs.extend(find_all(string, substring))
    return all_locs


all_seqs = []

for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    hits = find_all_from_list(rec.seq, is_one)
    seq = [(1 if i in hits else 0) for i in range(len(rec))]
    all_seqs.append(seq)
