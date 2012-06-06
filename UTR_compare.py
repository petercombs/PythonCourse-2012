import pysam
import bx.bbi.bigwig_file as bigwig
import sys
import progressbar
from Utils import parse_gtf, find_nearest_genes, filter_gtf
import argparse

from Bio import SeqIO

GENOME_SIZE = 4639675  # E. coli K12 MG1655 from Ensembl
DOWNSTREAM = 50
DOWNSTREAM_KEEP = 100

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf-fname', '-G')
    parser.add_argument('--no-drug-fname', '-0')
    parser.add_argument('--drug-fname', '-d')
    parser.add_argument('--genome', '-g', type=open)
    parser.add_argument('--outfh', '-o', type=argparse.FileType('w'))

    parser.add_argument('--downstream_start', '-S', type=int, default=0)
    parser.add_argument('--downstream_end', '-E', type=int, default=100)
    parser.add_argument('--downstream_region', '-D', type=int, default=50)
    parser.add_argument('--ratio-high', '-R', type=float, default=1e100)
    parser.add_argument('--ratio_low', '-r', type=float, default=10.0)
    parser.add_argument('--keep-neg', dest='keep_pos', default=1,
                        action='store_const', const=-1)
    parser.add_argument('--auto-name-out', '-A', default=False,
                        action="store_true")
    args =  parser.parse_args()
    print args.keep_pos
    if args.auto_name_out:
        name = 'hidrug' if (args.keep_pos > 0) else 'drug_nonresp'
        name += '_utrs_'
        name += '%d_%d.fa' % (args.downstream_start, args.downstream_end)
        print "making file: ", name
        args.outfh = open(name, 'w')
    args.Genome = {r.name : r
              for r in SeqIO.parse(args.genome, 'fasta')}

    args.nonfh = False
    return args



if __name__ == "__main__":
    args = parse_args()

    gtf_data = parse_gtf(args.gtf_fname)
    nearest_to_left, nearest_to_right = find_nearest_genes(gtf_data,
                                                           GENOME_SIZE)

    no_drug_reads = bigwig.BigWigFile(open(args.no_drug_fname))
    drug_reads = bigwig.BigWigFile(open(args.drug_fname))

    filter_gtf(gtf_data)
    all_ratios = {}
    outputted_data = []
    pbar = progressbar.ProgressBar(widgets=['Pileups: ',
                                           progressbar.Percentage(), ' ',
                                           progressbar.Bar(), ' ',
                                           progressbar.ETA()],
                                   maxval=len(gtf_data))
    for gene in pbar(gtf_data):
        gene_data = gene.other.split(';')
        gene_data = {dat.split()[0] : dat.split()[1] for dat in gene_data if dat}
        gene_id = gene_data.get('gene_name', 'unknown_gene').strip('"')
        if gene.strand == '+':
            if gene.end + args.downstream_region > nearest_to_right[gene.end + 3]:
                continue
            drug_pileup = drug_reads.get_as_array(gene.chrom, gene.end,
                                            gene.end+args.downstream_region)
            no_drug_pileup = no_drug_reads.get_as_array(gene.chrom, gene.end,
                                                  gene.end+args.downstream_region)
            drug_counts = sum(drug_pileup)
            no_drug_counts = sum(no_drug_pileup)
            if no_drug_counts:
                all_ratios[gene.other] = drug_counts / no_drug_counts
                if (args.ratio_low <  all_ratios[gene.other] < args.ratio_high
                    and args.outfh):
                    seq = args.Genome[gene.chrom][gene.end + args.downstream_start :
                                             gene.end+args.downstream_end]
                    seq.id = gene_id
                    if gene_id not in outputted_data:
                        seq.description = '%d-%d' % (gene.start, gene.end)
                        SeqIO.write(seq, args.outfh, 'fasta')
                        outputted_data.append(gene_id)
        elif gene.strand == '-':
            if gene.start - args.downstream_region < nearest_to_left[gene.start - 3]:
                continue
            drug_pileup = drug_reads.get_as_array(gene.chrom,
                                            gene.start - args.downstream_region,
                                            gene.start)
            no_drug_pileup = no_drug_reads.get_as_array(gene.chrom,
                                                  gene.start - args.downstream_region,
                                                  gene.start)

            drug_counts = sum(drug_pileup)
            no_drug_counts = sum(no_drug_pileup)
            if no_drug_counts:
                all_ratios[gene.other] = drug_counts / no_drug_counts
                if (args.ratio_low <  all_ratios[gene.other] < args.ratio_high
                    and args.outfh):
                    seq = args.Genome[gene.chrom][gene.start - args.downstream_end :
                                             gene.start - args.downstream_start]
                    seq = seq.reverse_complement()
                    seq.id = gene_id
                    if gene_id not in outputted_data:
                        seq.description = '%d-%d' % (gene.start, gene.end)
                        SeqIO.write(seq, args.outfh, 'fasta')
                        outputted_data.append(gene_id)
                elif all_ratios[gene.other] < 1 and args.nonfh:
                    seq = args.Genome[gene.chrom][gene.start - args.downstream_end :
                                             gene.start - args.downstream_start]
                    seq = seq.reverse_complement()
                    seq.id = gene_id
                    seq.description = '%d-%d' % (gene.start, gene.end)
                    SeqIO.write(seq, nonfh, 'fasta')

    args.outfh.close()
