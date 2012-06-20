""" CalcFalloff -- Calculates the read falloff around a gene

This is code for finding the average(?) expression across all genes as a
function of position. We'd like to see the upstream area (which should be
relatively low), the gene itself (relatively high, and nearly constant across
the gene), and then the downstream, which ought to fall off quickly in the
no-drug case, and fall off much more slowly in the +drug case.
"""

from __future__ import division
from argparse import ArgumentParser
import numpy as np
import progressbar
import sys
from Utils import parse_gtf, find_nearest_genes, filter_gtf
import bx.bbi.bigwig_file as bigwig


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-d', '--downstream', type=int, default=300)
    parser.add_argument('-u', '--upstream', type=int, default=300)
    parser.add_argument('-r', '--resolution', type=int, default=200)
    parser.add_argument('-c', '--clearance', type=int, default=50)
    parser.add_argument('-g', '--gtf-filename',
                   default='E_coli_k12.EB1_e_coli_k12.13.gtf')
    parser.add_argument('-S', '--genome-size', type=int, default=4639675)
    parser.add_argument('bigwig', type=FileType('rb'))

    global args
    args = parser.parse_args()
    return args

def find_upper_limits(end, strand, to_left, to_right):
    """ """
    if strand == '+':
        dist = args.downstream
    else:
        dist = args.upstream
    no_gene_pos = end + dist
    gene_pos = to_right[end] - args.clearance
    downstream = max(min(no_gene_pos, gene_pos),
                     end)
    return downstream

def find_lower_limits(start, strand, to_left, to_right):
    """ """
    if strand == "+":
        dist = args.upstream
    else:
        dist = args.downstream

    no_gene_pos = start - dist
    gene_pos = nearest_to_left[start] + args.clearance
    upstream= min(max(no_gene_pos, gene_pos),
                         start)
    return upstream

def calc_deltas(us, ue, ds, de, strand):
    delta_up = np.zeros(args.upstream)
    delta_down = np.zeros(args.downstream)

    if strand == '+':
        delta_up[args.upstream - ue + us:] += 1
        delta_down[:de - ds] += 1
    elif strand == '-':
        delta_up[args.upstream - de + ds:] += 1
        delta_down[:ue - us] += 1

    return delta_up, delta_down

if __name__ == "__main__":
    args = parse_args()
    gtf_data = parse_gtf(args.gtf_filename)
    nearest_to_left, nearest_to_right = find_nearest_genes(gtf_data,
                                                           args.genome_size)

    #reads = pysam.Samfile(sys.argv[1], 'rb')
    reads = bigwig.BigWigFile(args.bigwig)

    upstream = np.zeros(args.upstream) # number of reads (calculated with pileup)
    upstream_n = np.zeros(args.upstream) # number of genes used at each position
    downstream = np.zeros(args.downstream)
    downstream_n = np.zeros(args.downstream)

    gene_cov = np.zeros(args.resolution)
    gene_cov_n = 0

    dists = []

    filter_gtf(gtf_data)
    pbar = progressbar.ProgressBar(widgets=['Pileups: ',
                                           progressbar.Percentage(), ' ',
                                           progressbar.Bar(), ' ',
                                           progressbar.ETA()],
                                   maxval=len(gtf_data))
    all_ratios = {}
    for gene in pbar(gtf_data):
        start = gene.start
        end = gene.end
        strand = gene.strand
        chrom = gene.chrom # Should always be "Chromosome"

        if start == 4275492: continue
        # soxR-1 has sraL with an overlapping region with sraL.  I'm a little
        # surprised this seems to be the only problem...

        # Look to the left of the gene
        # (upstream for the + strand, downstream for the - strand)

        upstream_start = find_lower_limits(start, strand, nearest_to_left,
                                           nearest_to_right)
        upstream_end = start
        downstream_start = end
        downstream_end = find_upper_limits(end, strand, nearest_to_left,
                                           nearest_to_right)

        delta_up, delta_down = calc_deltas(upstream_start, upstream_end,
                                           downstream_start, downstream_end,
                                           strand)

        dists.append(upstream_end - upstream_start)
        total_downstream = 0.0


        #for col in reads.pileup(chrom, upstream_start, upstream_end):
        if upstream_start < upstream_end:
            for n, pos in zip(reads.get_as_array(chrom, upstream_start,
                                                   upstream_end),
                                range(upstream_start, upstream_end)):
                if not upstream_start <= pos < upstream_end:
                    continue

                if strand == '+':
                    i = args.upstream - start + pos
                    upstream[i] += n
                elif strand == '-':
                    i = start - pos - 1
                    total_downstream += n
                    downstream[i] += n
                    #if col.n > 50 and i > 50:
                        #print "-Something fishy...", start, col.pos, col.n


        # Look at the gene itself
        last_pct = -1
        total_reads = 0
        for n, pos in zip(reads.get_as_array(chrom, start, end),
                          range(start, end)):
            if not start <= pos < end:
                continue
            if last_pct == -1:
                # Genes with any expression whatsoever should get included
                gene_cov_n += 1
                upstream_n += delta_up
                downstream_n += delta_down
                last_pct = 0
            # Note the from future import __division__
            if strand == '+':
                hi = int(np.floor(args.resolution * (pos - start)/(end - start)))
            elif strand == '-':
                hi = int(np.floor(args.resolution * (end - pos)/(end - start)))

            total_reads += n
            if last_pct == hi:
                gene_cov[hi] += n
            else:
                gene_cov[last_pct:hi] += n
            last_pct = hi

        mean_reads = total_reads / float(end - start)


        # Look to the right of the gene
        # (downstream for the + strand, upstream for the - strand)

        dists.append(downstream_end - downstream_start)

        if downstream_start < downstream_end:
            for n, pos in zip(reads.get_as_array(chrom, downstream_start,
                                                 downstream_end) ,
                              range(downstream_start, downstream_end)):
                if not downstream_start < pos < downstream_end:
                    continue

                if strand == '+':
                    i = pos - end - 1
                    downstream[i] += n
                    total_downstream += n
                    #if n > 50 and i > 50:
                        #print "+Something fishy...", start, pos, n
                elif strand == '-':
                    i = args.upstream - pos + end
                    upstream[i] += n

        if sum(delta_down):
            all_ratios[gene.other] = total_downstream / sum(delta_down) / total_reads



