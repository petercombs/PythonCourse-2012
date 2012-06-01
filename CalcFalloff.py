""" CalcFalloff -- Calculates the read falloff around a gene

This is code for finding the average(?) expression across all genes as a
function of position. We'd like to see the upstream area (which should be
relatively low), the gene itself (relatively high, and nearly constant across
the gene), and then the downstream, which ought to fall off quickly in the
no-drug case, and fall off much more slowly in the +drug case.
"""

from __future__ import division
import pysam
import numpy as np
import progressbar
import sys
from Utils import parse_gtf, find_nearest_genes, filter_gtf


# Define constants
GENOME_SIZE = 4639675  # E. coli K12 MG1655 from Ensembl
RESOLUTION = 200
UPSTREAM = 300
DOWNSTREAM = 300
CLEARANCE = 50
gtf_filename = 'E_coli_k12.EB1_e_coli_k12.13.gtf'

def find_upper_limits(end, strand, to_left, to_right):
    """ """
    if strand == '+':
        dist = DOWNSTREAM 
    else:
        dist = UPSTREAM
    no_gene_pos = end + dist
    gene_pos = to_right[end] - CLEARANCE
    downstream = max(min(no_gene_pos, gene_pos),
                     end)
    return downstream

def find_lower_limits(start, strand, to_left, to_right):
    """ """
    if strand == "+":
        dist = UPSTREAM
    else:
        dist = DOWNSTREAM

    no_gene_pos = start - dist
    gene_pos = nearest_to_left[start] + CLEARANCE
    upstream= min(max(no_gene_pos, gene_pos),
                         start)
    return upstream

def calc_deltas(us, ue, ds, de, strand):
    delta_up = np.zeros(UPSTREAM)
    delta_down = np.zeros(DOWNSTREAM)

    if strand == '+':
        delta_up[UPSTREAM - ue + us:] += 1
        delta_down[:de - ds] += 1
    elif strand == '-':
        delta_up[UPSTREAM - de + ds:] += 1
        delta_down[:ue - us] += 1

    return delta_up, delta_down

if __name__ == "__main__":
    gtf_data = parse_gtf(gtf_filename)
    nearest_to_left, nearest_to_right = find_nearest_genes(gtf_data,
                                                           GENOME_SIZE)

    reads = pysam.Samfile(sys.argv[1], 'rb')

    upstream = np.zeros(UPSTREAM) # number of reads (calculated with pileup)
    upstream_n = np.zeros(UPSTREAM) # number of genes used at each position
    downstream = np.zeros(DOWNSTREAM)
    downstream_n = np.zeros(DOWNSTREAM)

    gene_cov = np.zeros(RESOLUTION)
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


        for col in reads.pileup(chrom, upstream_start, upstream_end):
            if not upstream_start <= col.pos < upstream_end:
                continue

            if strand == '+':
                i = UPSTREAM - start + col.pos
                upstream[i] += col.n
            elif strand == '-':
                i = start - col.pos - 1
                total_downstream += col.n
                downstream[i] += col.n
                #if col.n > 50 and i > 50:
                    #print "-Something fishy...", start, col.pos, col.n


        # Look at the gene itself
        last_pct = -1
        total_reads = 0
        for col in reads.pileup(chrom, start, end):
            if not start <= col.pos < end:
                continue
            if last_pct == -1:
                # Genes with any expression whatsoever should get included
                gene_cov_n += 1
                upstream_n += delta_up
                downstream_n += delta_down
                last_pct = 0
            # Note the from future import __division__
            if strand == '+':
                hi = int(np.floor(RESOLUTION * (col.pos - start)/(end - start)))
            elif strand == '-':
                hi = int(np.floor(RESOLUTION * (end - col.pos)/(end - start)))

            total_reads += col.n
            if last_pct == hi:
                gene_cov[hi] += col.n
            else:
                gene_cov[last_pct:hi] += col.n
            last_pct = hi

        mean_reads = total_reads / float(end - start)


        # Look to the right of the gene
        # (downstream for the + strand, upstream for the - strand)

        dists.append(downstream_end - downstream_start)

        for col in reads.pileup(chrom, downstream_start, downstream_end):
            if not downstream_start < col.pos < downstream_end:
                continue

            if strand == '+':
                i = col.pos - end - 1
                downstream[i] += col.n
                total_downstream += col.n
                #if col.n > 50 and i > 50:
                    #print "+Something fishy...", start, col.pos, col.n
            elif strand == '-':
                i = UPSTREAM - col.pos + end
                upstream[i] += col.n

        if sum(delta_down):
            all_ratios[gene.other] = total_downstream / sum(delta_down) / total_reads



