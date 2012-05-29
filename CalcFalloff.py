""" CalcFalloff -- Calculates the read falloff around a gene

This is code for finding the average(?) expression across all genes as a
function of position. We'd like to see the upstream area (which should be
relatively low), the gene itself (relatively high, and nearly constant across
the gene), and then the downstream, which ought to fall off quickly in the
no-drug case, and fall off much more slowly in the +drug case.
"""

from __future__ import division
import pysam
from collections import namedtuple
import numpy as np
import progressbar


# Define constants
Gene = namedtuple("Gene", ['chrom', 'start', 'end', 'strand', 'type', 'other'])
GENOME_SIZE = 4639675  # E. coli K12 MG1655 from Ensembl
RESOLUTION = 200
UPSTREAM = 100
DOWNSTREAM = 100
CLEARANCE = 200
gtf_filename = 'E_coli_k12.EB1_e_coli_k12.13.gtf'

def parse_gtf(filename):
    """ Parses the GTF in filename into a list of Genes

    At the moment, this strips out everything that's a comment or on the F
    chromosome, as well as anything that has the exact same extent as the
    previous entry.
    """
    genelist = []
    last_pos = ('', 0, 0)
    for line in open(filename):
        if line.startswith('#') or line.startswith('F'):
            continue

        data = line.split('\t')
        if (data[0], data[3], data[4]) != last_pos:
            genelist.append(Gene(chrom=data[0],
                                 start=int(data[3]),
                                 end=int(data[4]),
                                 strand=data[6],
                                 type=data[2],
                                 other=data[-1].strip()))
            last_pos = data[0], data[3], data[4]
    return genelist

def find_nearest_genes(gtf_data, genome_size):
    """For each base on the genome, find the nearest genes to the left and right

    This is useful for doing lookups, although it may be too large to scale to
    non-bacterial organisms.
    """
    nearest_to_left = [0]*genome_size
    nearest_to_right = [0]*genome_size

    old_start = 0
    old_end = 0
    first_start = 0
    first_end = 0
    gene = None

    pbar = progressbar.ProgressBar(maxval=len(gtf_data))

    for gene in pbar(gtf_data):
        if gene.chrom != 'Chromosome': continue
        for i in range(old_start, gene.start):
            nearest_to_right[i] = gene.start
        if first_end == 0:
            first_start = gene.start
            first_end = gene.end
        else:
            for i in range(old_end, gene.end):
                nearest_to_left[i] = old_end

        old_start = gene.start
        old_end = gene.end

    assert gene  # Should catch an empty gtf file

    for i in range(first_end):
        nearest_to_left[i] = gene.end
    for i in range(gene.end, len(nearest_to_left)):
        nearest_to_left[i] = gene.end
    for i in range(gene.start, len(nearest_to_right)):
        nearest_to_right[i] = first_start



    return nearest_to_left, nearest_to_right

def filter_gtf(gtf_data):
    """ Remove non-CDS entries from the GTF data"""
    for entry in gtf_data[:]:
        if entry.type != 'CDS':
            gtf_data.remove(entry)

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

    reads = pysam.Samfile('SRR400620.bam', 'rb')

    upstream = np.zeros(UPSTREAM) # number of reads (calculated with pileup)
    upstream_n = np.zeros(UPSTREAM) # number of genes used at each position
    downstream = np.zeros(DOWNSTREAM)
    downstream_n = np.zeros(DOWNSTREAM)

    gene_cov = np.zeros(RESOLUTION)
    gene_cov_n = 0

    dists = []

    filter_gtf(gtf_data)
    pbar = progressbar.ProgressBar(maxval=len(gtf_data))
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


        for col in reads.pileup(chrom, upstream_start, upstream_end):
            if not upstream_start <= col.pos < upstream_end:
                continue

            if strand == '+':
                i = UPSTREAM - start + col.pos
                upstream[i] += col.n
            elif strand == '-':
                i = start - col.pos - 1
                downstream[i] += col.n
                if col.n > 50 and i > 50:
                    print "-Something fishy...", start, col.pos, col.n


        # Look at the gene itself
        last_pct = -1
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

            if last_pct == hi:
                gene_cov[hi] += col.n
            else:
                gene_cov[last_pct:hi] += col.n
            last_pct = hi



        # Look to the right of the gene
        # (downstream for the + strand, upstream for the - strand)

        dists.append(downstream_end - downstream_start)

        for col in reads.pileup(chrom, downstream_start, downstream_end):
            if not downstream_start < col.pos < downstream_end:
                continue

            if strand == '+':
                i = col.pos - end - 1
                downstream[i] += col.n
                if col.n > 50 and i > 50:
                    print "+Something fishy...", start, col.pos, col.n
            elif strand == '-':
                i = UPSTREAM - col.pos + end
                upstream[i] += col.n




