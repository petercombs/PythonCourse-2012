from __future__ import division
import pysam
from collections import namedtuple
import numpy as np

def parse_gtf(filename):
    Gene = namedtuple("Gene", ['chrom', 'start', 'end', 'strand', 'other'])
    genelist = []
    for line in open(filename):
        if line.startswith('#') or line.startswith('F'):
            continue

        data = line.split('\t')
        if data[2] != 'CDS':
            continue
        genelist.append(Gene(chrom=data[0],
                             start=int(data[3]), end=int(data[4]),
                             strand=data[6], other=data[-1]))
    return genelist

def find_nearest_genes(gtf_data, genome_size):
    nearest_to_left = [0]*genome_size
    nearest_to_right = [0]*genome_size

    old_start = 0
    old_end = 0
    first_start = 0
    first_end = 0
    gene = None
    for gene in gtf_data:
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

if __name__ == "__main__":
    genome_size = 4639675  # E. coli K12 MG1655 from Ensembl
    RESOLUTION = 200
    UPSTREAM = 100
    DOWNSTREAM = 100
    CLEARANCE = 200
    gtf_filename = 'E_coli_k12.EB1_e_coli_k12.13.gtf'

    gtf_data = parse_gtf(gtf_filename)
    nearest_to_left, nearest_to_right = find_nearest_genes(gtf_data,
                                                           genome_size)

    f = pysam.Samfile('tophat_out/accepted_hits.bam', 'rb')

    upstream = np.zeros(UPSTREAM) # number of reads (calculated with pileup)
    upstream_n = np.zeros(UPSTREAM) # number of genes used at each position
    downstream = np.zeros(DOWNSTREAM)
    downstream_n = np.zeros(DOWNSTREAM)

    gene_cov = np.zeros(RESOLUTION)
    gene_cov_n = len(gtf_data)
    # Every gene will get covered completely, so n would just be the length of
    # gtf_data

    for gene in gtf_data:
        start = gene.start
        end = gene.end
        strand = gene.strand
        chrom = gene.chrom # Should always be "Chromosome"

        # Look to the left of the gene 
        # (upstream for the + strand, downstream for the - strand)

        dist = UPSTREAM if strand == '+' else DOWNSTREAM
        upstream_start = min(max(nearest_to_left[start] + CLEARANCE,
                                 start - dist),
                             start)
        upstream_end = start

        if strand == '+':
            upstream_n[UPSTREAM - start + upstream_start:UPSTREAM] += 1
        elif strand == '-':
            downstream_n[0:start - upstream_start] += 1

        for col in f.pileup(chrom, upstream_start, upstream_end):
            if not upstream_start <= col.pos < upstream_end:
                continue

            if strand == '+':
                i = UPSTREAM - start + col.pos
                upstream[i] += col.n
            elif strand == '-':
                i = start - col.pos - 1
                downstream[i] += col.n


        # Look at the gene itself
        last_pct = 0
        for col in f.pileup(chrom, start, end):
            if not start <= col.pos < end:
                continue
            # Note the from future import __division__
            if strand == '+':
                hi = int(np.floor(RESOLUTION * (col.pos - start)/(end - start)))
            elif strand == '-':
                hi = int(np.floor(RESOLUTION * (end - col.pos)/(end - start)))

            gene_cov[last_pct:hi] += col.n
            last_pct = hi



        # Look to the right of the gene
        # (downstream for the + strand, upstream for the - strand)

        dist = DOWNSTREAM if strand == '+' else UPSTREAM

        downstream_start = end
        downstream_end = max(min(nearest_to_right[end] - CLEARANCE,
                                 end + dist),
                             end)

        if strand == '+':
            downstream_n[0:downstream_end - end] += 1
        elif strand == '-':
            upstream_n[UPSTREAM - downstream_end + end: UPSTREAM]

        for col in f.pileup(chrom, downstream_start, downstream_end):
            if not downstream_start < col.pos < downstream_end:
                continue

            if strand == '+':
                i = col.pos - end - 1
                downstream[i] += col.n
            elif strand == '-':
                i = UPSTREAM - col.pos + end
                upstream[i] += col.n




