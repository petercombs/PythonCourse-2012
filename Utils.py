from collections import namedtuple
import progressbar

Gene = namedtuple("Gene", ['chrom', 'start', 'end', 'strand', 'type', 'other'])

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
        if entry.type != 'CDS' and entry.type != 'operon':
            gtf_data.remove(entry)


