### NOTE: Rather than trying to go line by line through the operon calls, why
###  not try going line by line through the genes, and look up whether adjacent genes
###  are in the same operon if necessary

def get_gene_info(gene_info):
    """ Make a dictionary with start and stop information
    """
    start_stop_strand = {}

    #Skip the header line
    gene_info.readline()

    for line in gene_info:
        data = line.split('\t')
        start = int(data[4])
        stop = int(data[5])
        strand = data[6]
        sysName = data[7]
        name = data[8]
        start_stop_strand[name] = (start, stop, strand)
        start_stop_strand[sysName] = (start, stop, strand)
        if name == 'rrmJ':
            print "Wacky hijinks!"
            start_stop_strand['rlmE'] = (start, stop, strand)

    return start_stop_strand

def make_operon_name(gene_list):
    return '-'.join(gene_list)

def operon_calls_to_gtf(operon_calls_file, sss_dict):
    out_file = open('operons.gff', 'w')

    # Skip the header line
    operon_calls_file.readline()
    current_operon_names = []
    current_operon_start = 0
    current_operon_stop = 0
    last_gene_name = ''
    last_gene_start = 0
    last_gene_stop = 0
    last_gene_strand = ''
    for line in operon_calls_file:
        data = line.split('\t')
        name1 = data[4]
        name2 = data[5]
        start1, stop1, strand1 = sss_dict[name1]
        start2, stop2, strand2 = sss_dict[name2]
        is_operon = (data[6] == "TRUE")


        if name1 != last_gene_name and last_gene_name:
            annot = 'operonName "' + last_gene_name + '";'
            write_gff_line(out_file, last_gene_start, last_gene_stop, last_gene_strand, annot)
            current_operon_names = []
            current_operon_start = 0
            current_operon_stop = 0

        current_operon_names.append(name1)
        if not current_operon_start:
            current_operon_start = start1
            current_operon_stop = stop1

        if is_operon:
            if strand1 == '+':
                current_operon_stop = stop2
            elif strand1 == '-':
                current_operon_start = start2
        else:
            if strand1 == '-':
                current_operon_names.reverse()
            operon_name = make_operon_name(current_operon_names)
            annot = 'operonName "' + operon_name + '";'
            write_gff_line(out_file, current_operon_start, current_operon_stop,
                           strand1, annot)
            current_operon_names = []
            current_operon_start = start2
            current_operon_stop = stop2

        last_gene_name = name2
        last_gene_start = start2
        last_gene_stop = stop2
        last_gene_strand = strand2
    out_file.close()


def write_gff_line(file, start, stop, strand, annotation):
    start, stop = sorted([start, stop])
    gff_line_format_list = ["Chromosome", "MicrobesOnline", "operon",
                            "%d", "%d", # Start and Stop
                            ".", "%s", # Score and Strand
                            ".", "%s"] # Frame and annotation
    gff_line_format = "\t".join(gff_line_format_list) + "\n"
    gff_line = gff_line_format % (start, stop, strand, annotation)
    file.write(gff_line)




def main():
    gene_info_file = open('geneInfo.txt')
    operonCalls_file = open('operonCalls.txt')
    start_stop_strand = get_gene_info(gene_info_file)
    operon_calls_to_gtf(operonCalls_file, start_stop_strand)

if __name__ == "__main__":
    main()
