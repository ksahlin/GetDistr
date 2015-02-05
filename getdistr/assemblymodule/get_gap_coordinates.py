

import re

def ReadInContigseqs(contigfile):
    cont_dict = {}
    k = 0
    temp = ''
    accession = ''
    for line in contigfile:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            cont_dict[accession] = ''
            k += 1
        elif line[0] == '>':
            cont_dict[accession] = temp
            temp = ''
            accession = line[1:].strip().split()[0]
        else:
            temp += line.strip()
    cont_dict[accession] = temp

    return(cont_dict)


def print_gaps(gap_coordinates, outfile):
	for ref in gap_coordinates:
		for coord in gap_coordinates[ref]:
			print >> outfile, '{0}\t{1}\t{2}'.format(ref,coord[0],coord[1])

def get_gap_coordinates(fasta_file,outpath):
	ctg_dict = ReadInContigseqs(open(fasta_file,'r'))

	#p = re.compile("[Nn]+")
	# at least 20 N's to consider to be a gap
	p = re.compile("[Nn]{20,}")
	gap_coordinates = {}

	for acc,ctg in ctg_dict.iteritems():
		gap_coordinates[acc] = set()
		for m in p.finditer(ctg):
			gap_coordinates[acc].add((m.start() ,m.end()))

	outfile = open(outpath, 'w')
	print_gaps(gap_coordinates, outfile)

	return gap_coordinates