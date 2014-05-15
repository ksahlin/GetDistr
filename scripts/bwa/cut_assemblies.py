import argparse

def calc_Nxx_contig_length(Nxx,ctgs):
    lengths = map(lambda x: len(ctgs[x]),ctgs )
    lengths.sort(reverse=True)
    N100 = sum(lengths)
    stop = N100/(100/float(Nxx))
    curr_sum = 0
    for length in lengths:
        curr_sum += length
        if curr_sum >= stop:
        	print 'N'+str(Nxx)+':', length
        	return length


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

def write_filtered_contigs(outfile,ctgs,length):

	for acc,ctg in ctgs.iteritems():
		if len(ctg) < length:
		#	print len(ctg)
			print >>outfile, '>'+acc+'\n'+ctg


def main(args):
	ctgs = ReadInContigseqs(open(args.infile,'r'))
	length =calc_Nxx_contig_length(args.n,ctgs)
	write_filtered_contigs(open(args.outfile,'w'),ctgs,length)


if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Remove all contigs larger than the NXX value for a given assembly.")
    parser.add_argument('infile', type=str, help='Assembly fasta file ')
    parser.add_argument('n', type=int, help='The X of NX value (integer 0-100) ')
    parser.add_argument('outfile', type=str, help='Path and filename to output location of reduced assembly. ')


    args = parser.parse_args()
    main(args)