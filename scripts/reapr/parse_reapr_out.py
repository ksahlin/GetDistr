
import argparse


def main(args):
	outfile = open(args.outfile,'a')
	print >> outfile, 'type\t#errors'
	for line in open(args.infile.strip(),'r'):
		if line[0] == '#':
			continue
		else:
			scf_name = line.split()[0]
			nr_errors = line.split()[12]
			errors_length = line.split()[13]
		print >> outfile, scf_name+'\t'+nr_errors + '\t'+errors_length
if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Parse reapr's results.")
    parser.add_argument('infile', type=str, help='Reaprs infile. ')
    parser.add_argument('outfile', type=str, help='outfile. ')

    args = parser.parse_args()

    main(args)