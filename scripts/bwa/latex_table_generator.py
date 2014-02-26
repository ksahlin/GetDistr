import os
import csv
import argparse





# def supplementary_table(line):

# 	line.split()

# def main_table(line):
# 	line.split()

def main(args):
	if args.suppl:
		outfile = os.path.join(args.outpath,'supplemental_table.csv')
	else:
		outfile = os.path.join(args.outpath,'main_table.csv')

	# read tab-delimited file
	with open(args.infile,'rb') as fin:
	    cr = csv.reader(fin, delimiter='\t')
	    filecontents = [line for line in cr]

	# write comma-delimited file (comma is the default delimiter)
	with open(outfile,'wb') as fou:
	    cw = csv.writer(fou, quotechar='', quoting=csv.QUOTE_NONE)
	    cw.writerows(filecontents)

	# for line in args.infile:
	# 	if args.suppl:
	# 		supplementary_table(outfile,line)
	# 	else:
	# 		main_table(outfile,line)





if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Formats outfile of getdistr and bwa experiment to a CSV file\
    	for input to http://www.tablesgenerator.com/latex_tables .")
    parser.add_argument('infile', type=str, help='Simulation file. ')
    parser.add_argument('outpath', type=str, help='Path to output location. ')
    parser.add_argument('suppl', type=bool, help='Main or supplemental table (1 for main 0 for supplemental). ')


    args = parser.parse_args()
    main(args)