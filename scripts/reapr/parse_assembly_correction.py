import re
import gzip
import argparse
from genomics_tools.file_formats.fasta import fasta_iter

def is_true_positive(pred_scaf,pred_start, pred_stop, true_breakpoints ):
	for start, stop in true_breakpoints[pred_scaf]:
		## in true expansion interval or contraction interval of a gap
		if start <= pred_start <= stop or start <= pred_stop <= stop:
			return True
		# a contracted break point in predicted interval of reapr
		elif pred_start <= start <= pred_stop or pred_start <= stop <= pred_stop:
			return True  

	return False

def compare_misassemblies(scafs, infile, true_breakpoints, scaffold_tp_fp):
	for line in infile:
		values = line.strip().split()
		if values:
			scaf_name = values[0]
			
			if values[2] == 'FCD' or values[2] == 'FCD_gap':
				start,stop = int(line.strip().split()[3]), int(line.strip().split()[4])
				if is_true_positive(scaf_name, start, stop,true_breakpoints):
					scaffold_tp_fp[scaf_name][0] += 1
				else:
					scaffold_tp_fp[scaf_name][1] += 1

	return scaffold_tp_fp	


def get_true_breakpoints(infile,true_breakpoints):
	for line in infile:
		scaf_name = line.strip().split()[0]
		start,stop = int(line.strip().split()[3]), int(line.strip().split()[4])
		true_breakpoints[scaf_name].add((start,stop))
	return true_breakpoints

def initialize_containers(args):
	scaffold_tp_fp = {}
	true_breakpoints = {}
	for acc, scf in fasta_iter(open(args.scafs, 'r')):
		scaffold_tp_fp[acc] = [0,0]
		true_breakpoints[acc] = set()
	return true_breakpoints, scaffold_tp_fp

def main(args):
	#outfile = open(args.outfile,'a')
	true_breakpoints, scaffold_tp_fp =  initialize_containers(args)
	true_breakpoints = get_true_breakpoints(open(args.true_gff, 'r'), true_breakpoints)

	if args.tools_gff[-3:] == '.gz':
		tool_results = gzip.open(args.tools_gff, 'rb')
	else:
		tool_results = open(args.tools_gff, 'r')
	scaffold_tp_fp = compare_misassemblies( args.scafs, tool_results, true_breakpoints, scaffold_tp_fp)

	for scaf_name, (TP,FP) in scaffold_tp_fp.iteritems():
		#scf_gap100_errorsize_minus50
		#scf_gap0_errorsize_50
		result = re.search('gap[0-9]+',scaf_name)
		if result:
			gap = result.group(0)[3:]
		else:
			continue
		
		result = re.search('errorsize[0-9]+',scaf_name)
		if result:
			error = result.group(0)[9:]
		else:
			result = re.search('minus[0-9]+',scaf_name)
			error = '-'+result.group(0)[5:]

		print '{0}\t{1}\t{2}\t{3}'.format(gap, error, TP,FP)
		#print '{0}\t{1}\t{2}'.format(scaf_name, TP,FP)
# 	print >> outfile, 'type\t#errors'
# 	for line in open(args.infile.strip(),'r'):
# 		if line[0] == '#':
# 			continue
# 		else:
# 			scf_name = line.split()[0]
# 			nr_errors = line.split()[12]
# 			errors_length = line.split()[13]
# 		print >> outfile, scf_name+'\t'+nr_errors + '\t'+errors_length

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Parse reapr's results.")
    parser.add_argument('true_gff', type=str, help="True misassemblies. ")    
    parser.add_argument('tools_gff', type=str, help="Tool's predicted misassemblies. ")
    parser.add_argument('scafs', type=str, help=" scaffold fasta file." )
    #parser.add_argument('outfile', type=str, help="Tools results. ")

    args = parser.parse_args()

    main(args)