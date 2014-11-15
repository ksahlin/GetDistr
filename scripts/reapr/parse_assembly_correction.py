
import gzip
import argparse

def is_true_positive(pred_start, pred_stop, true_breakpoints ):
	for start,stop in true_breakpoints:
		## in true expansion interval or contraction interval of a gap
		if start <= pred_start <= stop or start <= pred_stop <= stop:
			return True
		# a contracted break point in predicted interval of reapr
		elif pred_start <= start <= pred_stop or pred_start <= stop <= pred_stop:
			return True  

	return False

def compare_misassemblies(infile, true_breakpoints):
	scaffold_tp_fp = {}
	for line in infile:
		print line
		values = line.strip().split()
		scaf_name = values[0]

		if scaf_name not in scaffold_tp_fp:
			scaffold_tp_fp[scaf_name] = [0,0] # TP, FP
		
		if values[2] == 'FCD':
			start,stop = int(line.strip().split()[3]), int(line.strip().split()[4])
			if is_true_positive:
				scaffold_tp_fp[scaf_name][0] += 1
			else:
				scaffold_tp_fp[scaf_name][1] += 1

	return scaffold_tp_fp	


def get_true_breakpoints(infile):
	true_breakpoints = set()
	for line in infile:
		start,stop = int(line.strip().split()[3]), int(line.strip().split()[4])
		true_breakpoints.add((start,stop))
	return true_breakpoints

def main(args):
	outfile = open(args.outfile,'a')
	true_breakpoints = get_true_breakpoints(open(args.true_gff, 'r'))
	if args.tools_gff[-3:] == '.gz':
		tool_results = open(args.tools_gff, 'rb')
	else:
		tool_results = open(args.tools_gff, 'r')
	scaffold_tp_fp = compare_misassemblies( tool_results, true_breakpoints)

	for scaf_name, (TP,FP) in scaffold_tp_fp.iteritems():
		print '{0}\t{1}\t{2}'.format(scaf_name, TP,FP)
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
    #parser.add_argument('reapr_outfile', type=str, help=" Reapr's utfile." )
    parser.add_argument('outfile', type=str, help="Tools results. ")

    args = parser.parse_args()

    main(args)