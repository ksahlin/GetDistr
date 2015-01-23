import re
import gzip
import argparse
from genomics_tools.file_formats.fasta import fasta_iter

def is_true_positive(pred_scaf,pred_start, pred_stop, true_breakpoints ):
	for start, stop in true_breakpoints[pred_scaf]:
		## in true expansion interval or contraction interval of a gap
		if (start <= pred_start <= stop) or (start <= pred_stop <= stop):
			return True
		# a contracted break point in predicted interval 
		elif (pred_start <= start <= pred_stop) or (pred_start <= stop <= pred_stop):
			return True  

	return False

def compare_misassemblies(infile, true_breakpoints, scaffold_tp_fp):
	for line in infile:

		values = line.strip().split()

		if values:
			start, stop = int(values[1]), int(values[2])
			ref_name = values[0]

			if is_true_positive(ref_name, start, stop, true_breakpoints):
				scaffold_tp_fp[ref_name][0] += 1
			else:
				scaffold_tp_fp[ref_name][1] += 1

	return scaffold_tp_fp	


def get_true_breakpoints(true_breakpoints,variant_size):
	
	for i  in range(100):
		variant_start = 10000 + i*(8000 + min(0, variant_size)) #0 if insertion else negative
		if variant_size > 0: # insertion in donor
			variant_stop = variant_start + 1
		else: #deletion in donor
			variant_stop = variant_start + abs(variant_size)

		#print i, variant_start, variant_stop
		#scaf_name = line.strip().split()[0]
		#start,stop = int(line.strip().split()[3]), int(line.strip().split()[4])
		true_breakpoints['sequence_0-donor'].add((variant_start, variant_stop))
	return true_breakpoints

def initialize_containers():
	variants_TP_FP = {}
	true_breakpoints = {}
	for acc in ['sequence_0-donor','newer_1-donor','sequence_2-donor','sequence_3-donor','sequence_4-donor','newer_2-donor','new1-donor']:
		variants_TP_FP[acc] = [0,0]
		true_breakpoints[acc] = set()
	return true_breakpoints, variants_TP_FP

def main(args):
	#outfile = open(args.outfile,'a')
	true_breakpoints, variants_TP_FP =  initialize_containers()
	true_breakpoints = get_true_breakpoints(true_breakpoints, args.variant_size)

	variants_TP_FP = compare_misassemblies( open(args.clever_breaks,'r'), true_breakpoints, variants_TP_FP)
	print variants_TP_FP
	# for ref_name, (TP,FP) in variants_TP_FP.iteritems():
	# 	result = re.search('gap[0-9]+',scaf_name)
	# 	if result:
	# 		gap = result.group(0)[3:]
	# 	else:
	# 		continue
		
	# 	result = re.search('errorsize[0-9]+',scaf_name)
	# 	if result:
	# 		error = result.group(0)[9:]
	# 	else:
	# 		result = re.search('minus[0-9]+',scaf_name)
	# 		error = '-'+result.group(0)[5:]

	# 	print '{0}\t{1}\t{2}\t{3}'.format(gap, error, TP,FP)

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Parse reapr's results.")    
    parser.add_argument('clever_breaks', type=str, help="Tool's predicted SVs. ")
    parser.add_argument('variant_size', type=int, help=" scaffold fasta file." )
    #parser.add_argument('outfile', type=str, help="Tools results. ")

    args = parser.parse_args()

    main(args)