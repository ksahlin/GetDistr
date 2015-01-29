import re
import argparse
import os, sys
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
		variant_start = 10000 + i*(20000 + min(0, variant_size)) #0 if insertion else negative
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
	for acc in ['sequence_0-donor','newer_1-donor','sequence_2-donor','sequence_3-donor','sequence_4-donor','newer_2-donor','new1-donor','sequence_0','newer_1','sequence_2','sequence_3','sequence_4','newer_2','new1']:
		variants_TP_FP[acc] = [0,0]
		true_breakpoints[acc] = set()
	return true_breakpoints, variants_TP_FP

def main(args):
	#outfile = open(args.outfile,'a')
	#rootdir = 'C:/Users/sid/Desktop/test'
	prev_variant_size = None
	for subdir, dirs, files in os.walk(args.rootdir):

	    for file in files:
			file_path = os.path.join(subdir, file)
			# print os.path.join(subdir, file)
			result = re.search('[0-9]+',file_path)
			variant_size = int(result.group())
			result = re.search('del',file_path)
			if result:
				variant_size = -variant_size

			true_breakpoints, variants_TP_FP =  initialize_containers()
			true_breakpoints = get_true_breakpoints(true_breakpoints, variant_size)
			variants_TP_FP = compare_misassemblies( open(file_path,'r'), true_breakpoints, variants_TP_FP)

			TPs = 0
			FPs = 0
			for ref_name, (TP,FP) in variants_TP_FP.iteritems():
				TPs += TP
				FPs += FP

			if variant_size == prev_variant_size:
				print ' {0} & {1} \\\ '.format( TPs,FPs)
			else:
				sys.stdout.write('{0} & {1} & {2} &'.format(variant_size, TPs,FPs))
			prev_variant_size = variant_size


			# print os.path.join(subdir, file)
			# result = re.search('[0-9]+',os.path.join(subdir, file))
			# variant_size = int(result.group())

			# true_breakpoints, variants_TP_FP_original,variants_TP_FP_corrected =  initialize_containers()
			# true_breakpoints = get_true_breakpoints(true_breakpoints, variant_size)

			# variants_TP_FP_original = compare_misassemblies( open(file,'r'), true_breakpoints, variants_TP_FP_original)
			# variants_TP_FP_corrected = compare_misassemblies( open(args.clever_corr_breaks,'r'), true_breakpoints, variants_TP_FP_corrected)

			# TPs_orig = 0
			# FPs_orig = 0
			# TPs_corr = 0
			# FPs_corr = 0
			# for (ref_name_orig, (TP_orig,FP_orig)) , (ref_name_corr, (TP_corr,FP_corr)) in zip(variants_TP_FP_original.iteritems(), variants_TP_FP_corrected.iteritems()) :
			# 	TPs_orig += TP_orig
			# 	FPs_orig += FP_orig
			# 	TPs_corr += TP_corr
			# 	FPs_corr += FP_corr

			# print '{0} $ {1} $ {2} $ {3} $ {4} '.format(variant_size, TPs_orig,FPs_orig,TPs_corr,FPs_corr)



if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Parse reapr's results.")   
    parser.add_argument('rootdir', type=str, help="root result folder ")

    #parser.add_argument('clever_breaks', type=str, help="Tool's predicted SVs. ")
    #parser.add_argument('clever_corr_breaks', type=str, help="Tool's predicted SVs. ")
    #parser.add_argument('variant_size', type=int, help=" scaffold fasta file." )
    #parser.add_argument('outfile', type=str, help="Tools results. ")

    args = parser.parse_args()

    main(args)