
import argparse

import calc_pvalues
import lib_est
import get_bp_stats

import os

class Parameters(object):
	"""docstring for Parameters"""
	def __init__(self, mu, sigma, min_isize, max_isize, read_length, pval):
		super(Parameters, self).__init__()
		self.mu = mu
		self.sigma = sigma
		self.min_isize = min_isize
		self.max_isize = max_isize
		self.read_length = read_length
		self.pval = pval

def collect_libstats(infile):
	info_file = open(os.path.join(infile,'info.txt'),'r')
	param = Parameters()
	vals = filter( lambda line: line[0] != '#', info_file.readlines())[0:4]
	print vals
	[mean,stddev] =  vals[0].strip().split()
	[min_lib_isize,max_lib_isize] = vals[1].strip().split()
	[read_length] = vals[2].strip().split()
	[adjusted_mean, adjusted_stddev] = vals[3].strip().split()
	param.mean, param.stddev, param.adjusted_mean,param.adjusted_stddev = float(mean), float(stddev), float(adjusted_mean), float(adjusted_stddev)
	param.min_isize, param.max_isize = int(min_lib_isize), int(max_lib_isize)
	param.pval = p_val_threshold
	param.read_length = float(read_length)
	print param.mean, param.stddev, param.adjusted_mean,param.adjusted_stddev, param.min_isize, param.max_isize, param.read_length
	return param

def main_pipline(args):
	"""
		Algorithm a follows:
			1 Estimate library parameters (lib_est) and print to library_info.txt.
			2 Parse bamfile and get mean and stddev over each position in assembly (get_bp_stats)
				Print to bp_stats.csv
			3 Get gap coordinates in assembly (get_gap_coordinates)
				print to gap_coordinates.csv
			4 Calculate pvalues based on expected span mean and stddev (calc_pvalues)
				print ctg_accesion, pos, pvalues to p_values.csv
			5 Cluster p-values into significant cliques and print significant
				locations on GFF format.


	"""

	# 1
	lib_out = os.path.join(args.outfolder,'library_info.txt')
	lib_est.LibrarySampler(args.bampath,lib_out)
 
 	# 2
	param = collect_libstats(lib_out)
	get_bp_stats(args.bampath, param, os.path.join(args.outfolder,'bp_stats.txt'))

	#3
	cluster_pvals(args.outfolder, args.assembly_file, args.pval, args.window_size)



if __name__ == '__main__':

	# create the top-level parser
	parser = argparse.ArgumentParser()#prog="Infer variants with simple p-value test using theory of GetDistr - proof of concept.")
	#parser.add_argument('--foo', action='store_true', help='help for foo arg.')
	subparsers = parser.add_subparsers(help='help for subcommand')

	# create the parser for the "pipeline" command
	pipeline = subparsers.add_parser('pipeline', help='Run the entire pipeline')
	pipeline.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	pipeline.add_argument('assembly_file', type=str, help='Fasta file with assembly/genome. ')
	pipeline.add_argument('outfolder', type=str, help='Outfolder. ')
	pipeline.add_argument('window_size', type=int, help='Window size ')
	pipeline.add_argument('pval', type=float, help='p-value threshold for calling a variant. ')
	pipeline.set_defaults(which='pipeline')


	# create the parser for the "lib_est" command	
	lib_est = subparsers.add_parser('lib_est', help='Estimate library parameters')
	lib_est.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	get_bp_stats.set_defaults(which='lib_est')
	
	# create the parser for the "get_bp_stats" command
	get_bp_stats = subparsers.add_parser('get_bp_stats', help='Scan bam file and calculate pvalues for each base pair')
	get_bp_stats.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	get_bp_stats.add_argument('outfolder', type=str, help='Outfolder. ')
	get_bp_stats.set_defaults(which='get_bp_stats')


	# create the parser for the "cluster" command
	cluster = subparsers.add_parser('cluster_pvals', help='Takes a pvalue file and clusters them into significan regions')
	cluster.add_argument('assembly_file', type=str, help='Fasta file with assembly/genome. ')
	cluster.add_argument('outfolder', type=str, help='Folder with p-value fila and "info.txt" file contianing parameters "scan_bam" output. ')
	cluster.add_argument('window_size', type=int, help='Window size ')
	cluster.add_argument('pval', type=float, help='p-value threshold for calling a variant. ')
	cluster.set_defaults(which='cluster')

	# parser.add_argument('mean', type=int, help='mean insert size. ')
	# parser.add_argument('stddev', type=int, help='Standard deviation of insert size ')

	
	args = parser.parse_args()
	#print args

	if args.which == 'pipeline':
		main_pipline(args)
	elif args.which == 'scan_bam':
		scan_bam(args.bampath, args.assembly_file, args.outfolder)
	elif args.which == 'cluster':
		cluster_pvals(args.outfolder, args.assembly_file ,args.pval, args.window_size)


	#main(args)