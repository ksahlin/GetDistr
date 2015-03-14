
import argparse

from getdistr.assemblymodule import lib_est
from getdistr.assemblymodule import get_bp_stats
from getdistr.assemblymodule import get_gap_coordinates
from getdistr.assemblymodule import filter_bam
from getdistr.assemblymodule import cluster_p_vals

import os,sys

class Parameters(object):
	"""docstring for Parameters"""
	def __init__(self):
		super(Parameters, self).__init__()
		self.mu = None
		self.sigma = None
		self.adjusted_mu = None
		self.adjusted_sigma = None
		self.min_isize = None
		self.max_isize = None
		self.read_length = None
		self.pval = None
		self.total_basepairs = None
		self.ess_ratio = None

		self.scaf_lengths = {}
		self.outfolder = None
		self.plots = None
		self.stats_output = None

		self.nr_reads = 0
		self.nr_mapped = 0
		self.nr_proper_mapped = 0
		# self.lib_min = None
		# self.lib_max = None

def collect_libstats(args,outfolder,param):
	info_file = open(os.path.join(outfolder,'library_info.txt'),'r')
	vals = filter( lambda line: line[0] != '#', info_file.readlines())[0:]
	print vals
	[mean,stddev] =  vals[0].strip().split()
	[min_lib_isize,max_lib_isize] = vals[1].strip().split()
	[read_length] = vals[2].strip().split()
	[adjusted_mean, adjusted_stddev] = vals[3].strip().split()
	param.mu, param.sigma, param.adjusted_mu,param.adjusted_sigma = float(mean), float(stddev), float(adjusted_mean), float(adjusted_stddev)
	param.min_isize, param.max_isize = int(min_lib_isize), int(max_lib_isize)
	param.read_length = float(read_length)
	print param.mu, param.sigma, param.adjusted_mu, param.adjusted_sigma, param.min_isize, param.max_isize, param.read_length
	param.max_window_size = int(param.mu) if param.mu <= 1000 else int(param.mu)/2
	param.ess_ratio = float(vals[4].strip().split()[0])
	param.total_basepairs = int(vals[5].strip().split()[0])
	param.pval = 0.05/ param.total_basepairs # bonferroni correction

	
	for line in vals[6:]:
		ref,length = line.strip().split('\t')
		length = int(length)
		param.scaf_lengths[ref] = length


def filter_bamfile(args,param):
	param.lib_min = args.lib_min
	param.lib_max = args.lib_max
	bam_out = os.path.join(args.outfolder,'bam_filtered.bam')
	filter_bam.within_reference(args.bampath, bam_out, args.n, param )

def get_lib_est(bam_file,param):
	#bam_in = os.path.join(args.outfolder,'bam_filtered.bam')
	lib_est.LibrarySampler(bam_file, param)

def bp_stats(args,param):
	bam_in = os.path.join(args.outfolder,'bam_filtered.bam')
	collect_libstats(args,args.outfolder, param)
	get_bp_stats.parse_bam(bam_in, param)

def gap_coordinates(args,param):
	gaps_out = os.path.join(args.outfolder,'gap_coordinates.txt')
	get_gap_coordinates.get_gap_coordinates(args.assembly_file, gaps_out)

def p_value_cluster(args,param):
	bp_file = os.path.join(args.outfolder,'bp_stats.txt')
	gap_file = os.path.join(args.outfolder,'gap_coordinates.txt')
	collect_libstats(args,args.outfolder,param)
	cluster_p_vals.main(bp_file, gap_file,param)

def main_pipline(args,param):
	"""
		Algorithm a follows:
			1 Filter bam file to only consider interesting reads
			2 Estimate library parameters (lib_est) and print to library_info.txt.
				(mu, sigam, adjusted mu, sigma, ESS etc.)
			3 Parse bamfile and get mean and stddev over each position in assembly (get_bp_stats)
				Print to bp_stats.csv
			4 Get gap coordinates in assembly (get_gap_coordinates)
				print to gap_coordinates.csv 
			5 Calculate pvalues based on expected span mean and stddev (calc_pvalues)
				print ctg_accesion, pos, pvalues to p_values.csv
			5' Cluster p-values into significant cliques and print significant
				locations on GFF format.


	"""

	if not os.path.exists(args.outfolder):
		os.makedirs(args.outfolder)

	# 1
	bam_out = os.path.join(args.outfolder,'bam_filtered.bam')
	filter_bamfile(args,param)

	# 2
	lib_est.LibrarySampler(bam_out,param)
 
 	# 3
	collect_libstats(args,args.outfolder,param)
	get_bp_stats.parse_bam(bam_out, param)

	# 4
	#get_gap_coordinates.
	gap_coordinates(args,param)

	# 5-5' 
	p_value_cluster(args,param)



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
	pipeline.add_argument('--n', dest='n', type=int, default=1, help='Neighborhood size. ')	
	pipeline.add_argument('--lib_min', dest='lib_min', type=int, default=200, help='Minimum insert size (if in doubt, just set lib_min = 2*read_length). ')	
	pipeline.add_argument('--lib_max', dest='lib_max', type=int, default=200000, help='Maximum insert size (tight bound is not necessary, choose a larger value rather than smaller). ')	
	pipeline.add_argument('--plots', dest="plots", action='store_true', help='Plot pval distribution.')
	pipeline.set_defaults(which='pipeline')

	# create the parser for the "filter" command	
	filter_parser = subparsers.add_parser('filter', help='Filters bam file for better uniform coverage.')
	filter_parser.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	filter_parser.add_argument('outfolder', type=str, help='Outfolder. ')
	filter_parser.add_argument('--n', dest='n', type=int, default=1, help='Neighborhood size. ')	
	filter_parser.add_argument('--lib_min', dest='lib_min', type=int, default=200, help='Minimum insert size (if in doubt, just set lib_min = 2*read_length). ')	
	filter_parser.add_argument('--lib_max', dest='lib_max', type=int, default=200000, help='Maximum insert size (tight bound is not necessary, choose a larger value rather than smaller). ')	
	filter_parser.add_argument('--plots', dest="plots", action='store_true', help='Plot isize distribution.')
	filter_parser.set_defaults(which='filter')

	# create the parser for the "lib_est" command	
	lib_est_parser = subparsers.add_parser('lib_est', help='Estimate library parameters')
	lib_est_parser.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	lib_est_parser.add_argument('outfolder', type=str, help='Outfolder. ')
	lib_est_parser.add_argument('--plots', dest="plots", action='store_true', help='Plot pval distribution.')
	lib_est_parser.set_defaults(which='lib_est')
	
	# create the parser for the "get_bp_stats" command
	get_bp_stats_parser = subparsers.add_parser('get_bp_stats', help='Scan bam file and calculate pvalues for each base pair')
	get_bp_stats_parser.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	get_bp_stats_parser.add_argument('outfolder', type=str, help='Outfolder. ')
	get_bp_stats_parser.add_argument('--plots', dest="plots", action='store_true', help='Plot pval distribution.')
	get_bp_stats_parser.set_defaults(which='get_bp_stats')

	# create the parser for the "get_gaps" command
	get_gaps_stats_parser = subparsers.add_parser('get_gaps', help='Contig assembly file for gap positions')
	get_gaps_stats_parser.add_argument('assembly_file', type=str, help='Fasta file with assembly/genome. ')
	get_gaps_stats_parser.add_argument('outfolder', type=str, help='Outfolder. ')
	get_gaps_stats_parser.add_argument('--plots', dest="plots", action='store_true', help='Plot pval distribution.')
	get_gaps_stats_parser.set_defaults(which='get_gaps')

	# create the parser for the "pvalue_cluster" command
	pvalue_cluster_parser = subparsers.add_parser('cluster_pvals', help='Takes a pvalue file and clusters them into significan regions')
	pvalue_cluster_parser.add_argument('bp_file', type=str, help='bp_stats.txt file generated by "get_bp_stats" command. ')	
	pvalue_cluster_parser.add_argument('gap_file', type=str, help='gap_coordinates.txt file generated by "get_bp_stats" command. ')	
	pvalue_cluster_parser.add_argument('outfolder', type=str, help='Folder with p-value fila and "info.txt" file contianing parameters "scan_bam" output. ')
	pvalue_cluster_parser.add_argument('--plots', dest="plots", action='store_true', help='Plot pval distribution.')
	pvalue_cluster_parser.set_defaults(which='cluster_pvals')

	# parser.add_argument('mean', type=int, help='mean insert size. ')
	# parser.add_argument('stddev', type=int, help='Standard deviation of insert size ')

	
	args = parser.parse_args()

	if args.which == 'pipeline' or args.which == 'get_bp_stats' or args.which == 'lib_est':
		try:
		    open(args.bampath)
		except IOError as e:
		    sys.exit("couldn't find BAM file: " + args.bampath + " check that the path is correct and that the file exists")
		try:
		    open(args.bampath + '.bai')
		except IOError as e:
		    print "couldn't find index file: ", args.bampath + '.bai', " check that the path is correct and that the bam file is sorted and indexed"
		    sys.exit(0)
	#print args
	param = Parameters()
	param.plots = args.plots
	param.outfolder = args.outfolder
	param.stats_output = os.path.join(param.outfolder,'stats_output.txt')
	# param.lib_min = args.lib_min
	# param.lib_max = args.lib_max


	if not os.path.exists(param.outfolder):
		os.makedirs(param.outfolder)
	if param.plots:
		param.plotfolder = os.path.join(param.outfolder,'plots')

		if not os.path.exists(param.plotfolder):
			os.makedirs(param.plotfolder)


	if args.which == 'pipeline':
		main_pipline(args,param)
	elif args.which == 'filter':
		filter_bamfile(args,param)
	elif args.which == 'lib_est':
		get_lib_est(args.bampath,param)
	elif args.which == 'get_bp_stats':
		bp_stats(args,param)
	elif args.which == 'get_gaps':
		gap_coordinates(args,param)
	elif args.which == 'cluster_pvals':
		p_value_cluster(args,param)
	else:
		print 'invalid call'


