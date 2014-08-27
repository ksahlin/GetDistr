'''
Created on Sep 18, 2013

@author: ksahlin
'''

import argparse
from mathstats.normaldist.normal import MaxObsDistr
# class PositionStats:
# 	def get_spanning_reads(self, position, read_length, read_path):
# 	   	observations = [ ]
# 	   	for read in pysam.Samfile( read_path, "rb" ):
# 	   		if not read.is_proper_pair or read.is_read2: # or len( read.cigar ) > 1:
# 	   			continue

# 	       	if read.tlen > 0:
# 	        	if (read.positions[ 0 ] <= position) and ((read.positions[ 0 ] + read.tlen) >= position + read_length):
# 					observations.append( abs( read.tlen ) )
# 	       	else:
# 	           	if (read.positions[ -1 ] >= (position + 1)) and ((read.positions[ - 1 ] + read.tlen) < position):
# 					observations.append( abs( read.tlen ) )

# 		return observations





import sys
import pysam
from scipy.stats import norm, t, chi2
import math
from itertools import ifilter

from getdistr import model

class Parameters(object):
	"""docstring for Parameters"""
	def __init__(self):
		super(Parameters, self).__init__()
		self.mean = None
		self.stddev = None
		self.d = None
		self.pval = None
		self.genome_length = None
		self.nobs = None

	def get_pval_threshold(self):
		mean_sd = t.ppf( 0.975, self.nobs - 1 ) * self.stddev / math.sqrt( self.nobs )
		sd_sd = math.sqrt( (self.nobs - 1) * self.stddev**2 / chi2.ppf( 0.025, self.nobs - 1 ) )

		mean_conservative = self.mean + mean_sd
		sd_conservative = sd_sd

		z = norm.ppf( self.pval / 2.0, mean_conservative, sd_conservative )
		return norm.cdf( z, self.mean, self.stddev )

class ReadContainer(object):
	"""docstring for ReadContainer"""
	def __init__(self, position):
		super(ReadContainer, self).__init__()
		self.position = position
		self.reads = []
	def __str__(self):
		"Prints reads on fasta format"
		pass
	def add_read(self,read):
		self.reads.append(read)

	def print_bam(self):
		"Print reads on bam format "

	def calc_observed_insert(self):
		n = len(self.reads)
		if n == 0:
			return 0

		isize_obs = 0
		modulo = 0
		for read1,read2 in zip(self.reads[:-1],self.reads[1:]):
			if modulo == 1:
				modulo = 0
				continue
			else: 
		   		assert read1.tlen == - read2.tlen
		   		assert read1.is_proper_pair == read2.is_proper_pair
		   		isize_obs += abs(read1.tlen)
				modulo = 1

		return isize_obs / float(len(self.reads)/2)

	
	def calc_pvalue(self,expected_insert,sigma):
		n = len(self.reads)
		if n == 0:
			return 0.5 # return a non significant p-value
	    #insert size
		mean_insert_obs = self.calc_observed_insert()
		# tot_insert = 0
		# for read in self.reads:
		# 	tot_insert += abs(read.tlen)
		# mean_insert_obs = tot_insert/float(n)  
		#print mean_insert_obs

		z = (mean_insert_obs - expected_insert)/(float(sigma)/math.sqrt(n)) # z-statistic
		p_value_upper_area = norm.sf(z) 
		#print z, p_value
		return p_value_upper_area



class BreakPointContainer(object):
	"""docstring for BreakPointContainer"""
	def __init__(self,param):
		super(BreakPointContainer, self).__init__()
		self.clusters = {}
		self.index = 1
		self.clusterinfo = {}
		self.param = param

	def add_bp_to_cluster(self, pos, p_val, nr_obs, mean_obs, sv_type_observed, dist_thresh):
		new_cluster = 1
		for (sv_type,i), cluster in self.clusters.iteritems():
			min_pos = min(cluster)
			max_pos = max(cluster)

			if sv_type_observed != sv_type:
				continue

			if pos <= max_pos + dist_thresh and pos >= min_pos - dist_thresh:
				new_cluster = 0
				cluster.append(pos)
				self.clusterinfo[(sv_type,i)].append((p_val,nr_obs,mean_obs))
				break

		if new_cluster:
			self.clusters[(sv_type_observed, self.index)] = [pos]
			self.clusterinfo[(sv_type_observed, self.index)] = [(p_val,nr_obs,mean_obs)]
			self.index += 1


		if len(self.clusters) == 0:
			self.clusters[(sv_type_observed, self.index)] = [pos]
			self.clusterinfo[(sv_type_observed, self.index)] = [(p_val,nr_obs,mean_obs)]
			self.index += 1

	def get_final_bp_info(self):
		self.final_bps = []
		for region in self.clusters:
			n = len(self.clusters[region])
			median_basepair = self.clusters[region][n/2]
			median_info =  self.clusterinfo[region][n/2]
			cluster_length_span = max(self.clusters[region]) - min(self.clusters[region]) + 1
			self.final_bps.append((region[0],median_basepair,median_info,cluster_length_span,n))

	def __str__(self):
		output_string= '#sv clusters:\n#<type>\t<pos>\t<cluster range>\t<nr sign. pvals in cluster>\t<info on called postion(pval,nr_obs,obs isize)>\n'
		
		for sv_type,median_basepair,median_info,cluster_length_span,n in self.final_bps:
			if median_basepair > self.param.mean and median_basepair < self.param.genome_length - self.param.mean:
				output_string += '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(sv_type,median_basepair,cluster_length_span,n,median_info) 
		return output_string


def AdjustInsertsizeDist(mean_insert, std_dev_insert, insert_list):
    k = MaxObsDistr(len(insert_list), 0.95)
    filtered_list = list(filter((lambda x : (x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert)), insert_list))
    if len(insert_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)


def ParseBAMfile(bamfile,param):
	with pysam.Samfile(bamfile, 'rb') as bam:
	# bam = pysam.Samfile(bamfile, 'rb')
	#bam2 = pysam.Samfile(bamfile, 'rb')


		ref_lengths = bam.lengths

		##
		# create a container holding all interesting reads 
		# for a given position in the reference sequence
		counter =0

	   	counter2= 0 

		container = {}
		for i in range(ref_lengths[0]): # assumes only one chromosome, i.e. a single reference strand
			container[i] = ReadContainer(i)
		param.genome_length = int(ref_lengths[0])

		isize_list = []
		#isize_sum_sq = 0
		#isize_n = 0

		bam_filtered = ifilter(lambda r: r.flag <= 200, bam)
	   	
		for read1 in bam_filtered:
			read2 = next(bam_filtered)
			#print read1
			#print read2
			#print read1.tlen,read2.tlen
	   		assert read1.tlen == - read2.tlen
	   		assert read1.is_proper_pair == read2.is_proper_pair

		# Make sure that bamfile is sorted so that mates come after eachotern
		# Throw assertion error here otherwise.
		# use for read1,read2 in zip(bam1,bam2,2) to iterate over the samfiles simultaneosly
		# to use last aligned coordinate of mate instead of read_length (in case of softclipps)!
	   		#print read.tlen,read.is_proper_pair,read.tlen, read.aend,read.pnext
	   		counter += 1
			if counter % 10000 == 0:
				print '#Processed ', counter, 'reads.'
				#break

	   		#if not read1.is_proper_pair: # or len( read.cigar ) > 1:
	   		#	counter2 +=1
	   		#	continue

			if read1.is_unmapped or read2.is_unmapped:
				continue

			if any(x[ 0 ] in ( 1, 2 ) for x in read1.cigar + read2.cigar):
				continue


	   		## add do insert size distribution calculation if proper pair
			if len( read1.cigar ) == 1 and len( read2.cigar ) == 1:
				isize_list.append(abs(read1.tlen))
	   		#isize_sum_sq += read1.tlen**2
	   		#isize_n += 1


			#print inner_start_pos, inner_end_pos
			if read1.tlen > 0:
				inner_start_pos = read1.aend
				inner_end_pos = read2.pos
				for pos in range(inner_start_pos,inner_end_pos):
					container[pos].add_read(read1)
					container[pos].add_read(read2)
			else:
				inner_start_pos = read2.aend
				inner_end_pos = read1.pos

				for pos in range(inner_start_pos,inner_end_pos):
					container[pos].add_read(read1)
					container[pos].add_read(read2)

		n_isize = float(len(isize_list))
		mean_isize = sum(isize_list)/n_isize
		std_dev_isize =  (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), isize_list))) / (n_isize - 1)) ** 0.5
		print '#Mean before filtering :', mean_isize
		print '#Stddev before filtering: ', std_dev_isize
		extreme_obs_occur = True
		while extreme_obs_occur:
			extreme_obs_occur, filtered_list = AdjustInsertsizeDist(mean_isize, std_dev_isize, isize_list)
			n_isize = float(len(filtered_list))
			mean_isize = sum(filtered_list) / n_isize
			std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n_isize - 1)) ** 0.5
			isize_list = filtered_list

		print '#Mean converged:', mean_isize
		print '#Std_est converged: ', std_dev_isize

	param.nobs = n_isize
	param.mean = mean_isize
	param.stddev = std_dev_isize 

	# print counter2, counter,param.mean,param.stddev
	return container




def get_sv_clusters(container,param):
	tot_mean = 0
	nr_obs =0
	est = model.NormalModel(int(round(param.mean,0)),int(round(param.stddev,0)),100,s_inner=0)
	# exp_stddev =  lagg till test av forvantad standard avvikelse
	exp_insert = float(est.expected_mean(1,param.genome_length))
	#exp_insert =  mu+ (sigma**2)/float(mu + 1)
	print '#Average predicted mean over a base pair under p_0: {0}'.format(exp_insert)

	##
	# need to select a consensus loci for a breakpoint
	# using only insert size 
	sv_container = BreakPointContainer(param)

	pval_threshold = param.get_pval_threshold()
	print "#Adjusted threshold: {0}".format(pval_threshold)
	for bp in range(param.genome_length):
		#print container[bp]
		#print container[bp].calc_insert_pvalue
		p_val_upper_area = container[bp].calc_pvalue(exp_insert,param.stddev)
		#print p_val
		# tot_insert = 0
		# if len(container[bp].reads) > 0:
		# 	for read in container[bp].reads:
		# 		tot_insert += abs(read.tlen)
		# 	tot_mean += tot_insert/len(container[bp].reads)

		#print bp, p_val_upper_area
		obs_mean = container[bp].calc_observed_insert()
		nr_reads_over_bp = len(container[bp].reads)
		if  nr_reads_over_bp > 10:
			nr_obs += 1
			tot_mean += obs_mean


		if p_val_upper_area < pval_threshold:
			sv_container.add_bp_to_cluster(bp,  p_val_upper_area, nr_reads_over_bp, obs_mean, 'deletion', param.d)
			#print 'Significant large position insert size reads: Pos: ',bp, ' p-value = ', p_val_upper_area, 'nr of reads:', len(container[bp].reads), 'obs insert: ', obs_mean,'inferred insertion in donor here'

		if p_val_upper_area > 1 - pval_threshold:
			sv_container.add_bp_to_cluster(bp, 1- p_val_upper_area, nr_reads_over_bp, obs_mean, 'insertion', param.d)
			#print 'Significant small position insert size reads: Pos: ',bp, ' p-value = ', p_val_upper_area, 'nr of reads:', len(container[bp].reads), 'obs insert: ', obs_mean,'inferred deletion in donor here'

		

	print '#Average observed mean over a base pair:{0} (we expect this to be similar to the predicted one if not too many SVs)'.format(tot_mean/nr_obs)

	return sv_container

# def output_breaks(sv_container):
# 	print '#sv clusters:'
# 	for sv_type, cluster in sv_container.clusters.iteritems():
# 		print sv_type, cluster
# 		print sv_container.clusterinfo[sv_type]
# 		# print 




def main(args):
	param = Parameters()
	param.d, param.pval = args.d, args.pval
	container = ParseBAMfile(args.bampath,param)
	sv_container = get_sv_clusters(container,param)

	print '#Estimated library params: mean:{0} sigma:{1}'.format(param.mean,param.stddev)
	print '#Genome length:{0}'.format(param.genome_length)

	sv_container.get_final_bp_info()
	print(sv_container)
	# output_breaks(sv_container)



if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Infer variants with simple p-value test using theory of GetDistr - proof of concept.")
	parser.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	parser.add_argument('pval', type=float, help='p-value threshold for calling a variant. ')
	parser.add_argument('d', type=int, help='distance threshold for clustering close p-values. All significant p-values\
	closer than d will be clustered ')

	# parser.add_argument('mean', type=int, help='mean insert size. ')
	# parser.add_argument('stddev', type=int, help='Standard deviation of insert size ')


	args = parser.parse_args()
	main(args)



        
