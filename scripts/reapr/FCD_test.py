'''
Created on Sep 18, 2013

@author: ksahlin
'''

import argparse
import os
from mathstats.normaldist.normal import MaxObsDistr
from scipy.stats import ks_2samp,norm
import random
import re
import math

import pysam
#import math
from itertools import ifilter
import model
import bisect

from Queue import Queue

import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF

EMPIRICAL_BINS = 200
SAMPLE_SIZE = 50000  # for estimating true full read pair distribution

def is_proper_aligned_unique_innie(read):
    return (read.is_reverse and not read.mate_is_reverse and read.is_read1 and read.tlen < 0 and read.rname == read.mrnm) or \
                (not read.is_reverse and read.mate_is_reverse and read.is_read1 and read.tlen > 0 and read.rname == read.mrnm ) \
                and not read.mate_is_unmapped and read.mapq > 10 and not read.is_secondary

def ReadInContigseqs(contigfile,contig_filter_length):
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

    # for contig in cont_dict.keys():
    #     print 'Initial length:', len(cont_dict[contig])
    if contig_filter_length:
        singled_out = 0
        for contig in cont_dict.keys():
            if len(cont_dict[contig]) < contig_filter_length:
                del cont_dict[contig]
                singled_out += 1
    return(cont_dict)


class Parameters(object):
	"""docstring for Parameters"""
	def __init__(self):
		super(Parameters, self).__init__()
		self.mean = None
		self.stddev = None
		self.d = None
		self.pval = None
		self.scaf_lengths = {}
		self.nobs = None
		self.true_distr = None
		self.corrected_pval = None

	def get_pval_threshold(self):
		mean_sd = t.ppf( 0.975, self.nobs - 1 ) * self.stddev / math.sqrt( self.nobs )
		sd_sd = math.sqrt( (self.nobs - 1) * self.stddev**2 / chi2.ppf( 0.025, self.nobs - 1 ) )

		mean_conservative = self.mean + mean_sd
		sd_conservative = sd_sd

		z = norm.ppf( self.pval / 2.0, mean_conservative, sd_conservative )
		return norm.cdf( z, self.mean, self.stddev )

	def sample_distribution(self,bamfile,outfile):
		isize_list = []
		#i = 0
		bam_filtered = ifilter(lambda r: is_proper_aligned_unique_innie(r), bamfile)
		#while i <= sample_size:
		for sample_nr,read in enumerate(bam_filtered):
	   		## add do insert size distribution calculation if proper pair
			if is_proper_aligned_unique_innie(read):
				isize_list.append(abs(read.tlen))
				#sample_nr+=1
			if sample_nr > SAMPLE_SIZE:
				break
		print >> outfile, 'Insert size sample size:', sample_nr
		bamfile.reset()

		n_isize = float(len(isize_list))
		mean_isize = sum(isize_list)/n_isize
		std_dev_isize =  (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), isize_list))) / (n_isize - 1)) ** 0.5
		print >> outfile,'#Mean before filtering :', mean_isize
		print >> outfile,'#Stddev before filtering: ', std_dev_isize
		extreme_obs_occur = True
		while extreme_obs_occur:
			extreme_obs_occur, filtered_list = AdjustInsertsizeDist(mean_isize, std_dev_isize, isize_list)
			n_isize = float(len(filtered_list))
			mean_isize = sum(filtered_list) / n_isize
			std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n_isize - 1)) ** 0.5
			isize_list = filtered_list

		print >> outfile,'#Mean converged:', mean_isize
		print >> outfile,'#Std_est converged: ', std_dev_isize

		self.nobs = n_isize
		self.mean = mean_isize
		self.stddev = std_dev_isize 
		self.full_ECDF = ECDF(isize_list)
		self.adjustedECDF_no_gap = None
		self.get_true_normal_distribution(random.sample(isize_list, min(10000,sample_nr)),outfile)


	def get_true_normal_distribution(self,sample,outfile):
		read_len = 100
		softclipps = 0
		#self.temp_cdf = {}
		#self.temp_cdf[norm.cdf(2*(read_len-softclipps), self.mean, self.stddev) * 1] = 2*(read_len-softclipps)  # weight one for one placement
		cdf_list =[norm.cdf(2*(read_len-softclipps), self.mean, self.stddev) * 1]
		for x in range(2*(read_len-softclipps) +1,int(self.mean + 6*self.stddev)):
			increment_area = norm.pdf(x, self.mean, self.stddev) * (x-(2*(read_len-softclipps)-1))
			#self.temp_cdf[sum(cdf_list) + increment_area] = x
			cdf_list.append( cdf_list[-1] + increment_area)
		tot_cdf = cdf_list[-1]
		cdf_list_normalized = map(lambda x: x /float(tot_cdf),cdf_list)

		# Now create a weighted sample
		self.true_distr = []
		for i in range(1000):
			obs = random.uniform(0, 1)
			pos = bisect.bisect(cdf_list_normalized, obs) - 1
			#print obs, pos
			self.true_distr.append(pos + 2*(read_len-softclipps))


		n = len(self.true_distr)
		self.adjusted_mean = sum(self.true_distr)/float(len(self.true_distr))
		self.adjusted_stddev = (sum(list(map((lambda x: x ** 2 - 2 * x * self.adjusted_mean + self.adjusted_mean ** 2), self.true_distr))) / (n - 1)) ** 0.5

		print >> outfile,'Corrected mean:{0}, corrected stddev:{1}'.format(self.adjusted_mean, self.adjusted_stddev)
		# for x in range(2*(read_len-softclipps) +1,self.mean + 6*self.stddev):

		# norm.pdf()
		# min_weight, max_weight = observations[0], observations[-1]
		
		# total_weight = 


	#def cdf(x):
	#	return self.cdf


	def get_weight(self,x,gap_coordinates,r,s):
		if not gap_coordinates:
			return x - (2*(r-s)-1)
		total_restriction_positions_left = 0
		total_restriction_positions_right = 0

		for start,stop in gap_coordinates:
			if  x < start:
				continue #total_restriction_positions_right += 0

			elif  (r-s)-1 <= start  <= x:
				#print '1'
				total_restriction_positions_right += min(stop + (r-s)-1, x) - start  + (r-s)-1
			
			elif 0 <= start <= (r-s)-1:
				#print '2'
				total_restriction_positions_right += min(stop + (r-s)-1, x) 

			elif stop < -x:
				total_restriction_positions_left += (r-s)

			elif -x <= stop < -( (r-s)-1):
				#print '3'
				total_restriction_positions_left += stop - max(start - (r-s)-1 ,-x)  + (r-s)-1
				
			elif -( (r-s)-1) <= stop <= 0:
				#print '4'
				total_restriction_positions_left += -(max(start,-x))

			elif start <0 and stop > 0:
				#print '5'
				total_restriction_positions_right +=  stop + (r-s)
				total_restriction_positions_left += - start + (r-s)

		#print 'tot restrict:', total_restriction_positions
		total_restriction_positions_right = max(total_restriction_positions_right,s)
		total_restriction_positions_left = max(total_restriction_positions_left,s)
		weight = x - total_restriction_positions_right - total_restriction_positions_left

		return max( 0 , weight)



			# elif stop <= 0 and -stop < x < -start:
			# 	pass
			# elif start <0 and stop > 0:
			# 	pass
		
		return 

	def get_correct_ECDF(self,outfile, gap_coordinates):
		if not gap_coordinates and self.adjustedECDF_no_gap:
			return self.adjustedECDF_no_gap


		read_len = 100
		softclipps = 0


		x_min = max(2*(read_len-softclipps) , int(self.mean - 5*self.stddev) )
		x_max = int(self.mean + 5*self.stddev)
		stepsize =  (x_max - x_min) / EMPIRICAL_BINS
		cdf_list = [ self.full_ECDF( x_min) * self.get_weight(2*(read_len-softclipps), gap_coordinates, read_len, softclipps)  ] #[ self.full_ECDF( 2*(read_len-softclipps)) * self.get_weight(2*(read_len-softclipps), gap_coordinates, read_len, softclipps) ]


		for x in range( x_min + stepsize , x_max, stepsize):
			increment_area = self.get_weight(x,gap_coordinates, read_len, softclipps) * (self.full_ECDF(x) - self.full_ECDF(x-stepsize))
			#increment_area = norm.pdf(x, self.mean, self.stddev) * (x-(2*(read_len-softclipps)-1))
			cdf_list.append( cdf_list[-1] + increment_area)

		#print 'stepsize:', stepsize
		#print 'BINS:',len(cdf_list)
		tot_cdf = cdf_list[-1]
		cdf_list_normalized = map(lambda x: x /float(tot_cdf),cdf_list)

		# Now create a weighted sample
		self.true_distr = []
		for i in range(1000):
			obs = random.uniform(0, 1)
			pos = bisect.bisect(cdf_list_normalized, obs) - 1
			#print obs, pos
			self.true_distr.append(pos*stepsize + x_min)

		# initialization of no gap true distribution
		if not gap_coordinates:
			print 'getting initial gap free distr.'
			self.adjustedECDF_no_gap = self.true_distr

		n = len(self.true_distr)
		self.adjusted_mean = sum(self.true_distr)/float(len(self.true_distr))
		self.adjusted_stddev = (sum(list(map((lambda x: x ** 2 - 2 * x * self.adjusted_mean + self.adjusted_mean ** 2), self.true_distr))) / (n - 1)) ** 0.5

		#print 'Corrected mean:{0}, corrected stddev:{1}, gap_coordinates: {2}'.format(self.adjusted_mean, self.adjusted_stddev, gap_coordinates)

		#print >> outfile,'Corrected mean:{0}, corrected stddev:{1}, gap_coordinates: {2}'.format(self.adjusted_mean, self.adjusted_stddev, gap_coordinates)
		return self.true_distr

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


	def write_pval_to_file(self,outfile,ref_name):
		print >> outfile, '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(ref_name, self.position, self.pval,len(self.isize_list), self.mean_isize, self.std_dev_isize)

	def calc_observed_insert(self):
		#print self.reads
		#print 'list:',map(lambda x: x.tlen, self.reads )
		self.mean_isize = 0
		self.std_dev_isize = 0
		self.isize_list = map(lambda x: x.tlen, self.reads )
		#print len( self.isize_list)
		if len( self.isize_list) <= 5:
			self.mean_isize = -1
			return

		n_isize = float(len(self.isize_list))
		self.mean_isize = sum(self.isize_list)/n_isize
		self.std_dev_isize =  (sum(list(map((lambda x: x ** 2 - 2 * x * self.mean_isize + self.mean_isize ** 2), self.isize_list))) / (n_isize - 1)) ** 0.5


	
	# def calc_pvalue(self,expected_insert,sigma):
	# 	n = len(self.reads)
	# 	if n == 0:
	# 		return 0.5 # return a non significant p-value
	#     #insert size
	# 	mean_insert_obs = self.calc_observed_insert()
	# 	# tot_insert = 0
	# 	# for read in self.reads:
	# 	# 	tot_insert += abs(read.tlen)
	# 	# mean_insert_obs = tot_insert/float(n)  
	# 	#print mean_insert_obs

	# 	z = (mean_insert_obs - expected_insert)/(float(sigma)/math.sqrt(n)) # z-statistic
	# 	p_value_upper_area = norm.sf(z) 
	# 	#print z, p_value
	# 	return p_value_upper_area

	def calc_ks_test(self,true_distribution):
		if len(self.isize_list) >= 5:
			KS_statistic, self.pval = ks_2samp(self.isize_list, true_distribution)
			return KS_statistic, self.pval 
		else:
			self.pval = -1
			return -1, -1


class BreakPointContainer(object):
	"""docstring for BreakPointContainer"""
	def __init__(self,param):
		super(BreakPointContainer, self).__init__()
		self.clusters = {}
		self.index = 1
		self.clusterinfo = {}
		self.param = param

	def add_bp_to_cluster(self, scf, pos, p_val, nr_obs, mean_obs, sv_type_observed, window_size):
		new_cluster = 1
		for (sv_type, ref, i), cluster in self.clusters.iteritems():
			if scf != ref:
				continue
			min_pos = min(cluster)
			max_pos = max(cluster)

			if sv_type_observed != sv_type:
				continue

			if pos <= max_pos + window_size and pos >= min_pos - window_size:
				new_cluster = 0
				cluster.append(pos)
				self.clusterinfo[(sv_type,scf,i)].append((p_val,nr_obs,mean_obs))
				break

		if new_cluster:
			self.clusters[(sv_type_observed, scf, self.index)] = [pos]
			self.clusterinfo[(sv_type_observed, scf, self.index)] = [(p_val,nr_obs,mean_obs)]
			self.index += 1


		if len(self.clusters) == 0:
			self.clusters[(sv_type_observed, scf, self.index)] = [pos]
			self.clusterinfo[(sv_type_observed, scf, self.index)] = [(p_val,nr_obs,mean_obs)]
			self.index += 1

	def get_final_bp_info(self,window_size):
		self.final_bps = []
		for region in self.clusters:
			start_pos = self.clusters[region][0] - window_size
			end_pos = self.clusters[region][-1]
			#n = len(self.clusters[region])
			#median_basepair = self.clusters[region][n/2]
			region_pvalues = map(lambda x:x[0], self.clusterinfo[region])
			avg_region_pvalue = sum(region_pvalues)/len(region_pvalues)
			region_nr_obs = map(lambda x:x[1], self.clusterinfo[region])
			avg_region_nr_obs = sum(region_nr_obs)/len(region_nr_obs)
			region_mean_obs = map(lambda x:x[2], self.clusterinfo[region])
			avg_region_mean_obs = sum(region_pvalues)/len(region_mean_obs)
			# if self.clusters[region][0] >= (start_pos  + end_pos)/2:
			# 	median_info =  self.clusterinfo[region][0]
			# else:
			# 	median_pos = int((start_pos  + end_pos)/2)
			# 	i = 0
			# 	curr_pos = 0
			# 	while curr_pos < median_pos:
			# 		curr_pos = self.clusters[region][i]
			# 		i+=1

			# 	median_info =  self.clusterinfo[region][i]				

			self.final_bps.append( (region[1], 'GetDistr', 'FCD', start_pos, end_pos, avg_region_pvalue,'.','.','type:{0};avg_nr_span_obs:{1};mean_obs_isize:{2}'.format(region[0], avg_region_nr_obs, avg_region_mean_obs) ) )
			#self.final_bps.append( (region[1], 'GetDistr', 'FCD', start_pos, end_pos, median_info[0],'.','.','type:{0};avg_nr_span_obs:{1};mean_obs_isize:{2}'.format(region[0], median_info[1], median_info[2]) ) )
			

	def __str__(self):
		"""
			Prints GFF/GTF file of expansion and contraction regions breakpoint regions
		"""
		output_string= 'seqname\tsource\tfeature\tstart\tend\tscore(p_val)\tstrand\tframe\tattribute\n'
		
		for seqname, source, feature, start, end, avg_p_val, strand, frame, attribute in self.final_bps:
			if start > self.param.mean and end < self.param.scaf_lengths[seqname] - self.param.mean:
				output_string += '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(seqname, source, feature, int(start), int(end), avg_p_val, strand, frame, attribute) 
		return output_string

def AdjustInsertsizeDist(mean_insert, std_dev_insert, insert_list):
    k = MaxObsDistr(len(insert_list), 0.95)
    filtered_list = list(filter((lambda x : (x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert)), insert_list))
    if len(insert_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)


def calc_p_values(bamfile,outfile,param, info_file,assembly_dict):

	p_values = []
	with pysam.Samfile(bamfile, 'rb') as bam:

		#sample true distribution
		param.sample_distribution(bam, info_file)

		# start scoring
		#reference_tids = map(lambda x: bam.gettid(x),bam.references )
		reference_lengths = map(lambda x: int(x), bam.lengths)
		scaf_dict = dict(zip(bam.references, reference_lengths))
		bam_filtered = ifilter(lambda r: r.flag <= 200, bam)
		current_scaf = -1
		print >> info_file, scaf_dict
		param.scaf_lengths = scaf_dict
		for i,read in enumerate(bam_filtered):
			if read.is_unmapped:
				continue
			else:
				current_coord = read.pos
				current_ref = bam.getrname(read.tid)


			if (i + 1) %100000 == 0:
				# print i
				print >> info_file, 'Processing coord:{0}'.format(current_coord)
			# if (i+1) % 30000 == 0:
			# # # 	print 'extra!'
			# 	break

			# initialize read container for new scaffold
			if current_ref != current_scaf:
				print >> info_file, current_ref
				container = []
				scaf_length = scaf_dict[current_ref]

				# store positions in reverse to reduce complexity when removing items from lits in
				# python. We remove the last item and so on
				for i in range(scaf_length,0,-1):
					container.append(ReadContainer(i))
				current_scaf = current_ref 
			# the read pairs we want to use for calculating FCD
			# also, we follow reaprs suggestion and remove read pairs that are further away than
			# 1.5*param.mean from the position of interest. First we can safely remove anything
			# bigger than 3*param.mean (because at least one read of this read pair
			# is going to be further away than 1.5*mean from the position of interest in this read pair )
			if is_proper_aligned_unique_innie(read) and abs(read.tlen) <= 3*param.mean:
				#if current_scaf == 'scf_gap0_errorsize75' and (read.pos > 2352 or read.mpos > 2350):
				#	print current_scaf, read.pos, read.mpos
				if read.aend >= scaf_length or read.aend < 0 or read.mpos +read.rlen > scaf_length or read.pos < 0:
					#print 'Read coordinates outside scaffold length for {0}:'.format(current_scaf), read.aend, read.aend, read.mpos +read.rlen, read.pos 
					continue
	
				# Here we only add the observations that are not
				# further away than the pairs 1.5*param.mean < obs < 3*param.mean 
				if abs(read.tlen) > 1.5*param.mean:
					excluded_region_size = int(abs(read.tlen) - 1.5*param.mean)
					#print 'LOOOOOOL'
				else:
					excluded_region_size = 0

				if read.tlen > 0:
					inner_start_pos = read.aend
					inner_end_pos = read.mpos

					for pos in range(inner_start_pos + excluded_region_size, inner_end_pos - excluded_region_size):
						container[scaf_length - pos].add_read(read)
						#container[pos].add_read(read2)
				else:
					inner_start_pos = read.mpos +read.rlen
					inner_end_pos = read.pos

					for pos in range(inner_start_pos + excluded_region_size, inner_end_pos - excluded_region_size):
						container[scaf_length - pos].add_read(read)
						#container[pos].add_read(read2)

			# write positions out to file
			if  current_coord > scaf_length - len(container):

				while scaf_length - len(container) < current_coord:
					#print 'lloollzz',current_coord, scaf_length - len(container)
					
					#print container[-1].position, current_ref
					# get true distribution
					if container[-1].position % 10000 == 0:
						print 'position', container[-1].position
					sequence_in_window = assembly_dict[ current_ref ][container[-1].position - int(1.5*param.mean) : container[-1].position + int(1.5*param.mean) ]
					p = re.compile("[Nn]+")
					gap_coordinates = []
					for m in p.finditer(sequence_in_window):
						gap_coordinates.append((m.start() - int(1.5*param.mean) ,m.end() - int(1.5*param.mean) ))
						#print m.start(), m.end(), m.group()
						#print gap_coordinates
					true_distribution = param.get_correct_ECDF(outfile, gap_coordinates)
					container[-1].calc_observed_insert()
					KS_statistic, two_side_p_val = container[-1].calc_ks_test(true_distribution) 


					# do ks_2_sample
					#container[-1].calc_observed_insert()
					#KS_statistic, two_side_p_val = container[-1].calc_ks_test(param.true_distr) 
					if two_side_p_val > 0:
						p_values.append(two_side_p_val)
					# write chromosome, coord, p_val to file
					container[-1].write_pval_to_file(outfile,current_ref)
					del container[-1]


	#plt.hist(p_values,bins=50)
	#plt.show()


# def get_misassembly_clusters(container,param):
# 	# tot_mean = 0
# 	# nr_obs =0
# 	# est = model.NormalModel(param.mean,param.stddev,100,s_inner=0)
# 	# # exp_stddev =  lagg till test av forvantad standard avvikelse
# 	# exp_insert = float(est.expected_mean(1,param.scaf_lengths))
# 	# #exp_insert =  mu+ (sigma**2)/float(mu + 1)
# 	# print '#Average predicted mean over a base pair under p_0: {0}'.format(exp_insert)
# 	##
# 	# need to select a consensus loci for a breakpoint
# 	# using only insert size 
# 	sv_container = BreakPointContainer(param)
# 	p_values = []
# 	#pval_threshold = param.get_pval_threshold()
# 	#print "#Adjusted threshold: {0}".format(pval_threshold)
# 	for bp in range(param.scaf_lengths):


# 		#p_val_upper_area = container[bp].calc_pvalue(exp_insert,param.stddev)

# 		# do KS-2sample test
# 		avg_isize = container[bp].calc_observed_insert()
# 		KS_statistic, two_side_p_val = container[bp].calc_ks_test(param.true_distr) 
# 		if two_side_p_val >=0:
# 			p_values.append(two_side_p_val)
# 		else: 
# 			continue

# 		if bp% 1000==0:
# 			print bp, two_side_p_val
# 		# obs_mean = container[bp].calc_observed_insert()
# 		# nr_reads_over_bp = len(container[bp].reads)
# 		# if  nr_reads_over_bp > 10:
# 		# 	nr_obs += 1
# 		# 	tot_mean += obs_mean


# 		if two_side_p_val < param.pval:
# 			if avg_isize > param.adjusted_mean:
# 				sv_container.add_bp_to_cluster(bp,  two_side_p_val, len(container[bp].isize_list), container[bp].mean_isize, 'expansion', param.d)
# 			else:
# 				sv_container.add_bp_to_cluster(bp,  two_side_p_val, len(container[bp].isize_list), container[bp].mean_isize, 'contraction', param.d)

# 	#print p_values
# 	plt.hist(p_values,bins=50)
# 	plt.show()
		
# 	return sv_container


def is_significant(window,pval):
	"""
	The window has a misassembly according to reapers thresholds, that is 
	, at least 80% of the bases in the window has a p-value under pval.
	"""
	significant = []
	for pos, pos_p_val, n_obs, mean, stddev in window:
		if 0 < pos_p_val < pval:
			significant.append(pos_p_val)

	if len(significant)/float(len(window)) >= 0.8:
		return True
	else:
		return False

def get_misassemly_regions(pval_file,param, info_file):
	window = []
	sv_container = BreakPointContainer(param)

	for line in pval_file:
		[scf, pos, pos_p_val, n_obs, mean, stddev] = line.strip().split()
		[pos, pos_p_val, n_obs, mean, stddev] = map(lambda x: float(x), [pos, pos_p_val, n_obs, mean, stddev] )
		window.append( (pos, pos_p_val, n_obs, mean, stddev) )
		if pos >= param.window_size:
			w = window[-param.window_size:]
			if is_significant(w, param.pval):
				avg_pval = sum(map(lambda x: x[1],w))/param.window_size
				avg_obs = sum(map(lambda x: x[2],w))/param.window_size
				avg_mean = sum(map(lambda x: x[3],w))/param.window_size

				if avg_mean > param.adjusted_mean:
					sv_container.add_bp_to_cluster(scf, pos, avg_pval, n_obs, mean, 'expansion', param.window_size)
				else:
					sv_container.add_bp_to_cluster(scf, pos, avg_pval, n_obs, mean, 'contraction', param.window_size)
				#print 'start', len(window) - window_size, 
				#print 'end', scf, pos, pos_p_val, n_obs ,mean, stddev

		if pos % 100000 == 0:
			#print 'Evaluating pos {0}'.format(pos)
			print >> info_file, 'Evaluating pos {0}'.format(pos)

	return sv_container


def main(args):
	param = Parameters()
	if not os.path.exists(args.outfolder):
		os.makedirs(args.outfolder)
	gff_file = open(os.path.join(args.outfolder,'estimated_misassm.gff'),'w')
	info_file = open(os.path.join(args.outfolder,'info.txt'),'w')
	pval_file_out = open(os.path.join(args.outfolder,'p_values.txt'),'w')

	if args.window_size >= 1000:
		param.window_size = args.window_size/2 
	else:
		param.window_size = args.window_size

	param.pval = args.pval

	assembly_dict = ReadInContigseqs(open(args.assembly_file,'r'),param.window_size)
	calc_p_values(args.bampath, pval_file_out, param, info_file,assembly_dict)
	#param.corrected_pval = param.get_pval_threshold()
	pval_file_out.close()
	pval_file_in = open(os.path.join(args.outfolder,'p_values.txt'),'r')
	sv_container =  get_misassemly_regions(pval_file_in,param, info_file) #open(args.pval_file,'r')

	print >> info_file, '#Estimated library params: mean:{0} sigma:{1}'.format(param.mean,param.stddev)
	print >> info_file,'#Genome length:{0}'.format(param.scaf_lengths)

	sv_container.get_final_bp_info(param.window_size)
	print >> gff_file, str(sv_container)
	# output_breaks(sv_container)



if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Infer variants with simple p-value test using theory of GetDistr - proof of concept.")
	parser.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	parser.add_argument('assembly_file', type=str, help='Fasta file with assembly. ')

	parser.add_argument('pval', type=float, help='p-value threshold for calling a variant. ')
	parser.add_argument('window_size', type=int, help='Window size ')
	parser.add_argument('outfolder', type=str, help='Outfolder. ')


	# parser.add_argument('mean', type=int, help='mean insert size. ')
	# parser.add_argument('stddev', type=int, help='Standard deviation of insert size ')


	args = parser.parse_args()
	main(args)



        
