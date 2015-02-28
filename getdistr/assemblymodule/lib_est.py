
import os
import math
from itertools import ifilter
import bisect
import random


import pysam
import numpy
from statsmodels.distributions.empirical_distribution import ECDF
from rpy2 import robjects
from rpy2.robjects import packages

#import assemblymodule.filter_bam as fb
from getdistr.assemblymodule import find_normal_parameters as fit
from mathstats.normaldist.normal import MaxObsDistr
from mathstats.normaldist.truncatedskewed import param_est



try:
	import matplotlib
	matplotlib.use('agg')
	import matplotlib.pyplot as plt
	import seaborn as sns
	sns.set_palette("husl", desat=.6)
except ImportError:
	pass


EMPIRICAL_BINS = 500
SAMPLE_SIZE = 2**32  # for estimating true full read pair distribution
NORMAL_QUANTILE_TWO_SIDED_95 = 1.96

def plot_isize(isizes,outfile):
    plt.hist(isizes,bins=100)
    plt.ylabel('frequency') 
    plt.xlabel('fragment size')  
    plt.title('Insert size distribution')
    plt.legend( )
    plt.savefig(outfile)
    plt.close()
    plt.clf()

def AdjustInsertsizeDist(mean_insert, std_dev_insert, insert_list):
    k = MaxObsDistr(len(insert_list), 0.95)
    filtered_list = list(filter((lambda x : (x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert)), insert_list))
    if len(insert_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)


def is_proper_aligned_unique_innie(read):
    mapped_ok = not (read.is_unmapped or read.mate_is_unmapped)
    orientation_ok = (read.is_reverse and not read.mate_is_reverse and read.tlen < 0) or \
            (not read.is_reverse and read.mate_is_reverse and read.tlen > 0)
    quality_ok = read.mapq > 10 and not read.is_secondary
    same_ref = read.rname == read.mrnm

    return mapped_ok and orientation_ok and quality_ok and same_ref

	# return not read.is_unmapped and not read.mate_is_unmapped and \
	# 			( (read.is_reverse and not read.mate_is_reverse  and read.tlen < 0 ) or \
	#             (not read.is_reverse and read.mate_is_reverse  and read.tlen > 0 ) ) \
	#             and  and read.mapq > 10 and not read.is_secondary and read.rname == read.mrnm 


class LibrarySampler(object):
	"""docstring for LibrarySampler"""
	def __init__(self, bampath,param):
		super(LibrarySampler, self).__init__()
		self.bamfile = pysam.Samfile(bampath, 'rb')
		self.bampath = bampath
		self.param = param
		self.lib_file = open(os.path.join(param.outfolder,'library_info.txt'), 'w')
		self.stats_file = open(os.path.join(param.outfolder,'stats.txt'), 'w')


		self.sample_distribution()
		self.bamfile.close()



	def get_correct_ECDF(self):
		read_len = int(self.read_length)
		softclipps = 0 #read_len #int(self.read_length*0.6)


		x_min = self.min_isize #max(2*(read_len-softclipps) , int(self.mean - 5*self.stddev) )
		x_max = self.max_isize #int(self.mean + 5*self.stddev)
		stepsize =  max(1,(x_max - x_min) / EMPIRICAL_BINS)
		#print (x_max - x_min)/float(EMPIRICAL_BINS)
		cdf_list = [ self.full_ECDF( x_min) * self.get_weight(int(round(x_min+stepsize/2.0,0)), read_len, softclipps)  ] #[ self.full_ECDF( 2*(read_len-softclipps)) * self.get_weight(2*(read_len-softclipps), gap_coordinates, read_len, softclipps) ]


		# create weigted (true) distribution

		for x in range( x_min + stepsize , x_max, stepsize):
			increment_area = self.get_weight(int((x+stepsize/2)), read_len, softclipps) * (self.full_ECDF(x) - self.full_ECDF(x-stepsize))
			cdf_list.append( cdf_list[-1] + increment_area)

		tot_cdf = cdf_list[-1]

		cdf_list_normalized = map(lambda x: x /float(tot_cdf),cdf_list)

		# Now create a weighted sample

		self.true_distr = [ bisect.bisect(cdf_list_normalized, random.uniform(0, 1)) * stepsize  + x_min for i in range(50000) ]

		# initialization of no gap true distribution
		self.adjustedECDF_no_gap = self.true_distr

		n = len(self.true_distr)
		self.adjusted_mean = sum(self.true_distr)/float(len(self.true_distr))
		self.adjusted_stddev = (sum(list(map((lambda x: x ** 2 - 2 * x * self.adjusted_mean + self.adjusted_mean ** 2), self.true_distr))) / (n - 1)) ** 0.5


	def sample_distribution(self):
		
		## pilot sample
		read_lengths = []
		# max_tlen = 0
		bam_filtered = ifilter(lambda r: is_proper_aligned_unique_innie(r), self.bamfile)
		isize_list = []
		for sample_nr,read in enumerate(bam_filtered):
	  		## add do insert size distribution calculation if proper pair
			if is_proper_aligned_unique_innie(read) and not  read.is_reverse:
				assert read.tlen > 0
				read_lengths.append(read.rlen)	
				isize_list.append(read.tlen)
				# if abs(read.tlen) > max_tlen:
				# 	max_tlen = abs(read.tlen)
			if sample_nr >= SAMPLE_SIZE:
				break

		# for read, mate_pos in fb.proper_read_isize(self.bamfile, self.param.lib_min, self.param.ligetdistr.assemblymodule.b_max):
	 #  		sample_nr += 1
	 #  		## add do insert size distribution calculation if proper pair
	 #  		if read.tlen >= 0:
		# 	#if is_proper_aligned_unique_innie(read) and read.is_read1:
		# 		read_lengths.append(read.rlen)	
		# 		isize_list.append(read.tlen)
		# 		# if abs(read.tlen) > max_tlen:
		# 		# 	max_tlen = abs(read.tlen)
		# 	if sample_nr >= SAMPLE_SIZE:
		# 		break

		self.bamfile.reset()
		#max_tlen = max_tlen+1000
		self.read_length = sum(read_lengths)/float(len(read_lengths))

		## sample proper reads
		
		# isize_list = []
		# for sample_nr,read in enumerate(proper_read_isize_iter(self.bampath, self.read_length, max_tlen)):
	 #   		isize_list.append(read)
		# 	if sample_nr > SAMPLE_SIZE:
		# 		break
		print >> self.lib_file, '#Insert size sample size:', sample_nr
		
		isize_list = filter(lambda x: 0 < x - 2*self.read_length,isize_list)
		n_isize = float(len(isize_list))
		mean_isize = sum(isize_list)/n_isize
		std_dev_isize =  (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), isize_list))) / (n_isize - 1)) ** 0.5
		print >> self.lib_file,'#Mean before filtering :', mean_isize
		print >> self.lib_file,'#Stddev before filtering: ', std_dev_isize
		extreme_obs_occur = True
		while extreme_obs_occur:
			extreme_obs_occur, filtered_list = AdjustInsertsizeDist(mean_isize, std_dev_isize, isize_list)
			n_isize = float(len(filtered_list))
			mean_isize = sum(filtered_list) / n_isize
			std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n_isize - 1)) ** 0.5
			isize_list = filtered_list

		self.min_isize, self.max_isize = min(isize_list), max(isize_list) 
		print >> self.lib_file,'#Mean converged:', mean_isize
		print >> self.lib_file,'#Std_est converged: ', std_dev_isize
		print >> self.lib_file,'{0}\t{1}'.format( mean_isize, std_dev_isize)
		print >> self.lib_file,'{0}\t{1}'.format(self.min_isize, self.max_isize )
		print >> self.lib_file,'{0}'.format(self.read_length)

		self.nobs = n_isize
		self.mean = mean_isize
		self.stddev = std_dev_isize 
		self.full_ECDF = ECDF(isize_list)
		self.adjustedECDF_no_gap = None
		self.adjustedECDF_no_gap = self.get_correct_ECDF()
		print >> self.lib_file,'#Corrected mean:{0}, corrected stddev:{1}'.format(self.adjusted_mean, self.adjusted_stddev)
		print >> self.lib_file,'{0}\t{1}'.format(self.adjusted_mean, self.adjusted_stddev)

		samples = min(SAMPLE_SIZE,len(isize_list))
		ess = self.effectiveSampleSize(isize_list[:samples])
		self.ess_ratio = ess / float(samples)
		print >> self.lib_file,'{0}'.format(self.ess_ratio)
		reference_lengths = map(lambda x: int(x), self.bamfile.lengths)
		ref_list = zip(self.bamfile.references, reference_lengths)
		total_basepairs = sum(reference_lengths)
		print >> self.lib_file,'{0}'.format(total_basepairs)
		for ref, length in ref_list:
			print >> self.lib_file,'{0}\t{1}'.format(ref, length)



		print >> self.stats_file, 'Proper reads sampled:', samples
		print >> self.stats_file, 'ESS of proper reads sampled:', ess
		print >> self.stats_file, 'ESS_ratio:', self.ess_ratio
		coverage = self.read_length*samples*2/float(total_basepairs)
		print >> self.stats_file, 'Mean coverage proper reads:{0}'.format( coverage )
		inner_span_coverage = coverage * (self.mean -2*self.read_length)/(2*self.read_length)
		print >> self.stats_file, 'Average theoretical inner span coverage:{0}'.format( inner_span_coverage )
		print >> self.stats_file, 'Mean full lib:{0}'.format(self.mean)
		print >> self.stats_file, 'Stddev full lib:{0}'.format(self.stddev)
		print >> self.stats_file, 'Emperically adjusted mean:{0}'.format(self.adjusted_mean)
		print >> self.stats_file, 'Emperically adjusted stddev:{0}'.format(self.adjusted_stddev)
		mu_naive = self.mean + self.stddev**2/float(self.mean - 2*self.read_length+1)
		sigma_naive = math.sqrt(self.stddev**2 - self.stddev**4/(self.mean -2*self.read_length +1)**2 )
		print >> self.stats_file, 'Naive adjusted mean:{0}'.format(mu_naive)
		print >> self.stats_file, 'Naive adjusted stddev:{0}'.format(sigma_naive)
		mu_sophisticated = param_est.mean_given_d(self.mean, self.stddev, self.read_length, total_basepairs, total_basepairs, 0)
		sigma_sophisticated = param_est.stddev_given_d(self.mean, self.stddev, self.read_length, total_basepairs, total_basepairs, 0)
		print >> self.stats_file, 'Sophisticated adjusted mean:{0}'.format(mu_sophisticated)
		print >> self.stats_file, 'Sophisticated adjusted stddev:{0}'.format(sigma_sophisticated)
		theoretical_margin_of_error = NORMAL_QUANTILE_TWO_SIDED_95*self.stddev / math.sqrt(inner_span_coverage)
		print >> self.stats_file, 'Theoretical margin of error two sided 95%:{0}'.format(theoretical_margin_of_error)
		self.stats_file.close()




		if self.param.plots:
			outfile = os.path.join(self.param.plotfolder, 'isize.eps')
			plot_isize(isize_list, outfile)	
			outfile = os.path.join(self.param.plotfolder, 'fitted_params_isize.eps')
			fit.main(isize_list, outfile)


	def get_weight(self,x,r,s):
		return x - (2*(r-s)-1)

	def lowest_bound_biological(self):
		pass

	def plot(self):
		pass

	def effectiveSampleSize(self,data) :
		""" 
			Effective sample size, as computed by EffectiveSize in coda, R.
			returns a python float, the ess
		"""
		robjects.packages.importr("coda")
		r = robjects.r
		data = numpy.array(data)
		normalizedData = data - data.mean()
		rData = robjects.IntVector(normalizedData)
		mcmc_r = r.mcmc(rData)	 
		return list(r.effectiveSize(mcmc_r))[0]

	  # samples = len(data)

	  # assert len(data) > 1,"no stats for short sequences"

	  # maxLag = min(samples//3, 1000)

	  # gammaStat = [0,]*maxLag
	  # #varGammaStat = [0,]*maxLag

	  # varStat = 0.0;

	  # if type(data) != numpy.ndarray :
	  #   data = numpy.array(data)

	  # normalizedData = data - data.mean()
	  
	  # for lag in range(maxLag) :
	  #   v1 = normalizedData[:samples-lag]
	  #   v2 = normalizedData[lag:]
	  #   v = v1 * v2
	  #   gammaStat[lag] = sum(v) / len(v)
	  #   #varGammaStat[lag] = sum(v*v) / len(v)
	  #   #varGammaStat[lag] -= gammaStat[0] ** 2

	  #   # print lag, gammaStat[lag], varGammaStat[lag]
	    
	  #   if lag == 0 :
	  #     varStat = gammaStat[0]
	  #   elif lag % 2 == 0 :
	  #     s = gammaStat[lag-1] + gammaStat[lag]
	  #     if s > 0 :
	  #        varStat += 2.0*s
	  #     else :
	  #       break
	      
	  # # standard error of mean
	  # # stdErrorOfMean = Math.sqrt(varStat/samples);

	  # # auto correlation time
	  # act = stepSize * varStat / gammaStat[0]

	  # # effective sample size
	  # ess = (stepSize * samples) / act
