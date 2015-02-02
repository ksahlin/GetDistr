from itertools import ifilter
from mathstats.normaldist.normal import MaxObsDistr
from statsmodels.distributions.empirical_distribution import ECDF
import bisect
import random
import pysam

EMPIRICAL_BINS = 500
SAMPLE_SIZE = 50000  # for estimating true full read pair distribution


def AdjustInsertsizeDist(mean_insert, std_dev_insert, insert_list):
    k = MaxObsDistr(len(insert_list), 0.95)
    filtered_list = list(filter((lambda x : (x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert)), insert_list))
    if len(insert_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)


def is_proper_aligned_unique_innie(read):
	return not read.is_unmapped and \
				(read.is_reverse and not read.mate_is_reverse and read.is_read1 and read.tlen < 0 ) or \
	            (not read.is_reverse and read.mate_is_reverse and read.is_read1 and read.tlen > 0 ) \
				or (read.is_reverse and not read.mate_is_reverse and read.is_read2 and read.tlen < 0 ) or \
				(not read.is_reverse and read.mate_is_reverse and read.is_read2 and read.tlen > 0 ) \
	            and not read.mate_is_unmapped and read.mapq > 10 and not read.is_secondary and read.rname == read.mrnm 

class LibrarySampler(object):
	"""docstring for LibrarySampler"""
	def __init__(self, bampath,outpath):
		super(LibrarySampler, self).__init__()
		self.bamfile = pysam.Samfile(bampath, 'rb')
		self.outfile = open(outpath, 'w')
		self.sample_distribution()
		self.bamfile.close()

	def sample_distribution(self):
		isize_list = []
		#i = 0
		bam_filtered = ifilter(lambda r: is_proper_aligned_unique_innie(r), self.bamfile)
		read_lengths = []
		for sample_nr,read in enumerate(bam_filtered):
	   		## add do insert size distribution calculation if proper pair
			if is_proper_aligned_unique_innie(read):
				isize_list.append(abs(read.tlen))
				read_lengths.append(read.rlen)

			if sample_nr > SAMPLE_SIZE:
				break
		print >> self.outfile, '#Insert size sample size:', sample_nr
		#bamfile.reset()

		self.read_length = sum(read_lengths)/float(len(read_lengths))
		isize_list = filter(lambda x: 0 < x - 2*self.read_length,isize_list)
		n_isize = float(len(isize_list))
		mean_isize = sum(isize_list)/n_isize
		std_dev_isize =  (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), isize_list))) / (n_isize - 1)) ** 0.5
		print >> self.outfile,'#Mean before filtering :', mean_isize
		print >> self.outfile,'#Stddev before filtering: ', std_dev_isize
		extreme_obs_occur = True
		while extreme_obs_occur:
			extreme_obs_occur, filtered_list = AdjustInsertsizeDist(mean_isize, std_dev_isize, isize_list)
			n_isize = float(len(filtered_list))
			mean_isize = sum(filtered_list) / n_isize
			std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n_isize - 1)) ** 0.5
			isize_list = filtered_list

		self.min_isize, self.max_isize = min(isize_list), max(isize_list) 
		print >> self.outfile,'#Mean converged:', mean_isize
		print >> self.outfile,'#Std_est converged: ', std_dev_isize
		print >> self.outfile,'{0}\t{1}'.format( mean_isize, std_dev_isize)
		print >> self.outfile,'{0}\t{1}'.format(self.min_isize, self.max_isize )
		print >> self.outfile,'{0}'.format(self.read_length)

		self.nobs = n_isize
		self.mean = mean_isize
		self.stddev = std_dev_isize 
		self.full_ECDF = ECDF(isize_list)
		self.adjustedECDF_no_gap = None
		self.adjustedECDF_no_gap = self.get_correct_ECDF([])
		print >> self.outfile,'#Corrected mean:{0}, corrected stddev:{1}'.format(self.adjusted_mean, self.adjusted_stddev)
		print >> self.outfile,'{0}\t{1}'.format(self.adjusted_mean, self.adjusted_stddev)

		def get_weight(self,x,r,s):
			return x - (2*(r-s)-1)

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

			self.true_distr = [ bisect.bisect(cdf_list_normalized, random.uniform(0, 1)) * stepsize  + x_min for i in range(1000) ]

			# initialization of no gap true distribution
			self.adjustedECDF_no_gap = self.true_distr

			n = len(self.true_distr)
			self.adjusted_mean = sum(self.true_distr)/float(len(self.true_distr))
			self.adjusted_stddev = (sum(list(map((lambda x: x ** 2 - 2 * x * self.adjusted_mean + self.adjusted_mean ** 2), self.true_distr))) / (n - 1)) ** 0.5

		def plot(self):
			pass

