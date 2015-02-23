from itertools import ifilter
from mathstats.normaldist.normal import MaxObsDistr
from statsmodels.distributions.empirical_distribution import ECDF
import bisect
import random
import pysam
import heapq
import numpy

EMPIRICAL_BINS = 500
SAMPLE_SIZE = 2**32  # for estimating true full read pair distribution


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
	def __init__(self, bampath,outpath,param):
		super(LibrarySampler, self).__init__()
		self.bamfile = pysam.Samfile(bampath, 'rb')
		self.bampath = bampath
		self.outfile = open(outpath, 'w')
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
			if is_proper_aligned_unique_innie(read) and read.is_read1:
				read_lengths.append(read.rlen)	
				isize_list.append(read.tlen)
				# if abs(read.tlen) > max_tlen:
				# 	max_tlen = abs(read.tlen)
			if sample_nr >= SAMPLE_SIZE:
				break
		self.bamfile.reset()
		#max_tlen = max_tlen+1000
		self.read_length = sum(read_lengths)/float(len(read_lengths))

		## sample proper reads
		
		# isize_list = []
		# for sample_nr,read in enumerate(proper_read_isize_iter(self.bampath, self.read_length, max_tlen)):
	 #   		isize_list.append(read)
		# 	if sample_nr > SAMPLE_SIZE:
		# 		break
		print >> self.outfile, '#Insert size sample size:', sample_nr
		
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
		self.adjustedECDF_no_gap = self.get_correct_ECDF()
		print >> self.outfile,'#Corrected mean:{0}, corrected stddev:{1}'.format(self.adjusted_mean, self.adjusted_stddev)
		print >> self.outfile,'{0}\t{1}'.format(self.adjusted_mean, self.adjusted_stddev)

		samples = min(SAMPLE_SIZE,len(isize_list))
		ess = self.effectiveSampleSize(isize_list[:samples] )
		print 'ESS:', ess
		self.ess_ratio = ess / float(samples)
		print >> self.outfile,'{0}'.format(self.ess_ratio)
		reference_lengths = map(lambda x: int(x), self.bamfile.lengths)
		ref_list = zip(self.bamfile.references, reference_lengths)
		total_base_pairs = sum(reference_lengths)
		print >> self.outfile,'{0}'.format(total_base_pairs)
		for ref, length in ref_list:
			print >> self.outfile,'{0}\t{1}'.format(ref, length)


	def get_weight(self,x,r,s):
		return x - (2*(r-s)-1)

	def plot(self):
		pass

	def effectiveSampleSize(self, data, stepSize = 1) :
		""" Effective sample size, as computed by BEAST Tracer.

		:param data: sequence of real values
		"""
		import subprocess
		proc = subprocess.Popen(['/path/to/RScript','/path/to/plottingfile.R'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = proc.communicate()
	  
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

	  return ess


# def read_pair_generator(bam,max_isize):
# 	read_pairs = {}
# 	read_pair_heap = []
# 	#visited = set()
# 	prev_read_ref = None
# 	bam_filtered = ifilter(lambda r: r.flag <= 255, bam)
# 	for read in bam_filtered:
# 		if read.tid != prev_read_ref  and prev_read_ref != None:
# 			while True:
# 				try:
# 					min_pos,r1,mate_pos = heapq.heappop(read_pair_heap)
# 					yield r1, mate_pos
# 				except IndexError:
# 					break
# 		prev_read_ref = read.tid


# 		if is_proper_aligned_unique_innie(read) and 0 <= read.tlen <= max_isize and not read.is_reverse:
# 			if (read.qname, read.is_reverse) in read_pairs:
# 				print 'bug, multiple alignments', read.qname
# 				# if '361218' == read.qname:
# 				# 	print 'lollong here'
# 				del read_pairs[(read.qname, read.is_reverse)]
# 				continue
# 			else:
# 				read_pairs[(read.qname, read.is_reverse)] = read
# 				# if '361218' == read.qname:
# 				# 	print 'pushing here'

# 		elif is_proper_aligned_unique_innie(read) and  - max_isize <= read.tlen < 0 and read.is_reverse:
# 			if (read.qname, read.is_reverse) in read_pairs :
# 				print 'bug, multiple reverse alignments',read.qname
# 				del read_pairs[(read.qname, read.is_reverse)]
# 				continue

# 			elif (read.qname, not read.is_reverse) in read_pairs:
# 				read_pairs[(read.qname, read.is_reverse)] = read
# 				#print 'gg',read.qname
# 				#if '361218' in read_pairs:
# 				#	print 'lollzzz'
# 				#visited.add(read.qname)
# 				read1 = read_pairs[(read.qname, not read.is_reverse)]	
# 				if read.tid != read1.tid:
# 					del read_pairs[(read.qname, not read.is_reverse)]
# 					del read_pairs[(read.qname, read.is_reverse)]
# 					continue
# 				assert read.mpos == read1.pos
# 				assert read.pos == read1.mpos
# 				# print 'Read has another forward alignment'
# 				# print read.pos, read.is_secondary, read.is_reverse
# 				# print read1.pos, read1.is_secondary, read1.is_reverse
# 				heapq.heappush(read_pair_heap, (read1.pos, read1, read.pos))
# 				heapq.heappush(read_pair_heap, (read.pos, read, read1.pos))
# 				while True:
# 					try:
# 						min_pos,r,mate_pos = heapq.heappop(read_pair_heap)
# 						#print 'index', r1.qname,r2.qname
# 						# print r1.pos, r2.pos
# 					except IndexError:
# 						print 'NOOO'
# 						break
# 					if read.pos - max_isize >= min_pos:
# 						#print 'p',read1.pos, min_pos
						
# 						#print 'here!', r1.pos,r1.flag,r1.mpos
# 						try:
# 							del read_pairs[(r.qname, r.is_reverse)]
# 							yield r, mate_pos
# 						except KeyError:
# 							pass
# 							print 'gah',read.is_reverse
# 							# print r1.qname, r2.qname
# 							# print r1.pos, r2.pos
							

# 					else:
# 						heapq.heappush(read_pair_heap, (min_pos, r, mate_pos))
# 						break


# 	# last reads
# 	while True:
# 		try:
# 			min_pos,r1,mate_pos = heapq.heappop(read_pair_heap)
# 			yield r1, mate_pos
# 		except IndexError:
# 			break


# def proper_read_isize_iter(bam_file_path, read_length, max_isize):
# 	with pysam.Samfile(bam_file_path, 'rb') as bam:
# 		reference_lengths = dict(zip(bam.references, map(lambda x: int(x), bam.lengths)))
# 		current_scaf = -1

# 		reads_fwd = 0
# 		reads_rev = 0

# 		already_sampled = set()
# 		duplicates = set()

# 		for i,(read,mpos) in enumerate(read_pair_generator(bam,max_isize)):
# 			if abs(read.tlen) <= 2*read_length:
# 				continue

# 			current_ref = bam.getrname(read.tid)

# 			if current_ref == -1:
# 				continue


# 			# print out bp stats for base pairs that we have passed
# 			if current_ref != current_scaf:
# 				reads_fwd = 0
# 				reads_rev = 0
# 				already_sampled = set()
# 				duplicates = set()
# 				current_scaf = current_ref 

# 			if read.is_reverse: #read.qname in visited and read.is_reverse:
# 				assert read.tlen < 0
# 				if (read.qname,mpos,read.pos) in duplicates:
# 					continue
# 				reads_rev +=1
				
# 			else: 

# 				if (read.pos,mpos) in already_sampled:
# 					duplicates.add((read.qname,read.pos,mpos))
# 					continue

# 				already_sampled.add((read.pos,mpos))
# 				yield read.tlen
# 				reads_fwd += 1
