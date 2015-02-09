
from itertools import ifilter
import math
import lib_est
import pysam
import sys
import os
import pickle

from statsmodels.distributions.empirical_distribution import ECDF

import heapq

try:
	import matplotlib
	matplotlib.use('agg')
	import matplotlib.pyplot as plt
except ImportError:
	pass


def plot_bp_specific_distr(infile, param):
	means = {}
	stddevs = {}
	for i in [2, 51, 201,501]:
		means[i]=[]
		stddevs[i] = []

	avg_mean = 0
	avg_stddev = 0
	avg_spancov = 0
	tot_pos = 0

	for line in infile:
		[ref,pos, n_obs,mean,sigma] = line.strip().split()
		n_obs = int(float(n_obs))
		mean = float(mean)
		sigma = float(sigma)
		
		if n_obs > 2:
			avg_mean += mean
			avg_stddev += sigma
			avg_spancov += n_obs
			tot_pos += 1

		if 2 < n_obs <= 50:
			means[2].append(mean)
			stddevs[2].append(sigma)

		elif 51 < n_obs <= 200:
			means[51].append(mean)
			stddevs[51].append(sigma)

		elif 201 < n_obs <= 500:
			means[201].append(mean)
			stddevs[201].append(sigma)

		elif 501 < n_obs:
			means[501].append(mean)
			stddevs[501].append(sigma)

	# print len(m_1), len(m_2), len(m_3),len(m_4)
	avg_mean = avg_mean / float(tot_pos)
	avg_stddev = avg_stddev / float(tot_pos)
	avg_spancov = avg_spancov /float(tot_pos)
	print avg_mean,avg_stddev, avg_spancov

	nr_obs, mu = zip(*filter(lambda x: means[x[0]] , means.iteritems()))
	nr_obs, sigma = zip(*filter(lambda x: stddevs[x[0]] , stddevs.iteritems()))
	nr_obs = list(nr_obs)
	nr_obs.sort()
	labels = []
	for low,high in zip(nr_obs[:-1],nr_obs[1:]):
		labels.append("{0}-{1} obs".format(low,high))
	labels.append(">{0} obs".format(high))
	plt.hist(mu, stacked=True, bins=100, log=True, label=labels)
	plt.ylabel('Frequency (log scale)')
	plt.xlabel('isize mean of mates spanning over position')
	title = "Bp specific mean insert size (avg. over genome = %.2f)" % (avg_mean)
	plt.title(title)
	plt.legend( )
	out = os.path.join(param.plotfolder, 'bp_specific_mean.pdf')
	plt.savefig(out)
	plt.close()

	plt.hist(sigma, stacked=True, bins=100, log=True, label=labels)
	plt.ylabel('Frequency (log scale)')
	plt.xlabel('isize standard deviation of mates spanning over position')
	title  = "Bp specific stddev of insert size (avg. over genome = %.2f)" % (avg_stddev)
	plt.title(title)
	plt.legend( )
	out = os.path.join(param.plotfolder, 'bp_specific_stddev.pdf')
	plt.savefig(out)
	plt.savefig(out)
	plt.close()




class ECDF_hist(object):
	"""docstring for ECDF_hist"""
	def __init__(self):
		super(ECDF_hist, self).__init__()
		self.means = []


	def make_ECDF(self):
		self.means_ECDF = ECDF(self.means)

	def get_quantiles(self):
		self.n = float(len(self.means))
		self.mean = sum(self.means)/self.n
		self.stddev = (sum(list(map((lambda x: x ** 2 - 2 * x * self.mean + self.mean ** 2), self.means))) / (self.n - 1)) ** 0.5
		self.p_value_table = {}
		self.means.sort()
		lowest_mean = int(self.means[0])
		highest_mean = int(self.means[-1])

		for mean in xrange(lowest_mean,highest_mean+1):
			cdf_val = self.means_ECDF(mean)
			self.p_value_table[mean] = cdf_val
		self.means = []
		self.means_ECDF = ECDF([0,0])


class Scanner(object):
	"""docstring for Scanner"""
	def __init__(self,name,outfile):
		super(Scanner, self).__init__()
		self.ref_name = name
		self.position = 0
		self.mu = 0
		self.var = 0
		self.o = 0
		self.o_sq = 0
		self.n_obs = 0.0
		self.outfile = outfile
		self.ecdf = ECDF_hist()

	def write_bp_stats_to_file(self,bp_index):
		#print  '{0}\t{1}\t{2}\t{3}\t{4}'.format(self.ref_name, bp_index, self.n_obs, self.mu, self.var)
		try:
			math.sqrt(self.var)
		except ValueError:
			print 'Negative stddev:'
			print '{0}\t{1}\t{2}\t{3}\t var:{4}'.format(self.ref_name, bp_index, self.n_obs, self.mu, self.var)
			return

		print >> self.outfile, '{0}\t{1}\t{2}\t{3}\t{4}'.format(self.ref_name, bp_index, int(self.n_obs), self.mu, math.sqrt(self.var))

	def add_obs(self, isize):
		self.o += isize
		self.o_sq += isize**2
		self.n_obs += 1 
		if self.n_obs < 2:
			self.mu = 0
			self.var = 0
			return
		self.mu = self.o / self.n_obs  #(self.n_obs * self.mu + isize)/ float(self.n_obs+1)
		#print self.mu,self.o,self.o_sq
		self.var = 1/(self.n_obs -1)* (self.o_sq - 2*self.mu*self.o + self.n_obs*self.mu**2) #(self.n_obs * self.var + isize**2)/ float(self.n_obs+1)

	def remove_obs(self, neg_isize):
		self.o += neg_isize
		self.o_sq -= neg_isize**2
		self.n_obs -= 1 
		if self.n_obs < 2:
			self.mu = 0
			self.var = 0
			return		
		self.mu = self.o / self.n_obs  #(self.n_obs * self.mu + isize)/ float(self.n_obs+1)
		self.var = 1/(self.n_obs -1)* (self.o_sq - 2*self.mu*self.o + self.n_obs*self.mu**2) #(self.n_obs * self.var + isize**2)/ float(self.n_obs+1)

	def update_pos(self,pos):
		if self.position < pos:
			for bp_index in range(self.position,pos):
				self.write_bp_stats_to_file(bp_index)
				self.ecdf.means.append(self.mu)
		self.position = pos

def read_pair_generator(bam,param):
	read_pairs = {}
	read_pair_heap = []
	#visited = set()
	prev_read_ref = None
	bam_filtered = ifilter(lambda r: r.flag <= 255, bam)
	for read in bam_filtered:
		if read.tid != prev_read_ref and not read.is_unmapped  and prev_read_ref != None:
			while True:
				try:
					min_pos,r1,mate_pos = heapq.heappop(read_pair_heap)
					yield r1, mate_pos
				except IndexError:
					break
		prev_read_ref = read.tid


		if lib_est.is_proper_aligned_unique_innie(read) and 0 <= read.tlen <= param.max_isize and not read.is_reverse:
			if (read.qname, read.is_reverse) in read_pairs:
				print 'bug, multiple alignments', read.qname
				del read_pairs[(read.qname, read.is_reverse)]
				continue
			else:
				read_pairs[(read.qname, read.is_reverse)] = read

		elif lib_est.is_proper_aligned_unique_innie(read) and  -param.max_isize <= read.tlen < 0 and read.is_reverse:
			if (read.qname, read.is_reverse) in read_pairs :
				print 'bug, multiple reverse alignments',read.qname
				del read_pairs[(read.qname, read.is_reverse)]
				continue

			elif (read.qname, not read.is_reverse) in read_pairs:
				read_pairs[(read.qname, read.is_reverse)] = read
				#print 'gg',read.qname
				#if '361218' in read_pairs:
				#	print 'lollzzz'
				#visited.add(read.qname)
				read1 = read_pairs[(read.qname, not read.is_reverse)]	
				if read.tid != read1.tid:
					del read_pairs[(read.qname, not read.is_reverse)]
					del read_pairs[(read.qname, read.is_reverse)]
					continue
				assert read.mpos == read1.pos
				assert read.pos == read1.mpos
				# print 'Read has another forward alignment'
				# print read.pos, read.is_secondary, read.is_reverse
				# print read1.pos, read1.is_secondary, read1.is_reverse
				heapq.heappush(read_pair_heap, (read1.pos, read1, read.pos))
				heapq.heappush(read_pair_heap, (read.pos, read, read1.pos))
				while True:
					try:
						min_pos,r,mate_pos = heapq.heappop(read_pair_heap)
						#print 'index', r1.qname,r2.qname
						# print r1.pos, r2.pos
					except IndexError:
						print 'NOOO'
						break
					if read.pos - param.max_isize >= min_pos:
						#print 'p',read1.pos, min_pos
						
						#print 'here!', r1.pos,r1.flag,r1.mpos
						try:
							del read_pairs[(r.qname, r.is_reverse)]
							yield r, mate_pos
						except KeyError:
							pass
							print 'gah',read.is_reverse
							# print r1.qname, r2.qname
							# print r1.pos, r2.pos
							

					else:
						heapq.heappush(read_pair_heap, (min_pos, r, mate_pos))
						break


	# last reads
	while True:
		try:
			min_pos,r1,r2 = heapq.heappop(read_pair_heap)
			yield r1, r2
		except IndexError:
			break


def parse_bam(bam_file,param):



	outfile = open(os.path.join(param.outfolder,'bp_stats.txt'),'w')

	with pysam.Samfile(bam_file, 'rb') as bam:
		reference_lengths = dict(zip(bam.references, map(lambda x: int(x), bam.lengths)))
		current_scaf = -1
		read_len = int(param.read_length)
		counter = 0

		reads_fwd = 0
		reads_rev = 0

		already_sampled = set()
		duplicates = set()
		scanner = Scanner('init',outfile)

		for i,(read,mpos) in enumerate(read_pair_generator(bam,param)):
			#print read.pos, mpos #, read2.pos

			if i %100000 == 0:
				print i

			if abs(read.tlen) <= 2*param.read_length:
				continue


			current_coord = read.pos + read_len
			current_ref = bam.getrname(read.tid)

			if current_ref == -1:
				continue

				
			if reference_lengths[current_ref] < param.max_isize:
				continue

			# print out bp stats for base pairs that we have passed
			if current_ref != current_scaf:
				print current_scaf
				#print 'visited to new ref', len(visited)
				print 'reads_fwd on scanned contig:',reads_fwd
				print 'reads_rev on scanned contig:',reads_rev
				reads_fwd = 0
				reads_rev = 0
				already_sampled = set()
				duplicates = set()

				# pass the emperical distribution on to the next object
				ecdf_ = scanner.ecdf
				scanner = Scanner(current_ref,outfile)
				scanner.ecdf = ecdf_

				scanner.update_pos(current_coord)
				scaf_length = reference_lengths[current_ref]
				current_scaf = current_ref 
			else:
				scanner.update_pos(current_coord)

			if read.is_reverse: #read.qname in visited and read.is_reverse:
				assert read.tlen < 0
				#print 'reads_fwd on scanned contig:',reads_fwd
				#print 'reads_rev on scanned contig:',reads_rev
				if (read.qname,mpos,read.pos) in duplicates:
					continue
				scanner.remove_obs(read.tlen)
				reads_rev +=1
				

			else: #if lib_est.is_proper_aligned_unique_innie(read) and (param.min_isize <= read.tlen <= param.max_isize):

				if read.aend >= scaf_length or read.aend < 0 or read.mpos +read.rlen > scaf_length or read.pos < 0:
					print 'Read coordinates outside scaffold length for {0}:'.format(current_scaf), read.aend, read.aend, read.mpos +read.rlen, read.pos 
					#continue
				# if read.tlen <0 :
				# 	print 'BUG', read.tlen

				if (read.pos,mpos) in already_sampled:
					duplicates.add((read.qname,read.pos,mpos))
					continue

				already_sampled.add((read.pos,mpos))
				scanner.add_obs(read.tlen)
				reads_fwd += 1
				counter += 1
		print 'Good read pair count: ', counter

		scanner.ecdf.make_ECDF()
		scanner.ecdf.get_quantiles()
		pickle.dump(scanner.ecdf, open(os.path.join(param.outfolder,'ecdf.pkl'),'w') )
	outfile.close()
	if param.plots:
		infile = open(os.path.join(param.outfolder,'bp_stats.txt'),'r')
		plot_bp_specific_distr(infile, param)

