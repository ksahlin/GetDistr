
from itertools import ifilter
import math
import lib_est
import pysam

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

	def write_bp_stats_to_file(self,bp_index):
		#print  '{0}\t{1}\t{2}\t{3}\t{4}'.format(self.ref_name, bp_index, self.n_obs, self.mu, self.var)
		print >> self.outfile, '{0}\t{1}\t{2}\t{3}\t{4}'.format(self.ref_name, bp_index, self.n_obs, self.mu, math.sqrt(self.var))

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

	def remove_obs(self, isize):
		self.o -= isize
		self.o_sq -= isize**2
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
		self.position = pos


def parse_bam(bam_file,libstats,out_path):
	outfile = open(out_path,'w')
	print 'here'
	with pysam.Samfile(bam_file, 'rb') as bam:
		reference_lengths = dict(zip(bam.references, map(lambda x: int(x), bam.lengths)))
		bam_filtered = ifilter(lambda r: r.flag <= 200, bam)
		current_scaf = -1
		prev_coord = (0,0)
		duplicate_count = 0
		visited = set()

		for i,read in enumerate(bam_filtered):

			current_coord = read.pos
			current_ref = bam.getrname(read.tid)

			coord1 = read.pos
			coord2 = read.mpos

			if (coord1, coord2) == prev_coord:
				duplicate_count += 1
				continue
			else:
				prev_coord = (coord1, coord2)


			if (i + 1) %100000 == 0:
				# print i
				print '#Processing read:{0}'.format(current_coord)
				
			if (i+1) % 100000 == 0:
				print '#removed {0} duplicates, processed {1} reads.'.format(duplicate_count,i)
		
			if reference_lengths[current_ref] < libstats.max_isize:
				continue

			# print out bp stats for base pairs that we have passed
			if current_ref != current_scaf:
				print current_ref
				scanner = Scanner(current_ref,outfile)
				# scanner.update_pos(reference_lengths[current_scaf])
				scaf_length = reference_lengths[current_ref]
				current_scaf = current_ref 
			else:
				scanner.update_pos(current_coord)


			# the read pairs we want to use for calculating FCD
			# also, we follow reaprs suggestion and remove read pairs that are further away than
			# 1.5*libstats.mean from the position of interest. First we can safely remove anything
			# bigger than 3*libstats.mean (because at least one read of this read pair
			# is going to be further away than 1.5*mean from the position of interest in this read pair )
			if read.qname in visited:
				visited.remove(read.qname) # we have visited both in pair
				scanner.remove_obs(abs(read.tlen))
				continue
			elif lib_est.is_proper_aligned_unique_innie(read) and (libstats.min_isize <= abs(read.tlen) <= libstats.max_isize):
				if read.aend >= scaf_length or read.aend < 0 or read.mpos +read.rlen > scaf_length or read.pos < 0:
					print 'Read coordinates outside scaffold length for {0}:'.format(current_scaf), read.aend, read.aend, read.mpos +read.rlen, read.pos 
					continue

				if read.tlen > 0:
					#inner_start_pos = read.aend
					#inner_end_pos = read.mpos
					scanner.add_obs(abs(read.tlen))
				else:
					print 'bug?', read.tlen, read.pos, read.mpos, read.qname
					scanner.add_obs(abs(read.tlen))
					#inner_start_pos = read.mpos +read.rlen
					#inner_end_pos = read.pos

				visited.add(read.qname)

