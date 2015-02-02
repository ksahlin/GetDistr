
from itertools import ifilter

import lib_est

class GenomeWindow(object):
	"""docstring for GenomeWindow"""
	def __init__(self):
		super(GenomeWindow, self).__init__()
		self.start_position = None
		self.mu = None
		self.sigma = None
		self.n_obs = 0

	def write_pval_to_file(self,outfile,ref_name):
		print >> outfile, '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(ref_name, self.position, self.pval,len(self.n_obs), self.mean_isize, self.std_dev_isize)

	def update_mean(self, read):
		pass

	def update_sigma(self, read):
		pass


	def calc_p_val(self,mu):
		if self.n_obs >= 5:
			return self.pval 
		else:
			self.pval = -1
			return -1


def parse_bam(bam,libstats,outfile):

	reference_lengths = map(lambda x: int(x), bam.lengths)
	

	# start scoring
	#reference_tids = map(lambda x: bam.gettid(x),bam.references )

	scaf_dict = dict(zip(bam.references, reference_lengths))
	bam_filtered = ifilter(lambda r: r.flag <= 200, bam)
	current_scaf = -1
	prev_coord = (0,0)
	duplicate_count = 0
	visited = set()
	#observed_span_distribution = []
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


		if scaf_dict[current_ref] < libstats.max_isize:
			continue
			
		# initialize read container for new scaffold
		if current_ref != current_scaf:
			print current_ref
			container = []
			scaf_length = scaf_dict[current_ref]

			# store positions in reverse to reduce complexity when removing items from lits in
			# python. We remove the last item and so on
			for i in range(scaf_length,0,-1):
				container.append(GenomeWindow(i))
			current_scaf = current_ref 
		# the read pairs we want to use for calculating FCD
		# also, we follow reaprs suggestion and remove read pairs that are further away than
		# 1.5*libstats.mean from the position of interest. First we can safely remove anything
		# bigger than 3*libstats.mean (because at least one read of this read pair
		# is going to be further away than 1.5*mean from the position of interest in this read pair )
		if read.qname in visited:
			#print read.tlen
			#print read.qname
			visited.remove(read.qname) # we have visited both in pair
			continue
		if is_proper_aligned_unique_innie(read) and (libstats.min_isize <= abs(read.tlen) <= libstats.max_isize):
			if read.aend >= scaf_length or read.aend < 0 or read.mpos +read.rlen > scaf_length or read.pos < 0:
				print 'Read coordinates outside scaffold length for {0}:'.format(current_scaf), read.aend, read.aend, read.mpos +read.rlen, read.pos 
				continue

			# Here we only add the observations that are not
			# further away than the pairs 1.5*libstats.mean < obs < 3*libstats.mean 

			excluded_region_size = 0

			if read.tlen > 0:
				inner_start_pos = read.aend
				inner_end_pos = read.mpos
				for pos in range(inner_start_pos + excluded_region_size, inner_end_pos - excluded_region_size):
					try:
						container[scaf_length - pos].isize_list.append(read.tlen) # .add_read(read)

					#container[pos].add_read(read2)
					except IndexError:
						print pos, read.tlen
						break
						pass
						#print read.tlen
						#print 'Tried adding read pair to position {0}. On scaffold {1} with length {2}, with container of size {3}'.format(pos, current_scaf, scaf_length, len(container)) 
			else:

				inner_start_pos = read.mpos +read.rlen
				inner_end_pos = read.pos
				for pos in range(inner_start_pos + excluded_region_size, inner_end_pos - excluded_region_size):
					try:
						container[scaf_length - pos].isize_list.append(abs(read.tlen)) #add_read(read)
					except IndexError:
						print pos, read.tlen
						break

					#container[pos].add_read(read2)
			visited.add(read.qname)

