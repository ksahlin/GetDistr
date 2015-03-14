from itertools import ifilter
#from mathstats.normaldist.normal import MaxObsDistr
import pysam
import heapq
import os,sys
import argparse

try:
	import matplotlib
	matplotlib.use('agg')
	import matplotlib.pyplot as plt
	import seaborn as sns
	sns.set_palette("husl", desat=.6)
except ImportError:
	pass




def is_proper_aligned_unique_innie(read):
    mapped_ok = not (read.is_unmapped or read.mate_is_unmapped)
    orientation_ok = (read.is_reverse and not read.mate_is_reverse and read.tlen < 0) or \
            (not read.is_reverse and read.mate_is_reverse and read.tlen > 0)
    quality_ok = read.mapq > 10 and not read.is_secondary
    same_ref = read.rname == read.mrnm

    return mapped_ok and orientation_ok and quality_ok and same_ref


def read_pair_generator(bam,max_isize,param):
	read_pairs = {}
	read_pair_heap = []
	#visited = set()
	prev_read_ref = None
	bam_filtered = ifilter(lambda r: r.flag <= 255, bam)
	for read in bam_filtered:
		
		param.nr_reads += 1
		if not read.is_unmapped:
			param.nr_mapped += 1

		if read.tid != prev_read_ref  and prev_read_ref != None:
			while True:
				try:
					min_pos,r1,mate_pos = heapq.heappop(read_pair_heap)
					yield r1, mate_pos
				except IndexError:
					break
		prev_read_ref = read.tid


		if is_proper_aligned_unique_innie(read) and 0 <= read.tlen <= max_isize and not read.is_reverse:
			if (read.qname, read.is_reverse) in read_pairs:
				print 'bug, multiple alignments', read.qname
				# if '361218' == read.qname:
				# 	print 'lollong here'
				del read_pairs[(read.qname, read.is_reverse)]
				continue
			else:
				read_pairs[(read.qname, read.is_reverse)] = read
				# if '361218' == read.qname:
				# 	print 'pushing here'

		elif is_proper_aligned_unique_innie(read) and  - max_isize <= read.tlen < 0 and read.is_reverse:
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
					if read.pos - max_isize >= min_pos:
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
			min_pos,r1,mate_pos = heapq.heappop(read_pair_heap)
			yield r1, mate_pos
		except IndexError:
			break


def proper_read_isize(bam, min_isize, max_isize, param):
	current_scaf = -1
	reads_fwd = 0
	reads_rev = 0

	already_sampled = set()
	duplicates = set()

	for i,(read,mpos) in enumerate(read_pair_generator(bam,max_isize,param)):
		if abs(read.tlen) <= min_isize:
			continue

		current_ref = bam.getrname(read.tid)

		if current_ref == -1:
			continue


		# print out bp stats for base pairs that we have passed
		if current_ref != current_scaf:
			reads_fwd = 0
			reads_rev = 0
			already_sampled = set()
			duplicates = set()
			current_scaf = current_ref 

		if read.is_reverse: #read.qname in visited and read.is_reverse:
			assert read.tlen < 0
			if (read.qname,mpos,read.pos) in duplicates:
				continue
			reads_rev +=1
			yield read, mpos
			
		else: 

			if (read.pos,mpos) in already_sampled:
				duplicates.add((read.qname,read.pos,mpos))
				continue

			already_sampled.add((read.pos,mpos))
			yield read,mpos
			reads_fwd += 1

class CorrelatedSample(object):
	"""docstring for CorrelatedSample"""
	def __init__(self,n):
		super(CorrelatedSample, self).__init__()
		self.n = n
		self.sample_community_heap = []
		self.sample_positions = set() 

	def remove_inactive(self, read_pos):
		while True:
			try:
				(pos1,pos2) = heapq.heappop(self.sample_community_heap)
			except IndexError:
				break
			if read_pos - pos1 < self.n:
				heapq.heappush(self.sample_community_heap, (pos1,  pos2))
				break
			self.sample_positions.remove((pos1,pos2))


	def add_sample(self, pos1, pos2):
		heapq.heappush(self.sample_community_heap, (pos1,  pos2))
		self.sample_positions.add((pos1,pos2))
		
	def is_correlated(self, mate_pos):
		#sort on smallest mate position
		for pos1,pos2 in heapq.nsmallest(len(self.sample_positions), self.sample_community_heap, key=lambda x: x[1]):
			(pos1,pos2) = heapq.heappop(self.sample_community_heap)
			if abs(mate_pos - pos2) < self.n:
				return True
		
		return False


def within_reference(bampath, outpath,n, param):
	"""
		Assumes that reads are aligned in PE orientation

	"""

	min_isize, max_isize = param.lib_min, param.lib_max

	correlated_check = CorrelatedSample(n)
	bamfile = pysam.Samfile(bampath, 'rb')
	outfile = pysam.Samfile(outpath, 'wb', template=bamfile)

	i = 0
	reads_rev = 0
	reads_fwd = 0
	read_read1 = 0
	read_read2 = 0
	printed_fwd = set()

	for read, mate_pos in proper_read_isize(bamfile, min_isize, max_isize, param):
		read_pos = read.pos

		# remove unactive samples (read1 is more than n bp away)
		correlated_check.remove_inactive(read_pos)

		if not read.is_reverse:
			if correlated_check.is_correlated(mate_pos):
				#print 'correlated'
				continue
				# read is forward
			else:
				# add to sample neighborhood
				correlated_check.add_sample(read_pos, mate_pos)
				# print to filtered bamfile
				outfile.write(read)
				printed_fwd.add(read.qname)
				reads_fwd += 1
				if read.is_read1:
					read_read1 += 1
				elif not read.is_read1:
					read_read2 += 1

		# read is reverse, both reads parsed
		else:
			if read.qname in printed_fwd:
				outfile.write(read)
				printed_fwd.remove(read.qname)
				reads_rev += 1

				if read.is_read1:
					read_read1 += 1
				elif not read.is_read1:
					read_read2 += 1
			#print 'not!'
		i += 1
		if i % 10000 == 1:
			print 'processing coordinate', read.pos, 'on ref:', read.tid
			print reads_fwd,reads_rev,read_read1,read_read2
	print reads_fwd,reads_rev,read_read1,read_read2
	pysam.index(outpath)
	outfile.close()



def between_reference(bampath, outpath,n, min_isize, max_isize, param):
	pass

if __name__ == '__main__':
	# create the top-level parser
	parser = argparse.ArgumentParser()#prog="Infer variants with simple p-value test using theory of GetDistr - proof of concept.")
	#parser.add_argument('--foo', action='store_true', help='help for foo arg.')
	subparsers = parser.add_subparsers(help='help for subcommand')

	# create the parser for the "filter" command	
	filter_parser_within = subparsers.add_parser('filter_within', help='Filters bam file for paired reads within references for better uniform coverage.')
	filter_parser_within.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	filter_parser_within.add_argument('outfolder', type=str, help='Outfolder. ')
	filter_parser_within.add_argument('--n', dest='n', type=int, default=1, help='Neighborhood size. ')	
	filter_parser_within.add_argument('--lib_min', dest='lib_min', type=int, default=200, help='Minimum insert size (if in doubt, just set lib_min = 2*read_length). ')	
	filter_parser_within.add_argument('--lib_max', dest='lib_max', type=int, default=200000, help='Maximum insert size (tight bound is not necessary, choose a larger value rather than smaller). ')	
	filter_parser_within.add_argument('--plots', dest="plots", action='store_true', help='Plot pval distribution.')
	filter_parser_within.set_defaults(which='filter_within')

	# create the parser for the "filter" command	
	filter_parser_between = subparsers.add_parser('filter_between', help='Filters bam file for paired reads on different references for better uniform coverage.')
	filter_parser_between.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	filter_parser_between.add_argument('outfolder', type=str, help='Outfolder. ')
	filter_parser_between.add_argument('--n', dest='n', type=int, default=1, help='Neighborhood size. ')	
	filter_parser_between.add_argument('--lib_max', dest='lib_max', type=int, default=20000, help='Maximum insert size (tight bound is not necessary, choose a larger value rather than smaller). ')	
	filter_parser_between.set_defaults(which='filter_between')
	
	args = parser.parse_args()

	if args.which == 'filter_within' or args.which == 'filter_between':
		try:
		    open(args.bampath)
		except IOError as e:
		    sys.exit("couldn't find BAM file: " + args.bampath + " check that the path is correct and that the file exists")
		try:
		    open(args.bampath + '.bai')
		except IOError as e:
		    print "couldn't find index file: ", args.bampath + '.bai', " check that the path is correct and that the bam file is sorted and indexed"
		    sys.exit(0)

	if not os.path.exists(args.outfolder):
		os.makedirs(args.outfolder)

	if args.plots:
		args.plotfolder = os.path.join(args.outfolder,'plots')

		if not os.path.exists(args.plotfolder):
			os.makedirs(args.plotfolder)

	if args.which == 'filter_between':
		outfile = os.path.join(args.outfolder,'bam_filtered.bam')
		between_reference(args.bampath, outfile, args.n,args.lib_min,args.lib_max)
	elif args.which == 'filter_within':
		outfile = os.path.join(args.outfolder,'bam_filtered.bam')
		within_reference(args.bampath, outfile, args.n,args.lib_min,args.lib_max,args)
	else:
		print 'invalid call'




