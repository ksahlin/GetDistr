
from collections import deque


def median(l):
    half = len(l) / 2
    l.sort()
    if len(l) % 2 == 0:
        return (l[half-1] + l[half]) / 2.0
    else:
        return l[half]

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
			min_pos = cluster['min']
			max_pos = cluster['max']

			if sv_type_observed != sv_type:
				continue

			if pos <= max_pos + window_size and pos >= min_pos - window_size:
				new_cluster = 0
				cluster['pos'].append(pos)
				if pos < min_pos:
					 cluster['min'] = pos
				if pos > max_pos:
					 cluster['max'] = pos

				self.clusterinfo[(sv_type,scf,i)].append((p_val,nr_obs,mean_obs,window_size))
				break

		if new_cluster:
			self.clusters[(sv_type_observed, scf, self.index)] = {'max':pos,'min':pos,'pos':[pos]}
			self.clusterinfo[(sv_type_observed, scf, self.index)] = [(p_val,nr_obs,mean_obs,window_size)]
			self.index += 1


		if len(self.clusters) == 0:
			self.clusters[(sv_type_observed, scf, self.index)] = {'max':pos,'min':pos,'pos':[pos]}
			self.clusterinfo[(sv_type_observed, scf, self.index)] = [(p_val,nr_obs,mean_obs,window_size)]
			self.index += 1

	def get_final_bp_info(self):
		self.final_bps = []
		for region in self.clusters:

			avg_window_size = sum(map(lambda x:x[-1], self.clusterinfo[region]))/len(self.clusterinfo[region])
			start_pos = self.clusters[region]['min'] #- window_size
			end_pos = self.clusters[region]['max'] + self.clusterinfo[region][-1][-1]

			region_pvalues = map(lambda x:x[0], self.clusterinfo[region])
			median_region_pvalue = median(region_pvalues)  #sum(region_pvalues)/len(region_pvalues)
			region_nr_obs = map(lambda x:x[1], self.clusterinfo[region])
			avg_region_nr_obs = sum(region_nr_obs)/len(region_nr_obs)
			region_mean_obs = map(lambda x:x[2], self.clusterinfo[region])
			avg_region_mean_obs = sum(region_mean_obs)/len(region_mean_obs)

			self.final_bps.append( (region[1], 'GetDistr', 'FCD', start_pos, end_pos, median_region_pvalue,'.','.','type:{0};avg_nr_span_obs:{1};mean_obs_isize:{2};window_size:{3}'.format(region[0], avg_region_nr_obs, avg_region_mean_obs, avg_window_size) ) )
			

	def __str__(self):
		"""
			Prints GFF/GTF file of expansion and contraction regions breakpoint regions
		"""
		output_string= 'seqname\tsource\tfeature\tstart\tend\tscore(p_val)\tstrand\tframe\tattribute\n'
		
		for seqname, source, feature, start, end, avg_p_val, strand, frame, attribute in self.final_bps:
			if start > self.param.mean and end < self.param.scaf_lengths[seqname] - self.param.mean:
				output_string += '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(seqname, source, feature, int(start), int(end), avg_p_val, strand, frame, attribute) 
		return output_string


class Window(object):
	"""docstring for Window"""
	def __init__(self,pval_thresh, max_size, read_length):
		super(Window, self).__init__()
		self.queue = deque() 
		self.avg_inner_mean = None
		self.read_length = read_length
		self.nr_significant = 0
		self.nr_in_window = 0
		self.pval_thresh = pval_thresh
		self.avg_pval = 0
		self.max_window_size = max_size

	def update(self, pos, pval,mean_isize):
		mean_inner_isize = mean_isize - 2*self.read_length
		# initialize window

		if self.avg_inner_mean == None and mean_inner_isize > 0:
			self.avg_inner_mean = mean_inner_isize 
			self.avg_pval = pval
			self.nr_in_window = 1
			if 0 <= pval < self.pval_thresh: 
				self.nr_significant = 1 
			else:
				self.nr_significant = 0
			self.queue.append((pos,pval,mean_inner_isize))

		# update with new value

		if 0 <= pval < self.pval_thresh: 
			self.nr_significant += 1 
		self.avg_inner_mean = (self.nr_in_window * self.avg_inner_mean + mean_inner_isize)/ float(self.nr_in_window+1)
		self.avg_pval = (self.nr_in_window * self.avg_pval + pval)/ float(self.nr_in_window+1)
		self.queue.append((pos,pval,mean_inner_isize))
		self.nr_in_window += 1

		# window is full
		# usually one position if any become significant, but many positions can become
		# significant at once if the average mena ofer the positions drastically decreases.
		# the window size natually decreases and adapts.
		while  self.nr_in_window >= self.avg_inner_mean or  self.nr_in_window >= self.max_window_size:
			if self.is_significant():
				pos, pval_left, mean = self.queue.popleft()

				self.avg_inner_mean = (self.nr_in_window * self.avg_inner_mean - mean)/ float(self.nr_in_window-1)
				self.avg_pval = max((self.nr_in_window * self.avg_pval - pval_left)/ float(self.nr_in_window-1),0)
				if self.avg_pval < 0:
					print 'Negative p-val', pos, self.nr_in_window,self.avg_inner_mean, pval_left, self.avg_pval
				self.nr_in_window -= 1
				if 0 <= pval_left < self.pval_thresh:
					self.nr_significant -=1

				yield pos
			else:
				pos, pval_left, mean = self.queue.popleft()

				self.avg_inner_mean = (self.nr_in_window * self.avg_inner_mean - mean)/ float(self.nr_in_window-1)
				self.avg_pval = max((self.nr_in_window * self.avg_pval - pval_left)/ float(self.nr_in_window-1),0)
				if self.avg_pval < 0:
					print 'Negative p-val non_sign', pos, self.nr_in_window,self.avg_inner_mean, pval_left, self.avg_pval
				self.nr_in_window -= 1
				if 0 <= pval_left < self.pval_thresh:
					self.nr_significant -=1		



	def is_significant(self):
		"""
		The window has a misassembly according to reapers thresholds, that is 
		, at least 80% of the bases in the window has a p-value under pval. Outfolder
		window size is defined as min(reaprs' window size, avg fragment length in window)
		This is because a contraction could not be found otherwise since it will reduce 
		the all insert sizes spanning over.
		"""
		# window is large, need to drop leftmost position and 
		# check if window is significant
	
		if self.nr_significant /float(self.nr_in_window) >= 0.8:
			return True
		else:
			return False




def read_in_gaps(gap_file_path):
	pass

def main(bp_file_path, gap_file_path,outfile_path, param):
	gff_file =  open(outfile_path,'w')

	gap_coordinates =  read_in_gaps(gap_file_path)
	
	current_seq = -1
	sv_container = BreakPointContainer(param)

	for line in open(bp_file_path,'r').readlines():
		scf, pos, n_obs, mean, stddev = line.strip().split()
		if float(n_obs) < 2 or float(mean) < 0:
			current_seq = -1
			continue

		if (scf != current_seq and pos >= param.max_window_size):
			current_seq = scf
			window = Window(param.pval, param.max_window_size, param.read_length)
			window.update(int(pos),float(pos_p_val), float(mean))
		
		else:
			for significant_position in  window.update(int(pos),float(pos_p_val), float(mean)):
				avg_window_mean = window.avg_inner_mean + 2*window.read_length
				if avg_window_mean > param.adjusted_mean:
					sv_container.add_bp_to_cluster(current_seq, int(significant_position), window.avg_pval, int(n_obs), avg_window_mean, 'expansion', min(window.nr_in_window, window.max_window_size))
				else:
					sv_container.add_bp_to_cluster(current_seq, int(significant_position), window.avg_pval, int(n_obs), avg_window_mean, 'contraction', min(window.nr_in_window, window.max_window_size))



	sv_container =  get_misassemly_regions(pval_file_in, param, info_file) #open(args.pval_file,'r')
	sv_container.get_final_bp_info()
	print >> gff_file, str(sv_container)
