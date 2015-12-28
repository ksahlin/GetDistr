
from collections import deque
import math
from scipy.stats import t,norm
import os
import pickle

try:
	import matplotlib
	matplotlib.use('agg')
	import matplotlib.pyplot as plt
	import seaborn as sns
	sns.set_palette("husl", desat=.6)
	sns.set(font_scale=3)
	plt.rcParams.update({'font.size': 22})
except:
	pass

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
			if start > self.param.mu and end < self.param.scaf_lengths[seqname] - self.param.mu:
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

def erfcc(x):
    """Complementary error function."""
    z = abs(x)
    k = 1. / (1. + 0.5 * z)
    r = k * math.exp(-z * z - 1.26551223 + k * (1.00002368 + k * (.37409196 +
        k * (.09678418 + k * (-.18628806 + k * (.27886807 +
        k * (-1.13520398 + k * (1.48851587 + k * (-.82215223 +
        k * .17087277)))))))))
    if (x >= 0.):
        return r
    else:
        return 2. - r

def normcdf(x, mu, sigma):
    k = x - mu;
    y = 0.5 * erfcc(-k / (sigma * math.sqrt(2.0)));
    if y > 1.0:
        y = 1.0;
    return y


def calc_cdf(mu_null, n, mean, stddev):
	one_sample_t = ( mean - mu_null)/ (stddev/math.sqrt(n))
	if n > 50:
		cdf_val = normcdf(one_sample_t,0. , 1.)
	else:
		cdf_val = t.cdf(one_sample_t, n - 1)
	
	return cdf_val	

def calc_pvalue(cdf_val):
	return 2*min( cdf_val, 1 - cdf_val)



def read_in_gaps(gap_file_path):
	gap_coordinates = {}
	for line in open(gap_file_path, 'r'):
		[scf, start, stop] = line.strip().split()
		gap_coordinates[scf] = (int(start), int(stop))
	return gap_coordinates

def plot_stats(outfile,values, title='',x='',y=''):
	"""
		plots an iterable of values, in this function it''s either p-values
		or CDF values
	"""
	plt.hist(values, bins=200)
	plt.ylabel(x)
	plt.xlabel(y)
	title_ = title
	plt.title(title_)
	plt.legend( )
	plt.savefig(outfile + '.eps')
	plt.close()
	plt.clf()

def plot_mean_chain(param,mean_isizes):
	isize_plot = os.path.join(param.plotfolder,'isize_chain.eps')
	x_axis = map(lambda x: x*500 ,range(len(mean_isizes)))
	plt.plot(x_axis, mean_isizes, '-')
	plt.ylabel('mean isize')
	plt.xlabel('reference position')
	title = "Mean spanning isize"
	plt.title(title)
	plt.legend( )
	plt.savefig(isize_plot)
	plt.close()
	plt.clf()



def main(bp_file_path, gap_file_path, param):
	if param.plots == True:
		good_pos = False # position where to start the isize chain plotting window
		mean_isize = []
		p_values_naive = []
		p_values_correct_bias1 = []
		p_values_correct_bias2 = []
		p_values_correct_both = []
		cdf_values_naive = []
		cdf_values_correct_bias1 = []
		cdf_values_correct_bias2 = []
		cdf_values_correct_both = []

	significant_regions = os.path.join(param.outfolder,'regions.gff')
	ecdf = pickle.load(open(os.path.join(param.outfolder,'ecdf.pkl'),'r'))
	#print ecdf.p_value_table
	gff_file =  open(significant_regions,'w')
	gap_coordinates =  read_in_gaps(gap_file_path)
	current_seq = -1
	sv_container = BreakPointContainer(param)

	print 'ESS_ratio used in plots:',param.ess_ratio
	
	for i, line in enumerate(open(bp_file_path,'r')):
		scf, pos, n_obs, mean, stddev = line.strip().split()
		pos, n_obs,mean,stddev = int(pos), float(n_obs), float(mean), float(stddev)
		# do not consider the very ends of a rererence or positions
		# with less than 2 effective spanning observations
		if (pos < param.adjusted_mu +3*param.adjusted_sigma) or (pos > param.scaf_lengths[scf] - param.adjusted_mu - 3*param.adjusted_sigma) or float(mean) < 0:
			current_seq = -1
			continue

		effective_samples = param.ess_ratio * n_obs

		# true cdf:
		if n_obs * param.ess_ratio >= 2:
			cdf_val = calc_cdf(param.adjusted_mu, effective_samples, mean, stddev)
			p_value = calc_pvalue(cdf_val)



		if param.plots == True:
			#plot a chain of avg insert size every 500 bp for the first million bp
			if mean > 1:
				# starting to plot the isize chain window at this position 
				if not good_pos:
					good_pos = 0.5
			# take 2000 samples i a region of 1M bp.
			if good_pos and i %500 == 0 and good_pos < 2001:
				mean_isize.append(mean)
				good_pos += 1

			# calculate p_values and plots for:

			
			if n_obs * param.ess_ratio >= 2:
				#correcting both biases
				p_values_correct_both.append(p_value)
				cdf_values_correct_both.append(cdf_val)
			
				# correcting only bias 2
				cdf_val_2 = calc_cdf(param.mu, effective_samples, mean, stddev)
				p_value_2 = calc_pvalue(cdf_val_2)
				p_values_correct_bias2.append(p_value_2)
				cdf_values_correct_bias2.append(cdf_val_2)

			if n_obs >= 2:
				# no corrected bias
				cdf_val_1_2 = calc_cdf(param.mu, n_obs, mean, stddev)
				p_value_1_2 = calc_pvalue(cdf_val_1_2)
				p_values_naive.append(p_value_1_2)
				cdf_values_naive.append(cdf_val_1_2)

				# correcting only bias 1
				cdf_val_1 = calc_cdf(param.adjusted_mu, n_obs, mean, stddev)
				p_value_1 = calc_pvalue(cdf_val_1)
				p_values_correct_bias1.append(p_value_1)
				cdf_values_correct_bias1.append(cdf_val_1)




		if n_obs * param.ess_ratio < 2:
			continue

		if (scf != current_seq and pos >= param.max_window_size):
			current_seq = scf
			window = Window(param.pval, param.max_window_size, param.read_length)
			window.update(pos, p_value, mean)
		
		elif pos >= param.max_window_size:
			for significant_position in  window.update(pos, p_value, mean):
				avg_window_mean = window.avg_inner_mean + 2*window.read_length
				if avg_window_mean > param.adjusted_mu:
					sv_container.add_bp_to_cluster(current_seq, int(significant_position), window.avg_pval, n_obs, avg_window_mean, 'expansion', min(window.nr_in_window, window.max_window_size))
				else:
					sv_container.add_bp_to_cluster(current_seq, int(significant_position), window.avg_pval, n_obs, avg_window_mean, 'contraction', min(window.nr_in_window, window.max_window_size))
		if i%10000 == 0:
			print 'Processing coord',i, p_value, mean, stddev, n_obs, param.adjusted_mu

	if param.plots == True:
		plot_mean_chain(param, mean_isize)
		outfile = os.path.join(param.plotfolder,'p_values_corr_1_2')
		plot_stats(outfile, p_values_correct_both, title='p-value distribution', x='p-value',y='Frequency')
		outfile = os.path.join(param.plotfolder,'cdf_values_corr_1_2')
		plot_stats(outfile, cdf_values_correct_both, title='CDF-value distribution', x='CDF-value',y='Frequency' )

		outfile = os.path.join(param.plotfolder,'p_values_corr_1')
		plot_stats(outfile, p_values_correct_bias1, title='p-value distribution', x='p-value',y='Frequency')
		outfile = os.path.join(param.plotfolder,'cdf_values_corr_1')		
		plot_stats(outfile, cdf_values_correct_bias1, title='CDF-value distribution', x='CDF-value',y='Frequency' )

		outfile = os.path.join(param.plotfolder,'p_values_corr_2')
		plot_stats(outfile, p_values_correct_bias2, title='p-value distribution', x='p-value',y='Frequency')
		outfile = os.path.join(param.plotfolder,'cdf_values_corr_2')
		plot_stats(outfile, cdf_values_correct_bias2, title='CDF-value distribution', x='CDF-value',y='Frequency' )

		outfile = os.path.join(param.plotfolder,'p_values_naive')
		plot_stats(outfile, p_values_naive, title='p-value distribution', x='p-value',y='Frequency')
		outfile = os.path.join(param.plotfolder,'cdf_values_naive')
		plot_stats(outfile, cdf_values_naive, title='CDF-value distribution', x='CDF-value',y='Frequency' )


	sv_container.get_final_bp_info()
	print >> gff_file, str(sv_container)
