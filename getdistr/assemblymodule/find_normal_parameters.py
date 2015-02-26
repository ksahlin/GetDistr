import argparse
from scipy.stats.distributions import norm
#import math
import numpy as np
from scipy.optimize import fmin

try:
	import matplotlib
	matplotlib.use('agg')
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab
	import seaborn as sns
	sns.set_palette("husl", desat=.6)
except ImportError:
	pass

def estimate_params_for_normal(x, low_bound , mu_initial, sigma_initial):
	"""
		Takes a vector x of truncated data with a known lower
		truncation bound and estimates the parameters of the 
		fit of an untruncated normal distribution.
		code from Chris Fonnesbeck's Python data analysis tutorial on Sense
		https://sense.io/prometheus2305/data-analysis-in-python/files/Statistical%20Data%20Modeling.py
	"""


	# normalize vector
	mu_initial = float(mu_initial)
	sigma_initial = float(sigma_initial)
	#x = np.random.normal(size=10000,loc=2000,scale= 2000)

	x = map(lambda y: (y-mu_initial )/sigma_initial ,x)
	a =  (low_bound - mu_initial)/sigma_initial # normalize lower bound
	

	#_ = plt.hist(x, bins=100)
	#plt.show()
	#plt.close()

	# We can construct a log likelihood for this function using the conditional
	# form	
	trunc_norm = lambda theta, a, x: -(np.log(norm.pdf(x, theta[0], theta[1])) - 
	                                      np.log(1 - norm.cdf(a, theta[0], theta[1]))).sum()

	# For this example, we will use another optimization algorithm, the
	# **Nelder-Mead simplex algorithm**. It has a couple of advantages: 
	# 
	# - it does not require derivatives
	# - it can optimize (minimize) a vector of parameters
	# 
	# SciPy implements this algorithm in its `fmin` function:

	# we have normalized data, given that the loer truncation point a
	# is pretty far out in the tail - the standard normal parameters are
	# a first good guess, i.e. 0,1
	initial_guess = np.array([0,1]) 
	sol = fmin(trunc_norm, initial_guess , args=(a, x))
	print sol
	mean_normalized,stddev_normalized = sol[0],sol[1]
	mean_est =( 1 + mean_normalized ) * mu_initial
	stddev_est = stddev_normalized * sigma_initial
	print mean_est,stddev_est
	return mean_est,stddev_est



def main(x,outfile,bam_file = False, bp_file=False):
	n = float(len(x))
	mu = sum(x)/n
	sigma = (sum(list(map((lambda t: t ** 2 - 2 * t * mu + mu ** 2), x))) / (n - 1)) ** 0.5

	mu_est, sigma_est = estimate_params_for_normal(x, 0, mu, sigma )

	#plt fit
	_ = plt.hist(x, bins=200, normed=True, alpha=0.6, color='g')
	xmin, xmax = plt.xlim()
	#mu, sigma = norm.fit(x)
	#print mu,sigma
	y = np.linspace(xmin, xmax, 100)
	plt.plot(y,mlab.normpdf(y,mu_est, sigma_est),'r')
	title = "Fit results: mu = %.2f,  std = %.2f" % (mu_est, sigma_est)
	plt.title(title)
	plt.savefig(outfile)
	plt.close()
	plt.clf()


if __name__ == '__main__':
	from src import CreateGraph_updated, GapCalculator
	parser = argparse.ArgumentParser()
	parser.add_argument('--bam_file', dest='bamfile',type=str, default=False, help='bamfile. ')
	parser.add_argument('--bp_file', dest='bp_file', type=str, default = False, help='Base pair stats file generated by GetDistr-assembly module. ')
	parser.add_argument('--out', dest='outfile', type=str, help='outfile. ')
	args = parser.parse_args()
	is_bam = args.bamfile
	is_bp = args.bp_file

	main(args.outfile, bam_file=is_bam, bp_file=is_bp)

