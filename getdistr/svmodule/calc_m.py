'''
Created on Sep 18, 2013

@author: ksahlin
'''

import argparse
from getdistr import coverage
import math

from scipy.stats import norm


def main(args):
	# calculate expected number of observations
	param = coverage.Param( args.mu, args.sigma, args.cov, args.r, args.s, args.s)
	n = coverage.mean_span_coverage(args.G/2, args.G/2, 0, param)
	print n
	corrected_significance = args.alpha/(2*args.G) # Bonferroni
	quantile_a = norm.ppf( 1 - corrected_significance )
	m = quantile_a*math.sqrt(args.sigma**2/n)
	print m



if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Calculate an approximate value of m (smallest variant) possible to detect given a library and significance level.")
	parser.add_argument('alpha', type=float, help='significance level (uncorrected for multiple tests), e.g. 0.05 ')
	parser.add_argument('G', type=float, help='Length of the genome (or number of significance tests). ')
	parser.add_argument('mu', type=float, help='Mean of library')
	parser.add_argument('sigma', type=float, help='Stddev of library')
	parser.add_argument('cov', type=float, help='Mean coverage of library')
	parser.add_argument('r', type=float, help='read length  of library')
	parser.add_argument('s', type=float, help='Maximum allowed softclipped bases')

	args = parser.parse_args()
	main(args)



        
