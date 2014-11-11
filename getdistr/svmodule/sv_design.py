#from math import exp
import argparse
from scipy.stats import poisson,norm, binom
import math
from decimal import Decimal, getcontext

from coverage import mean_span_coverage, Param


# poisson distribution
def P_breakpoints_in_interval(I, q, n):
	"""
		param:
		q		- breakpoint ratio
		I 		- Interval of lenght I
		n 		- number of break points 

		Calculates:
		k		- Number of expected breakpoints within an interval
		P( n bp in I) 	-	Probability of n breakpoints in I

		returns: 
		P( n bp in I)
	"""
	k = q * I

	# for i in range(0,n):


	# math.exp(-k)
	#print poisson.cdf(n, k)
	return poisson.cdf(n, k)



def P_number_true_obs_slow(args):
	param = Param(args.mean, args.stddev, args.cov, args.readlen, args.soft, args.soft) 
	E_true_links = 0
	# the diagonal of the square (i_1 = i_2) i.e. do not multiply by two
	p_at_least_one_bp_at_given_position = 1- P_breakpoints_in_interval(1, args.bp_ratio, 0)
	for i in range((args.readlen - args.soft),args.mean + 4*args.stddev +1):
		print i
		p = P_breakpoints_in_interval(i, args.bp_ratio, 0)
		if i != args.mean + 4*args.stddev:
			E_true_links += p**2 * p_at_least_one_bp_at_given_position**2 * mean_span_coverage(i, i, args.insertion_size, param)
		else:
			E_true_links += p**2 * mean_span_coverage(i, i, args.insertion_size, param)


	for i_1 in range((args.readlen - args.soft),args.mean + 4*args.stddev):
		for i_2 in range(i_1+1, args.mean + 4*args.stddev + 1):
			#print i_1,i_2
			p_1 = P_breakpoints_in_interval(i_1, args.bp_ratio, 0)
			p_2 = P_breakpoints_in_interval(i_2, args.bp_ratio, 0)
			# symmetrical probabilities
			if i_2 != args.mean + 4*args.stddev:

				E_true_links += 2 * p_1 * p_2 * p_at_least_one_bp_at_given_position**2 * mean_span_coverage(i_1, i_2, args.insertion_size, param)

			else:
				E_true_links += 2 * p_1 * p_2 * p_at_least_one_bp_at_given_position * mean_span_coverage(i_1, i_2, args.insertion_size, param)

	return E_true_links


def v(x,d,r,s):
	"""
		The volume of the geometrical shape or "3d"-weight. It is the 2d trapedzoid weight-function
		when adding a third dimension that is the width of the 2d trapetzoid. This third dimension
		is derived from looking att all possibe zises of I_1 and I_2. Number of trapetziods are 
		(x -d-r+s) for each x > d+r-s. The sum of weights of each such trapetziod at position x-d-r+s is
		sum_of_ints(1,x-d-r+s) = sum(1,n) = n(n+1)/2
	"""
	return x**2/2.0 -x*(d +r-s-1/2.0)/2.0 + 1/2.0*(d+r-s)*(d+r-s-1)


def vol(x,M):
	"""
	 M - maximum interval size (from 1 to say, mu + 4sigma)
	 x - insert size observation length
	"""
	#return max(0, M*M - 2*(M-(x-1)) * (x-1) + (M-(x-1))*(x-1)*x )
	return M*M - 2*(M-(x-1)) * (x-1) + (M-(x-1))*(x-1)*x

def vol2(x,M,c):
	"""
	 M - maximum interval size (from 1 to say, mu + 4sigma)
	 x - insert size observation length
	 c - 2(r-s)
	"""
	def sum_of_squares(n):
		return n*(n+1)*(2*n+1)/6.0

	def sum_of_natural_numbers(n):
		return n*(n+1)/2.0


	return sum_of_squares(M) - sum_of_squares(M-x+c-1) - sum_of_squares(x-c)/2.0 -sum_of_natural_numbers(x-c)/2.0 

def P_number_true_obs_fast(args):
	E_true_links = 0
	p_at_least_one_bp_at_given_position = 1- P_breakpoints_in_interval(1, args.bp_ratio, 0)
	k = 2 * args.readlen / float(args.cov)
	
	for i in range((args.insertion_size + args.readlen - args.soft)+1, args.mean + 4*args.stddev):
		# When internal breakpoint occur within mean + 4*sigma
		E_true_links += 2*(1/k) * (v(i,args.insertion_size, args.readlen, args.soft) - 1 ) * norm.pdf(i, args.mean,args.stddev)*poisson.pmf(0,args.bp_ratio*i)*p_at_least_one_bp_at_given_position**2
		# when no breakpoint occurs on one side
		E_true_links += 2*(1/k) * 1                                                        * norm.pdf(i, args.mean,args.stddev)*poisson.pmf(0,args.bp_ratio*i)*p_at_least_one_bp_at_given_position
		# when no breakpoint occurs on both sides
		E_true_links += (1/k) * 1 															* norm.pdf(i, args.mean,args.stddev)*poisson.pmf(0,args.bp_ratio*i)

		#print v(i,args.insertion_size, args.readlen, args.soft)
	# when no breakpoint occurs on one side
	i = args.mean + 4*args.stddev
	E_true_links += 2*(1/k)*v(i,args.insertion_size, args.readlen, args.soft)*norm.pdf(i, args.mean,args.stddev)*poisson.pmf(0,args.bp_ratio*i)*p_at_least_one_bp_at_given_position

	# when no breakpoint occurs on both sides
	i = args.mean + 4*args.stddev
	print v(i,args.insertion_size, args.readlen, args.soft)
	print 1/k
	E_true_links += (1/k)*v(i,args.insertion_size, args.readlen, args.soft)*norm.pdf(i, args.mean,args.stddev)*poisson.pmf(0,args.bp_ratio*i)

 
	return E_true_links


if __name__ == "__main__":
	arg_parser = argparse.ArgumentParser(description="Gives different optimal values for designing libraries and methods regarding structral variation detection. ")
	arg_parser.add_argument("bp_ratio", type=float, help="Breakpoint rate as a float number in (0,1) (number of\
		breakpoints divided by genome length.) ")
	arg_parser.add_argument("mean", type=int, help="Mean insert size.")
	arg_parser.add_argument("stddev", type=int, help="Standard deviation of insert size.")
	arg_parser.add_argument("cov", type=float, help="Mean coverage.")
	arg_parser.add_argument("readlen", type=int, help="Read length.")
	arg_parser.add_argument("soft", type=int, help="Number of softclipped bases allowed.")
	# arg_parser.add_argument("len1", type=int, help="Contig1 length.")
	# arg_parser.add_argument("len2", type=int, help="Contig2 length.")
	arg_parser.add_argument("insertion_size", type=int, help="Insertion size")



	args = arg_parser.parse_args()

	##
	# inner and outer softclipps are set to the same value as this is most likely for all applicatons
	#param = Param(args.mean, args.stddev, args.cov, args.readlen, args.soft, args.soft)

	print P_number_true_obs_slow(args)
	print P_number_true_obs_fast(args)

	print 'Vol:', vol(6,4)
	print 'Vol2:', vol2(500,4000,100)


