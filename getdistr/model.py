'''
Created on Sep 20, 2013

@author: ksahlin
'''

import sys
import warnings

from scipy import stats
import pysam

from mpmath import *
mp.dps = 50 # decimal digits of precision

from scipy.stats import poisson, nbinom,uniform

from coverage import mean_span_coverage, Param
from mathstats.normaldist.truncatedskewed import param_est
#import math

def normpdf(x, mu, sigma):
    #Get much better approximations with Decimal (simply more decimals)
    #getcontext().prec = 100
    #u = Decimal(str(x - mu)) / Decimal(str(abs(sigma)))
    u = mpf((x - mu)) / mpf(sigma)
    #y = float(str((1 / Decimal(str((math.sqrt(2 * math.pi) * abs(sigma))))) * Decimal(str(-u * u / 2)).exp()))
    y = (1 / mpf((sqrt(2 * pi) * sigma)) * exp(mpf(-u * u / 2)))
    return y

def w(alignment, r, a, s_inner,s_outer, b=None, infer_lib_mean=False, breakdancer=False):
    '''
        Calculate the weight for a given observation.

        Parameters:

        alignment       -- Position and number of softclipped bases
                        for an alignment. A tuple with three integers 
                        (observation, observed inner softclipped bases,
                            observed outer softclipped bases)
        r               -- The read length.
        a               --reference length 


        Returns the weight that this observation has.

    '''

    o, s_inner_obs, s_outer_obs = alignment[0],alignment[1], alignment[2]

    ## 
    # Check so that number of allowed softclipped bases 
    # (specified by user) is greater than the obsered 
    # softclipped bases (entered by the user)
    if 2*s_inner < s_inner_obs or 2*s_outer < s_outer_obs:
        warnings.warn("Warning: Observed more softclipped bases (inner softclipped: {0},\
                     outer softclipped: {1} bp) in the alignemnt tuple than what is \
                      specified to be allowed ({1} inner and {2} outer total bp) with \
                      parameters 's_inner' and 's_outer'. Please \
                      reconsider this parameter".format(s_inner_obs,s_outer_obs,s_inner,s_outer))
        sys.exit()


    w_fcn = []

    s_param_inner = 2*s_inner - s_inner_obs
    s_param_outer = 2*s_outer - s_outer_obs

    ##
    # Weight function in the case where we want to estimate the original library mean
    # (full original distribution) from alignments of paired reads onto a contig
    if infer_lib_mean:
        return(max(a - (o - s_param_outer) + 1, 0))


    ##
    # Weight function in the case where we want to calculate expected or inferred mean of
    # observations coming from our library with truncation and skewness 
    # (e.g. over insertions or scaffold gaps).


    # inner boundary
    w_fcn.append(max(o - 2 *(r-s_inner) + 1, 0))

    if breakdancer:
        if o > a:
            return 0
        else:
            return min(w_fcn)

    if b:
        # outer boundary
        w_fcn.append(max(a + b - (o -  s_param_outer) + 1, 0))
        # shorter ref
        # NOTE: If the smallest reference sequence here is smaller than
            # the read length, this formula is incorrect and we need keep track
        # of separate softclippes between the reads. Ref seq < read length is
        # however not recommented to use this model on as estimation have hight 
        # unceratianty so we therefore do not treat this case. 
        w_fcn.append(max(min(a, b) - (r - s_inner) + 1, 0))
    else:
        w_fcn.append(max(a - o + s_param_outer + 1, 0))

    return(min(w_fcn))


def get_possible_positions_dict(bam_file, softclipped, low, high):
    """
    This function could be highly optimized.
    """
    ref_lengths = sorted(map(lambda x: int(x),bam_file.lengths), reverse=True)
    print ref_lengths
    position_dict = {}
    print low, high
    if low < 0:
        low = 0 
    i = low
    while i <= high:
        position_dict[i] = 0
        j = 0
        while j < len(ref_lengths) and i <= ref_lengths[j]:
            position_dict[i] += w(i, 0, ref_lengths[j], s=softclipped, infer_lib_mean=True)
            j += 1

        i += 1

    #print position_dict
    return position_dict

def estimate_library_parameters(bam_path, list_of_obs, low, high, soft=None):
    sample_obs_sum = 0
    sample_obs_sum_sq = 0
    number_of_obs_sum = 0
    number_of_obs_sum_sq = 0
    with pysam.Samfile(bam_path, 'rb') as bam_file:
        w = get_possible_positions_dict(bam_file, soft, low, high)

    for o in list_of_obs:
        weight = float(w[o])
        sample_obs_sum += o / weight
        number_of_obs_sum += 1 / weight
        sample_obs_sum_sq  += o**2/weight #(o / weight)**2
        number_of_obs_sum_sq += 1 / weight**2
    # formulas for mean and variance
    mu = sample_obs_sum / number_of_obs_sum
    sigma = sqrt(sample_obs_sum_sq /number_of_obs_sum - mu**2 )      #/ number_of_obs_sum_sq - mu**2)
    return(mu,sigma)

# def estimate_library_stddev(list_of_obs, mu, weighted_sum, r, a, soft=None):
#     sample_obs_sum = 0
#     number_of_obs_sum = 0
#     for o in list_of_obs:
#         weight = float(w(o, r, a, s=soft, infer_lib_mean=True))
#         sample_obs_sum += o / weight
#         number_of_obs_sum += 1 / weight
#     print sample_obs_sum / number_of_obs_sum
#     return(sample_obs_sum / number_of_obs_sum)   

#     raise NotImplementedError




class NormalModel(object):
    ''' 
        Derives the most probable observed library distribution
        given a full library distribution that is normally distributed.

    '''

    def __init__(self, mu, sigma, r, s_inner=None, s_outer=None, breakdancer=False):

        self.s_inner = s_inner if s_inner != None else r / 2
        self.s_outer = s_outer if s_outer != None else r / 2

        # check that all given parameters are integers
        if not all( map( lambda param: type(param) == int, [mu, sigma, r, self.s_inner, self.s_outer])):
            warnings.warn("Warning: All parameters specified to NormalModel needs to be integers.\
                you specified: mu:{0}, sigma:{1}, read length: {2} converting to closest integer value".format(mu,sigma,r))           
        
        self.mu = int(mu)
        self.sigma = int(sigma)
        self.r = int(r)

        self.breakdancer = breakdancer


    def expected_mean(self, z, a, b=None):
        '''
        Returns expected mean fragment size length given a gap of length z and reference lengths a,b.
        It is assumed that you have my,sigma of the underlying library and the read length r.
        '''
        E_x_given_z = 0
        norm_const = 0
        for t in range(z + 2 * (self.r - self.s_inner), self.mu + 6 * self.sigma): #iterate over possible fragment sizes   ##range((self.mu - 5 * self.sigma) - y, self.mu + 6 * self.sigma - y): #
            norm_const += w((t - z,0,0) , self.r, a, self.s_inner, self.s_outer, b=b) * normpdf(t + 0.5, self.mu, self.sigma)  # +0.5 because we approximate a continuous distribution (avg function value of pdf given points i and i+1, just like integration)

        for y in range(2 * (self.r - self.s_inner), self.mu + 6 * self.sigma - z): # iterate over possible observation span
            weight = w((y,0,0) , self.r, a, self.s_inner, self.s_outer, b=b)
            w_times_f = weight * normpdf(z + y + 0.5, self.mu, self.sigma) # +0.5 because we approximate a continuous distribution (avg function value of pdf given points i and i+1, just like integration)
            E_x_given_z += (y + z) * w_times_f / norm_const
        return(E_x_given_z)

    def expected_standard_deviation(self, z, a, b=None):
        '''
        Returns expected standard deviation of fragment size length given a gap of length z and reference lengths a,b.
        It is assumed that you have my,sigma of the underlying library and the read length r.
        '''
        E_x_given_z = self.expected_mean(z, a, b)
        E_x_square_given_z = 0
        norm_const = 0
        for t in range(z + 2 * (self.r - self.s_inner), self.mu + 6 * self.sigma): #iterate over possible fragment sizes   ##range((self.mu - 5 * self.sigma) - y, self.mu + 6 * self.sigma - y): #
            norm_const += w((t - z,0,0) , self.r, a, self.s_inner,self.s_outer,  b) * stats.norm.pdf(t + 0.5, self.mu, self.sigma)  # +0.5 because we approximate a continuous distribution (avg function value of pdf given points i and i+1, just like integration)

        for y in range(2 * (self.r - self.s_inner), self.mu + 6 * self.sigma - z): # iterate over possible observation span
            weight = w((y,0,0) , self.r, a,  self.s_inner, self.s_outer, b = b)
            w_times_f = weight * stats.norm.pdf(z + y + 0.5, self.mu, self.sigma) # +0.5 because we approximate a continuous distribution (avg function value of pdf given points i and i+1, just like integration)
            E_x_square_given_z += (y + z) ** 2 * w_times_f / norm_const

        Var_x_given_z = (E_x_square_given_z - E_x_given_z ** 2) ** 0.5
        return(Var_x_given_z)



    def infer_mean_slow(self, list_of_obs, a, precision, b=None, coverage = False, n = False, coverage_model = False):
        '''
            Instance method of a NormalModel object. Infers the mean fragment size of a given set of 
            paired read observations O(n^2) complexity. This function is only included (or still remaining)
            for historical and verification purposes. It performes an exhaustive search for optimum.
            Use infer_mean_fast for O(n log n) in complexity.
            
            Keyword arguments:
            @ argument list_of_obs: A list of...
            @ argument a:           Reference sequence length
            @ argument precision:   Number of base pairs between every point estimate of the ML distribution. 
            @ argument coverage:    Mean coverage.
            @ argument coverage_model:  The assumed coverage distribution around the mean. Valid strings are 'Uniform',
                                        'NegBin' or 'Poisson'.

            Returns:
            Maximum likelihood value for fragment length X.
        '''
        warnings.warn("Warning: This function is deprecated. Use 'infer_mean_fast' or 'run_GapEst' for fast calculation.")

        likelihood_curve = self.get_likelihood_function(list_of_obs, a, precision, b, coverage, n, coverage_model)
        ml_gap = max(likelihood_curve, key=lambda x: x[1])
        #print likelihood_curve, ml_gap
        avg_obs = sum(list_of_obs) / len(list_of_obs)  # avg_obs is an integer (rounded to an even bp)
        #print avg_obs + ml_gap[0]
        return(avg_obs + ml_gap[0])

    def infer_mean_fast(self, list_of_obs, a, precision, b=None, coverage = False, n = False, coverage_model = False):
        '''
            Instance method of a NormalModel object. Infers the mean fragment size of a given set of 
            paired read observations using .
            
            Keyword arguments:
            @ argument list_of_obs: A list of...
            @ argument a:           Reference sequence length
            @ argument precision:   Number of base pairs between every point estimate of the ML distribution. 
            @ argument coverage:    Mean coverage.
            @ argument coverage_model:  The assumed coverage distribution around the mean. Valid strings are 'Uniform',
                                        'NegBin' or 'Poisson'.
        '''

        #do binary search among limited range of values
        
        z_upper= self.mu + 3 * self.sigma - (2 * (self.r - self.s_inner))
        z_lower=-3 * self.sigma
        z_u = int((z_upper+z_lower)/2.0 + precision)
        z_l = int((z_upper+z_lower)/2.0)
        while z_upper-z_lower>precision:
            
            likelihood_of_z_u = self.get_likelihood_value(z_u, list_of_obs, a, b, coverage,n , coverage_model)
            likelihood_of_z_l = self.get_likelihood_value(z_l, list_of_obs, a, b, coverage,n , coverage_model)
            if likelihood_of_z_u>likelihood_of_z_l:
                z_lower = z_u
                z_u = int((z_upper+z_u)/2.0 + precision)
                z_l = z_u -precision

            else:
                z_upper = z_l
                z_u = int((z_lower+z_l)/2.0 + precision)
                z_l = z_u - precision

        return(max(z_u,z_l))


    def run_GapEst(self, list_of_obs, a, b=None):
        mean_obs = float(sum(list_of_obs))/len(list_of_obs)
        if b == None:
            # GapEst only works with two reference sequences
            # Just split a into two sequences. Hopefully a is long
            # enough compared to the mean and std dev of the library
            # so that this does not affect the result too much.
            if a < 2*self.mu +6*self.sigma:
                warnings.warn("Warning: Unsafe to run GapEst on one reference sequence\
                (it is implemented for two reference sequences) that is shorter than \
                    2*mean + 6sigma of library insert size. Use infer_mean_fast instead.")
                sys.exit()  

            # reference is long enough for safe estimations            
            b = a/2
            a = a/2
            mean_obs = float(sum(list_of_obs))/len(list_of_obs)
            # gapest does not handle softclipping so we send in self.r - self.s_inner as the effective read length
            # note that this in not correct for short reference sequences, then we wold want to send in
            # self.s_inner. This is not treated here however and in general s_inner = s_outer
            return param_est.GapEstimator(self.mu, self.sigma, self.r - self.s_inner, mean_obs, a, b)
        else:
            return param_est.GapEstimator(self.mu, self.sigma, self.r - self.s_inner, mean_obs, a, b)

    def infer_variance(self):
        raise NotImplementedError

    def get_likelihood_value(self, z, list_of_alignmet_tuples, a, b=None, coverage = False,n = False, coverage_model = False):
        '''
        This function gives back a single likelihood value from the likelihood function
            parameters
            __________
            @param list_of_obs_tuples:         A list of observations
            @param a:                   Reference sequence length
            @param precision: 
            @ argument coverage:        Mean coverage.
            @ argument coverage_model:  The assumed coverage distribution around the mean. Valid strings are 'Uniform',
                                        'NegBin' or 'Poisson'.
        '''
        list_of_obs = map(lambda x : x[0], list_of_alignmet_tuples)

        if coverage:
            n = len(list_of_obs)
            if not coverage_model:
                warnings.warn("Warning: Coverage parameter is set to True but no 'coverage_model' parameter set. \
                defaulting to uniform coverage. ")
        ##
        # calculate the normalization constant for a given gap length z
        norm_const = 0
        for t in range(z + 2 * (self.r - self.s_inner), self.mu + 7 * self.sigma): #iterate over possible fragment sizes   ##range((self.mu - 5 * self.sigma) - y, self.mu + 6 * self.sigma - y): #
            #norm_const += w(t - z , self.r, a, self.s_inner, self.s_outer, b=b) * stats.norm.pdf(t + 0.5, self.mu, self.sigma)  # +0.5 because we approximate a continuous distribution (avg function value of pdf given points i and i+1, just like integration)
            norm_const += w((t - z, 0, 0) , self.r, a, self.s_inner, self.s_outer, b=b, breakdancer=self.breakdancer) * normpdf(t + 0.5, self.mu, self.sigma)
        # eventual softclipped outer reads, they will have the same observation "a" but be of different
        #lengths due to outer softclipped bases
        if a < self.mu + 7 * self.sigma:
            for i in range(1, min(2*self.s_outer, (self.mu + 7 * self.sigma - a ) )):
                norm_const += w((a, 0, i) , self.r, a, self.s_inner,self.s_outer, b=b) * normpdf(z+a+i + 0.5, self.mu, self.sigma)


        ##
        # calculate the nominator (relative frequency given a gap)
        # in log() format
        log_p_x_given_z = 0
        for o,s_inner_obs,s_outer_obs in list_of_alignmet_tuples:
            weight = w( (o,s_inner_obs,s_outer_obs) , self.r, a, self.s_inner, self.s_outer, b=b, breakdancer=self.breakdancer)
            lib_dist = normpdf(o + s_outer_obs + z + 0.5, self.mu, self.sigma)
            #print z, weight, lib_dist, norm_const
            log_p_x_given_z += log(weight) + log(lib_dist) - log(norm_const)

        if coverage:
            p_n_given_z = self.coverage_probability(n, a, self.mu, self.sigma, z, coverage, self.r, self.s_inner,self.s_outer, b, coverage_model)
            log_p_n_given_z = log(p_n_given_z)/n
            return log_p_x_given_z + log_p_n_given_z
        else:
            return log_p_x_given_z


    # Move out: for z in range(-3 * self.sigma, self.mu + 3 * self.sigma - (2 * (self.r - self.s)), precision):
    # and put it inside infer_mean function??
    # also add a function (separate out): "get_likelihood_function" so that plotting scripts still works

    def get_likelihood_function(self, list_of_obs, a, precision, b=None, coverage = False,n = False, coverage_model = False):
        '''
        This function gives back the likelihood values for Z (gap/unknown sequence length)
            parameters
            __________
            @param list_of_obs:         A list of observations
            @param a:                   Reference sequence length
            @param precision: 
            @ argument coverage:        Mean coverage.
            @ argument coverage_model:  The assumed coverage distribution around the mean. Valid strings are 'Uniform',
                                        'NegBin' or 'Poisson'.


        '''
        if coverage:
            n = len(list_of_obs)
            if not coverage_model:
                warnings.warn("Warning: Coverage parameter is set to True but no 'coverage_model' parameter set. \
                defaulting to uniform coverage. ")

        likelihood_curve = []

        ## For plotting!
        #o_temp =  list_of_obs[0] # <-- Only for plotting
        # z in range( 3 * self.sigma - o_temp, - o_temp + self.mu + int(5.5 * self.sigma) - (2 * (self.r - self.s)), precision): 
        # z in range( 3 * self.sigma - o_temp, - o_temp + self.mu + int(5 * self.sigma) - (2 * (self.r - self.s)), precision): 
        # <-- Use the above range for plotting with same intervals, this function gives back likelihood alues of 
        # X instead of Z which is what we want in the MLfcns plots.

        ##
        # This loop iterates over all possible gaps z, we want to see the ML estimation of
        # The interesting range is in general not above  mean + 3*stddev

        for z in range(-3 * self.sigma, self.mu + 3 * self.sigma - (2 * (self.r - self.s_inner)), precision): 
            likelihood_curve.append((z,self.get_likelihood_value(z, list_of_obs, a, b, coverage,n , coverage_model)))

        return(likelihood_curve)


    def coverage_probability(self,nr_obs, a, mean_lib, stddev_lib,z, coverage_mean, read_len, s_inner,s_outer, b=None, coverage_model = False):
        ''' Distribution P(o|c,z) for prior probability over coverage.
            This probability distribution is implemented as an poisson 
            distribution.

            Attributes:

            c       -- coverage
            mean    -- mean value of poisson distribution.

            Returns probability P(c)

        '''
        if not b: 
            # only one reference sequence.
            # We split the reference sequence into two equal 
            # length sequences to fit the model. 
            a = a/2
            b = a/2

        param = Param(mean_lib, stddev_lib, coverage_mean, read_len, s_inner,s_outer)
        lambda_ = mean_span_coverage(a, b, z, param)

        if coverage_model == 'Poisson':
            return poisson.pmf(nr_obs, lambda_, loc=0)

        elif coverage_model == 'NegBin':
            p = 0.01
            n = (p*lambda_)/(1-p)
            return nbinom.pmf(nr_obs, n, p, loc=0) 
        else:
            # This is equivalent to uniform coverage
            return 1 #uniform.pdf(nr_obs, loc=lambda_- 0.3*lambda_, scale=lambda_ + 0.3*lambda_ )


         




class GeneralModel(object): 
    def __init__(self, histogram, r, s=None):
        self.histogram = histogram
        self.r = r
        self.s = s if s != None else r / 2

    def expected_mean(self, z, a, b=None):

        return()

