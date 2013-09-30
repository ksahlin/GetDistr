'''
Created on Sep 20, 2013

@author: ksahlin
'''

from math import log

from scipy import stats

from decimal import Decimal, getcontext
import math

def normpdf(x, mu, sigma):
    #Get much better approximations with Decimal (simply more decimals)
    getcontext().prec = 100
    u = Decimal(str(x - mu)) / Decimal(str(abs(sigma)))
    y = float(str((1 / Decimal(str((math.sqrt(2 * math.pi) * abs(sigma))))) * Decimal(str(-u * u / 2)).exp()))
    return y

def w(o, r, a, b=None, s=None):
    s = s if s != None else r / 2
    w_fcn = []
    w_fcn.append(max(o - 2 * (r - s) + 1, 0))
    if b:
        w_fcn.append(max(a + b - (o - 2 * s) + 1, 0))
        w_fcn.append(max(min(a, b) - (r - s) + 1, 0))
    else:
        w_fcn.append(max(a - (o - 2 * s) + 1, 0))

    return(min(w_fcn))

def estimate_library_mean(list_of_obs, r, a, soft=None):
    sample_obs_sum = 0
    number_of_obs_sum = 0
    for o in list_of_obs:
        weight = float(w(o, r, a, s=soft))
        sample_obs_sum += o / weight
        number_of_obs_sum += 1 / weight

    print sample_obs_sum / number_of_obs_sum
    return(sample_obs_sum / number_of_obs_sum)

def estimate_library_stddev(list_of_obs, r, a, soft=None):
    raise NotImplementedError



def dw(o, r, a, b=None, s=None):
    raise NotImplementedError

def ddw(o, r, a, b=None, s=None):
    raise NotImplementedError

def f_general(o, r, a, b=None, s=None):

    raise NotImplementedError

def f_normal(x, mu, sigma):
    stats.norm.pdf(x, mu, sigma)
    raise NotImplementedError

def df_normal(o, r, a, b=None, s=None):
    raise NotImplementedError

class Weight(object):
    def __init__(self, o, r, a, b=None, s=None):
        self.o = o
        self.r = r
        self.a = a
        self.b = b
        self.s = s if s != None else r / 2
        return()

    def w(self):
        pass


    def dw(self):
        raise NotImplementedError


class NormalModel(object):
    def __init__(self, mu, sigma, r, s=None):
        self.mu = mu
        self.sigma = sigma
        self.r = r
        self.s = s if s != None else r / 2

    def expected_mean(self, z, a, b=None):
        E_x_given_z = 0
        norm_const = 0
        for t in range(z + 2 * (self.r - self.s), self.mu + 6 * self.sigma): #iterate over possible fragment sizes   ##range((self.mu - 5 * self.sigma) - y, self.mu + 6 * self.sigma - y): #
            norm_const += w(t - z , self.r, a, b, self.s) * stats.norm.pdf(t + 0.5, self.mu, self.sigma)  # +0.5 because we approximate a continuous distribution (avg function value of pdf given points i and i+1, just like integration)

        for y in range(2 * (self.r - self.s), self.mu + 6 * self.sigma - z): # iterate over possible observation span
            weight = w(y , self.r, a, b, self.s)
            w_times_f = weight * stats.norm.pdf(z + y + 0.5, self.mu, self.sigma) # +0.5 because we approximate a continuous distribution (avg function value of pdf given points i and i+1, just like integration)
            E_x_given_z += (y + z) * w_times_f / norm_const
        return(E_x_given_z)

    def expected_standard_deviation(self, z, a, b=None):
        E_x_given_z = self.expected_mean(z, a, b)
        E_x_square_given_z = 0
        norm_const = 0
        for t in range(z + 2 * (self.r - self.s), self.mu + 6 * self.sigma): #iterate over possible fragment sizes   ##range((self.mu - 5 * self.sigma) - y, self.mu + 6 * self.sigma - y): #
            norm_const += w(t - z , self.r, a, b, self.s) * stats.norm.pdf(t + 0.5, self.mu, self.sigma)  # +0.5 because we approximate a continuous distribution (avg function value of pdf given points i and i+1, just like integration)

        for y in range(2 * (self.r - self.s), self.mu + 6 * self.sigma - z): # iterate over possible observation span
            weight = w(y , self.r, a, b, self.s)
            w_times_f = weight * stats.norm.pdf(z + y + 0.5, self.mu, self.sigma) # +0.5 because we approximate a continuous distribution (avg function value of pdf given points i and i+1, just like integration)
            E_x_square_given_z += (y + z) ** 2 * w_times_f / norm_const

        Var_x_given_z = (E_x_square_given_z - E_x_given_z ** 2) ** 0.5
        return(Var_x_given_z)



    def infer_mean(self, list_of_obs, a, precision, b=None):
        '''
            Instance method of a NormalModel object. Infers the mean fragment size of a given set of 
            paired read observations.
            
            Keyword arguments:
            @ argument list_of_obs A list of...
            @ argument a Reference sequence length
            @ argument precision ... 
        '''

        likelihood_curve = []
        for z in range(-3 * self.sigma, self.mu + 3 * self.sigma - (2 * (self.r - self.s)), precision):
            ##
            # claculate the normalization constant for a given gap length z
            norm_const = 0
            for t in range(z + 2 * (self.r - self.s), self.mu + 6 * self.sigma): #iterate over possible fragment sizes   ##range((self.mu - 5 * self.sigma) - y, self.mu + 6 * self.sigma - y): #
                norm_const += w(t - z , self.r, a, b, self.s) * stats.norm.pdf(t + 0.5, self.mu, self.sigma)  # +0.5 because we approximate a continuous distribution (avg function value of pdf given points i and i+1, just like integration)


            ##
            # calculate the nominator (relative frequency given a gap)
            # in log() format
            log_p_x_given_z = 0
            for o in list_of_obs:
                weight = w(o , self.r, a, b, self.s)
                lib_dist = stats.norm.pdf(o + z + 0.5, self.mu, self.sigma)
                #print z, weight, lib_dist, norm_const
                log_p_x_given_z += log(weight) + log(lib_dist) - log(norm_const)

            likelihood_curve.append((z, log_p_x_given_z))

        ml_gap = max(likelihood_curve, key=lambda x: x[1])
        #print likelihood_curve, ml_gap
        avg_obs = sum(list_of_obs) / len(list_of_obs)  # avg_obs is an integer (rounded to an even bp)
        print avg_obs + ml_gap[0]
        return(avg_obs + ml_gap[0])

    def infer_variance(self):
        raise NotImplementedError



class GeneralModel(object):
    def __init__(self, histogram, r, s=None):
        self.histogram = histogram
        self.r = r
        self.s = s if s != None else r / 2

    def expected_mean(self, z, a, b=None):

        return()

