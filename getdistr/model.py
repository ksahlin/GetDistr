'''
Created on Sep 20, 2013

@author: ksahlin
'''


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
        w_fcn.append(max(a + b - o + 1, 0))
        w_fcn.append(max(min(a, b) - (r - s) + 1, 0))
    else:
        w_fcn.append(max(a - o + 1, 0))

    return(min(w_fcn))


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


class GeneralModel(object):
    def __init__(self, mu, sigma, r, s=None):
        self.mu = mu
        self.sigma = sigma
        self.r = r
        self.s = s if s != None else r / 2

    def expected_mean(self, z, a, b=None):
        E_x_given_z = 0

        norm_const = 0
        for t in range(z + 2 * (self.r - self.s), self.mu + 6 * self.sigma): #iterate over possible fragment sizes   ##range((self.mu - 5 * self.sigma) - y, self.mu + 6 * self.sigma - y): #
            norm_const += w(t - z , self.r, a, b, self.s) * stats.norm.pdf(t  , self.mu, self.sigma)

        for y in range(2 * (self.r - self.s), self.mu + 6 * self.sigma - z): # iterate over possible observation span
            weight = w(y, self.r, a, b, self.s)

            #norm_const = sum([weight * stats.norm.pdf(t + y, self.mu, self.sigma) for t in range(-(self.mu + 5 * self.sigma), self.mu + 5 * self.sigma) ])
            w_times_f = weight * stats.norm.pdf(z + y, self.mu, self.sigma)
            E_x_given_z += (y + z) * w_times_f / norm_const
            #print y, E_x_given_z, norm_const, w_times_f, (y + z)
        #sum([w(o, self.r, self.s) * stats.norm.pdf(o, mu, sigma) for o in  ])
        return(E_x_given_z)

    def expected_variance(self):
        raise NotImplementedError

    def infer_mean(self):
        raise NotImplementedError

    def infer_variance(self):
        raise NotImplementedError

class NormalModel(object):
    def __init__(self, o, r, a, b=None, s=None):
        self.o = o
        self.r = r
        self.a = a
        self.b = b
        self.s = s if s != None else r / 2
        return()

