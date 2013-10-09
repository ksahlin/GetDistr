'''
Created on Sep 25, 2013

@author: ksahlin
'''

import sys

import argparse

import model
from mathstats.normaldist.truncatedskewed import param_est as GC

# frag_obj = model.NormalModel(mean, stddev, read length, softclipped)
# frag_obj.infer_mean(list of obs, reference length, bp-precision)

# Infer clever inserts
for line in open(sys.argv[1]):
    o = list(line.split())
    o = [int(i) for i in o]
    print o
    mean_obs = sum(o) / len(o)
    print 'true: ', mean_obs + 300
    #avg_gap = GC.GapEstimator(500, 50, 100 , mean_obs, 10000, 10000)
    #print 'gapest:', avg_gap, avg_gap + mean_obs, 'true: ', mean_obs + 300
    frag_obj = model.NormalModel(500, 50, 100, 0)
    frag_size_est = frag_obj.infer_mean(o, 10000, 2)
    print frag_size_est - sum(o) / len(o)


