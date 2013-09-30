'''
Created on Sep 25, 2013

@author: ksahlin
'''

import sys

import argparse

import model


# frag_obj = model.NormalModel(mean, stddev, read length, softclipped)
# frag_obj.infer_mean(list of obs, reference length, bp-precision)

# Infer clever inserts
for line in open(sys.argv[1]):
    o = list(line.split())
    o = [int(i) for i in o]
    print o
    frag_obj = model.NormalModel(500, 50, 100, 2)
    frag_size_est = frag_obj.infer_mean(o, 10000, 2)
    print frag_size_est - sum(o) / len(o)


