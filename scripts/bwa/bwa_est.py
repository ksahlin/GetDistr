'''
Created on Sep 25, 2013

@author: ksahlin
'''

import sys

import argparse
import pysam

import model

# model.estimate_library_mean(list_of_obs, r, a, s=None)


if __name__ == '__main__':

    list_of_obs = []

    with pysam.Samfile(sys.argv[1], 'rb') as bam_file:
        for alignedread in bam_file:
            if alignedread.is_read1:  # only count an observation from a read pair once 
                list_of_obs.append(alignedread.tlen)



    mean_est = model.estimate_library_mean(list_of_obs, 100, 15000, soft=70)
    print(mean_est)




