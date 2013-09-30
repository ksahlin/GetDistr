'''
Created on Sep 25, 2013

@author: ksahlin
'''

import sys
import os

import argparse
import pysam

import model
import bam_file_gen

# model.estimate_library_mean(list_of_obs, r, a, s=None)

def read_input(sim_file):
    for line in sim_file:
        in_values = map(lambda x: int(x), line.strip().split())
        lib_mean, lib_std, cov, read_len, cont_len, nr_trials = in_values

def create_bam(args):
    bam_file_gen.main(args)

def get_bwa_results():
    pass

def get_getdistr_results(args):

    list_of_obs = []

    with pysam.Samfile(sys.argv[1], 'rb') as bam_file:
        for alignedread in bam_file:
            if alignedread.is_read1 and not alignedread.is_unmapped and not alignedread.mate_is_unmapped :  # only count an observation from a read pair once 
                list_of_obs.append(alignedread.tlen)

    mean_est = model.estimate_library_mean(list_of_obs, 100, int(sys.argv[2]), soft=80)
    mean_naive = sum(list_of_obs) / float(len(list_of_obs))
    print(mean_est, mean_naive)


def main(args):
        create_bam(args)
        get_bwa_results(os.path.join(args.outpath, 'bwa_out'))

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
    parser.add_argument('genome', type=str, help='Name of the reference sequence. ')
    #parser.add_argument('genomelen', type=int, help='Length of the reference sequence. ')
    #parser.add_argument('mean', type=int, help='mean insert size. ')
    #parser.add_argument('std_dev', type=int, help='Standard deviation of insert size ')
    parser.add_argument('read_length', type=int, help='Length of read fragments within a pair. ')
    #parser.add_argument('coverage', type=int, help='Coverage for read library. ')
    parser.add_argument('outpath', type=str, help='Path to output location. ')
    #parser.add_argument('c_len', type=int, help='Contig length. ')

    args = parser.parse_args()
    main(args)





