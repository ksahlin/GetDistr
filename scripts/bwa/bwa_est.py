'''
Created on Sep 25, 2013

@author: ksahlin
'''

import sys
import os

import argparse
import pysam
import subprocess

from getdistr import model
import bam_file_gen
import re

# model.estimate_library_mean(list_of_obs, r, a, s=None)

class Args(object):
    def __init__(self, in_values, outpath, genomelen):
        self.lib_mean, self.lib_std, self.coverage, self.read_length, self.cont_len, self.nr_experiments = in_values
        self.outpath = outpath
        self.genomelen = genomelen


##
# testing the yeild function
def read_input(sim_file):
    for line in sim_file:
        try:
            in_values = map(lambda x: int(x), line.strip().split())
        except:
            continue
        yield in_values

def create_bam(args):
    return(bam_file_gen.main(args))


def get_bwa_results(bwa_file):
    bwa_output = ''
    for line in open(bwa_file, 'r'):
        bwa_output += line
        #print line

    res = re.findall('[\d\.]+ \+/- [\d\.]+',bwa_output)
    mean,std_dev = map(lambda x: float(x), res[0].split('+/-'))
    print 'BWA OUUUUT:', mean,std_dev
    return mean,std_dev


def get_getdistr_results(bam_path, args):

    list_of_obs = []
    list_of_isize = []
    with pysam.Samfile(bam_path, 'rb') as bam_file:
        for alignedread in bam_file:
            if alignedread.is_read1 and not alignedread.is_unmapped and not alignedread.mate_is_unmapped :  # only count an observation from a read pair once 
                list_of_obs.append(alignedread.tlen)
                #if len(alignedread.cigar) > 1:
                print alignedread.cigarstring

    # Put allowed soft clips to 0, because BWA does not align outside boundaries of the reference.
    # i.e. reeds need to be fully contained (mapped) to the contig in this case.
    mean_est,std_dev_est = model.estimate_library_parameters(list_of_obs, 100, args.cont_len, soft=0)
    mean_naive = sum(list_of_obs) / float(len(list_of_obs))
    print(mean_est,std_dev_est, mean_naive, len(list_of_obs))
    print list_of_obs
    return mean_est,std_dev_est, len(list_of_obs)

def print_format(file_, parameters, bwa_means,bwa_stddevs,get_distr_means,get_distr_stddevs,nr_obs_list):
    n = parameters[-1] 
    avg_bwa_mean        = sum(bwa_means)/n
    avg_bwa_stddev       = sum(bwa_stddevs)/n
    avg_getdistr_mean   = sum(get_distr_means)/n
    avg_getdistr_stddev   = sum(get_distr_stddevs)/n
    avg_nr_obs_list        = sum(nr_obs_list)/n

    sample_stddev_bwa_mean = (sum(map(lambda x: (x-avg_bwa_mean)**2, bwa_means))/float(n))**0.5

    for i,item in enumerate(parameters):
        # if i+1 == len(parameters):
        #     file_.write(str(item)+'\n')
        # else:
        file_.write(str(item)+'\t')

    file_.write(str(nr_obs_list)+'\t'+str(avg_bwa_mean)+'\t'+str(avg_bwa_stddev)+'\t'
        +str(avg_getdistr_mean)+'\t'+str(avg_getdistr_stddev)+'\t'+str(sample_stddev_bwa_mean)+'\n')



def main(args):

    sim_out = open(os.path.join(args.outpath,'bwa_vs_getdristr_sim_out'),'w')
    sim_out.write('lib_mean\tlib_std\tcov\tread_len\tcont_len\tnr_trials\tavg_nr_mapped_pairs\t \
        bwa_mean_average\tbwa_stddev_average\tbwa_sample_std_dev_of_mean\tbwa_sample_std_dev_of_stddev\t \
        getdistr_mean_average\tgetdistr_stddev_average\tgetdistr_sample_stddev_of_mean\t \
        getdistr_sample_stddev_of_stddev\n')
    experiment_folder = os.path.join(args.outpath,'experiment_files')
    if not os.path.exists(experiment_folder):
        os.makedirs(experiment_folder)

    for in_values in read_input(open(args.sim_file, 'r')):
        successful_experiments = 1
        args_object = Args(in_values, args.outpath, args.genomelen)
        bwa_means = []
        bwa_stddevs = []
        get_distr_means = []
        get_distr_stddevs = []
        nr_obs_list = []

        while successful_experiments <= args_object.nr_experiments: # for exp in xrange(args.experiments):
            print 'Processing experiment ' + str(successful_experiments)
            print 'With setting: ', in_values
            try:
                bam_path = create_bam(args_object)
            except subprocess.CalledProcessError:
                continue

            bwa_mean, bwa_stddev = get_bwa_results(os.path.join(args_object.outpath, 'bwa_out'))
            #print bam_path
            bwa_means.append(bwa_mean)
            bwa_stddevs.append(bwa_stddev)

            getdistr_mean,getdistr_stddev,nr_obs = get_getdistr_results(bam_path, args_object)
            get_distr_means.append(getdistr_mean)
            get_distr_stddevs.append(getdistr_stddev)
            nr_obs_list.append(nr_obs)

            bam_file_name = os.path.split(bam_path)[-1]
            os.rename(bam_path,os.path.join(args_object.outpath,experiment_folder,bam_file_name+'-exp'+str(successful_experiments)+'.bam'))
            successful_experiments += 1

        experiment_parameters = [args_object.lib_mean, args_object.lib_std, args_object.coverage, args_object.read_length, args_object.cont_len, args_object.nr_experiments]
        print_format(sim_out, experiment_parameters, bwa_means,bwa_stddevs,get_distr_means,get_distr_stddevs,nr_obs_list)


if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
    parser.add_argument('sim_file', type=str, help='Main simulation file. ')
    #parser.add_argument('genome', type=str, help='Name of the reference sequence. ')
    #parser.add_argument('mean', type=int, help='mean insert size. ')
    #parser.add_argument('std_dev', type=int, help='Standard deviation of insert size ')
    #parser.add_argument('read_length', type=int, help='Length of read fragments within a pair. ')
    #parser.add_argument('coverage', type=int, help='Coverage for read library. ')
    parser.add_argument('outpath', type=str, help='Path to output location. ')
    #parser.add_argument('experiments', type=int, help='Number of experiment for each line in sim_in.txt file to run. ')
    parser.add_argument('genomelen', type=int, help='Length of the reference sequence. ')
    #parser.add_argument('c_len', type=int, help='Contig length. ')

    args = parser.parse_args()
    main(args)





