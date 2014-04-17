'''
Created on Sep 25, 2013

@author: ksahlin
'''

import sys
import os

import argparse
import pysam
import subprocess
import math

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
    """
    Gets the average BWA estimation of mean and standard deviation. 
    BWA estimates mean and stddev in batches, a total exact mean 
    can be derived here by linearity. However, to get the exact stddev,
    we need the original samples which we cannot get. The approximation will 
    however be OK if the means does not vary much across batches.
    """
    #TODO: gereralize this fumction to average over batches of several estimations
    bwa_output = ''
    for line in open(bwa_file, 'r'):
        bwa_output += line
        #print line

    res = re.findall('[\d\.]+ \+/- [\d\.]+',bwa_output)
    mean,std_dev = map(lambda x: float(x), res[0].split('+/-'))
    print 'BWA OUUUUT:', mean,std_dev
    return mean,std_dev

def get_quantile(list_,quantile):
    try:
        quantile_index = int(( len( list_ ) - 1) * quantile)
        return  list_[ quantile_index  ]
    except TypeError:   #<  There were an even amount of values
        # Make sure to type results of math.floor/ceil to int for use in list indices
        ceil = int( math.ceil( quantile_index ) )
        floor = int( math.floor( quantile_index ) )
        return (list_[ ceil ] + list_[ floor ] ) / 2

def MAD(list_,median):
    mad = 0
    absolute_deviations = map(lambda x: abs(median-x),list_)
    absolute_deviations.sort()
    mad = get_quantile(absolute_deviations,0.5)
    return mad



def flatten(lst):
    for elem in lst:
        if type(elem) in (tuple, list):
            for i in flatten(elem):
                yield i
        else:
            yield elem

def getdistr_result(bam_path,references,low, high):
    # Put allowed soft clips to 0, because BWA does not align outside boundaries of the reference.
    # i.e. reeds need to be fully contained (mapped) to the contig in this case.
    #weight_dict = model.get_possible_positions_dict(bam_file, 0, low, high)
    list_of_all_obs = [references[x]['o'] for x in references]
    list_of_all_obs = list(flatten(list_of_all_obs))
    list_of_obs = filter(lambda x: x<= high and x >= low,  list_of_all_obs)
    #print bam_file, references, list_of_obs 
    mean_est,std_dev_est = model.estimate_library_parameters(bam_path,list_of_obs, low, high, soft=0)
    # tot_mean_obs = 0
    # tot_stddev_obs = 0
    # tot_nr_obs = 0
    # mean_naive = 0
    # for reference,info in references.iteritems():
    #     list_of_obs = filter(lambda x: x<= high and x >= low,  info['o']) # to work with the same reads as BWA
    #     if len(list_of_obs) == 0:
    #         continue
    #     reference_length = info['l']
    #     temp_mean_est,temp_std_dev_est = model.estimate_library_parameters(list_of_obs, 100, reference_length, soft=0)
    #     tot_mean_obs += temp_mean_est*len(list_of_obs)
    #     tot_stddev_obs += temp_std_dev_est*len(list_of_obs)
    #     tot_nr_obs += len(list_of_obs)
    #     mean_naive += sum(list_of_obs) 

    # mean_est = tot_mean_obs / float(tot_nr_obs)
    # mean_naive = mean_naive / float(tot_nr_obs)
    # std_dev_est = tot_stddev_obs / float(tot_nr_obs)

    # #print(mean_est,std_dev_est, mean_naive, tot_nr_obs)
    return mean_est,std_dev_est


def get_getdistr_results(bam_path):
    """
    Work with the same reads as BWA:
    " BWA estimates the insert size distribution per 256*1024 read pairs.
    It first collects pairs of reads with both ends mapped with a single-end
    quality 20 or higher and then calculates median (Q2), lower and higher 
    quartile (Q1 and Q3). It estimates the mean and the variance of 
    the insert size distribution from pairs whose insert sizes are within
    interval [Q1-2(Q3-Q1), Q3+2(Q3-Q1)]. 

    (The maximum distance x for a pair 
    considered to be properly paired (SAM flag 0x2) is calculated by 
    solving equation Phi((x-mu)/sigma)=x/L*p0, where mu is the mean, 
    sigma is the standard error of the insert size distribution, 
    L is the length of the genome, p0 is prior of anomalous pair 
    and Phi() is the standard cumulative distribution function. 
    For mapping Illumina short-insert reads to the human genome, 
    x is about 6-7 sigma away from the mean. Quartiles, mean, variance 
    and x will be printed to the standard error output.)
    """

    #list_of_obs = []


    with pysam.Samfile(bam_path, 'rb') as bam_file:
        #contig dictionary
        references = dict(zip(bam_file.references, map(lambda x: {'l':int(x),'o':[]},bam_file.lengths)))
        #print references
        for alignedread in bam_file:
            # only count high quality mappings with both reads mapped. an observation from a read pair once 
            if alignedread.is_read2 and not alignedread.is_unmapped and not alignedread.mate_is_unmapped  and alignedread.mapq >= 20 and alignedread.flag in [147,163]:   
                #print references[bam_file.getrname(alignedread.tid)]['o']
                references[bam_file.getrname(alignedread.tid)]['o'].append( math.fabs(alignedread.tlen))

    all_observations = reduce( list.__add__, map(lambda x: x[1]['o'],references.iteritems()))
    all_observations.sort()

    # # Infer isize distribution with the same reads as BWA:
    # concatenate list of all observations:
    Q1 = get_quantile(all_observations,0.25)
    Q3 = get_quantile(all_observations,0.75)

    bwa_low = Q1-2*(Q3-Q1)
    bwa_high = Q3+2*(Q3-Q1)
    #print bwa_low, bwa_high
    mean_bwa, stddev_bwa = getdistr_result(bam_path,references,bwa_low, bwa_high)

    # Infer isize distribution with the same reads as picard:
    median = get_quantile(all_observations,0.5)
    mad = MAD(all_observations,median)
    #print mad
    mean_picard, stddev_picard = getdistr_result(bam_path,references,median - 10*mad, median + 10*mad)

    #print str(mean_bwa),str(stddev_bwa.real), str(mean_picard),str(stddev_picard.real)
    return str(mean_bwa),str(stddev_bwa.real), str(mean_picard),str(stddev_picard.real)

def print_format(file_, parameters, bwa_means,bwa_stddevs,get_distr_means,get_distr_stddevs,nr_obs_list):
    n = parameters[-1] 
    avg_bwa_mean        = sum(bwa_means)/n
    avg_bwa_stddev       = sum(bwa_stddevs)/n
    avg_getdistr_mean   = sum(get_distr_means)/n
    avg_getdistr_stddev   = sum(get_distr_stddevs)/n
    avg_nr_obs_list        = sum(nr_obs_list)/n

    sample_stddev_bwa_mean = (sum(map(lambda x: (x-avg_bwa_mean)**2, bwa_means))/float(n))**0.5
    sample_stddev_bwa_stddev = (sum(map(lambda x: (x-avg_bwa_stddev )**2, bwa_stddevs))/float(n))**0.5
    sample_stddev_getdistr_mean= (sum(map(lambda x: (x-avg_getdistr_mean )**2, get_distr_means))/float(n))**0.5
    sample_stddev_getdistr_stddev = (sum(map(lambda x: (x-avg_getdistr_stddev)**2, get_distr_stddevs))/float(n))**0.5

    for i,item in enumerate(parameters):
        # if i+1 == len(parameters):
        #     file_.write(str(item)+'\n')
        # else:
        file_.write(str(item)+'\t')

    file_.write(str(int(avg_nr_obs_list))+'\t'+str(int(avg_bwa_mean))+'\t'+str(int(avg_bwa_stddev))+'\t'
        +str(int(avg_getdistr_mean))+'\t'+str(int(avg_getdistr_stddev))+'\t'+str(int(sample_stddev_bwa_mean))
        +'\t'+str(int(sample_stddev_bwa_stddev)) + '\t'+str(int(sample_stddev_getdistr_mean))+ 
        '\t'+str(int(sample_stddev_getdistr_stddev))+'\n')



def main(args):

    sim_out = open(os.path.join(args.outpath,'bwa_vs_getdristr_sim_out'),'w')
    sim_out.write('lib_mean\tlib_std\tcov\tread_len\tcont_len\tnr_trials\tavg_nr_mapped_pairs\t \
        bwa_mean_average\tbwa_stddev_average\tgetdistr_mean_average\tgetdistr_stddev_average\t \
        bwa_sample_std_dev_of_mean\tbwa_sample_std_dev_of_stddev\tgetdistr_sample_stddev_of_mean\t \
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

            getdistr_mean,getdistr_stddev,nr_obs = get_getdistr_results(bam_path, args_object.cont_len)
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





