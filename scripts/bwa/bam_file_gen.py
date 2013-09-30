'''
Created on Sep 30, 2013

@author: ksahlin
'''


import sys
import os

import argparse
import random
import tempfile
import subprocess
import pysam

from simulate import simulate

def sam_to_bam(sam_path, bam_path):
    sam_file = pysam.Samfile(sam_path, "r")
    bam_file = pysam.Samfile(bam_path, "wb", template=sam_file)

    for alignment in sam_file:
        bam_file.write(alignment)

def get_reference(length):
    ref = ''.join(random.choice('AGCT') for i in range(length))
    return(ref)


def get_contig(ref, c_len):
    return(ref[ len(ref) / 2 - c_len : len(ref) / 2 + c_len])


def map(pe1_path, pe2_path, contig_path, args, threads=8):

    print 'Aligning with bwa mem'
    work_dir = tempfile.mkdtemp()
    genome_db = os.path.join(work_dir, "genome")
    bwa_output = os.path.join(work_dir, "output.sam")

    null = open("/dev/null")
    std_err = open(os.path.join(args.outpath, 'bwa_out'), 'w')
    subprocess.check_call([ "bwa", "index", "-p", genome_db, contig_path ], stderr=null)
    with open(bwa_output, "w") as bwa_file:
        subprocess.check_call([ "bwa", "mem", "-t", str(threads),
                                genome_db, pe1_path, pe2_path ],
                              stdout=bwa_file,
                              stderr=std_err)


    sam_to_bam(bwa_output, os.path.join(args.outpath, str(args.c_len) + '-' + str(args.std_dev)) + ".bam")
    #pysam.sort(bwa_output + ".bam", output_path)
    #pysam.index(output_path + '.bam')





def main():
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
    parser.add_argument('genome', type=str, help='Name of the reference sequence. ')
    parser.add_argument('genomelen', type=int, help='Length of the reference sequence. ')
    parser.add_argument('mean', type=int, help='mean insert size. ')
    parser.add_argument('std_dev', type=int, help='Standard deviation of insert size ')
    parser.add_argument('read_length', type=int, help='Length of read fragments within a pair. ')
    parser.add_argument('coverage', type=int, help='Coverage for read library. ')
    parser.add_argument('outpath', type=str, help='Path to output location. ')
    parser.add_argument('c_len', type=int, help='Contig length. ')

    args = parser.parse_args()


    ref = get_reference(args.genomelen)
    contig = get_contig(ref, args.c_len)

    # ref file
    ref_path = '/tmp/ref.fa'
    reffile = open(ref_path, 'w')
    print >> reffile, '>ref\n' + ref
    reffile.close()

    # contig file
    c_path = '/tmp/ctg.fa'
    c_file = open('/tmp/ctg.fa', 'w')
    print >> c_file, '>ctg\n' + ref
    c_file.close()



    simulate.pe_reads_haploid(args.outpath, args.mean, args.coverage, args.std_dev, args.read_length, ref_path)

    map(os.path.join(args.outpath, 'PE_1'), os.path.join(args.outpath, 'PE_2'), c_path, args)



if __name__ == '__main__':

    main()
