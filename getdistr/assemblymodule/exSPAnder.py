import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import argparse
from genomics_tools.file_formats import bam
from genomics_tools.simulate import genome,contigs,reads
from mapping import align

class Point(object):
    """docstring for Point"""
    def __init__(self, x, y):
        super(Point, self).__init__()
        self.x = x
        self.y = y

class Link(object):
    """docstring for Link"""
    def __init__(self, e,e_prime):
        super(Link, self).__init__()
        self.e = e
        self.e_prime = e_prime
        self.points = []

    def add_point(self):
        pass

    def density(self):
        #e, e_prime
        pass
    def expected(self):
        #e, e_prime
        pass

class Linkings(object):
    """docstring for Linkings"""
    def __init__(self):
        super(Linkings, self).__init__()
        self.linkings = {} # dict (e,e'):  
    def add_linking(self):
        pass

class Path(object):
    """docstring for Path"""
    def __init__(self,extension):
        super(Path, self).__init__()
        self.path = []
        self.linkings = []
        self.e = extension

    def add_linking(self):
        pass


    def score(self, e):
        for linking in self.paths:
            pass

    def support(self, e):

        pass
    def expected(self,path, e):
        pass

class ExSPAnder(object):
    """docstring for ExSPAnder"""
    def __init__(self, paths):
        super(ExSPAnder, self).__init__()
        self.paths = []

    def score_paths(self,e):
        for path in self.paths:
            path.score(e)





def plot():
    #-- Generate some data ----------------------------------------------------
    nx = 10
    x = np.linspace(0, 2*np.pi, 10)
    y = 2 * np.sin(x)

    groups = [('GroupA', (x[0], x[nx//3])),
              ('GroupB', (x[-2*nx//3], x[2*nx//3])),
              ('GroupC', (x[-nx//3], x[-1]))]

    #-- Plot the results ------------------------------------------------------
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Give ourselves a bit more room at the bottom
    plt.subplots_adjust(bottom=0.2)

    ax.plot(x,y, 'k^')

    # Drop the bottom spine by 40 pts
    ax.spines['bottom'].set_position(('outward', 40))

    # Make a second bottom spine in the position of the original bottom spine
    make_second_bottom_spine(label='Treatment')

    # Annotate the groups
    for name, xspan in groups:
        annotate_group(name, xspan)

    plt.xlabel('Dose')
    plt.ylabel('Response')
    plt.title('Experimental Data')

    plt.show()


def annotate_group(name, xspan, ax=None):
    """Annotates a span of the x-axis"""
    def annotate(ax, name, left, right, y, pad):
        arrow = ax.annotate(name,
                xy=(left, y), xycoords='data',
                xytext=(right, y-pad), textcoords='data',
                annotation_clip=False, verticalalignment='top',
                horizontalalignment='center', linespacing=2.0,
                arrowprops=dict(arrowstyle='-', shrinkA=0, shrinkB=0,
                        connectionstyle='angle,angleB=90,angleA=0,rad=5')
                )
        return arrow
    if ax is None:
        ax = plt.gca()
    ymin = ax.get_ylim()[0]
    ypad = 0.01 * np.ptp(ax.get_ylim())
    xcenter = np.mean(xspan)
    left_arrow = annotate(ax, name, xspan[0], xcenter, ymin, ypad)
    right_arrow = annotate(ax, name, xspan[1], xcenter, ymin, ypad)
    return left_arrow, right_arrow

def make_second_bottom_spine(ax=None, label=None, offset=0, labeloffset=20):
    """Makes a second bottom spine"""
    if ax is None:
        ax = plt.gca()
    second_bottom = mpl.spines.Spine(ax, 'bottom', ax.spines['bottom']._path)
    second_bottom.set_position(('outward', offset))
    ax.spines['second_bottom'] = second_bottom

    if label is not None:
        # Make a new xlabel
        ax.annotate(label, 
                xy=(0.5, 0), xycoords='axes fraction', 
                xytext=(0, -labeloffset), textcoords='offset points', 
                verticalalignment='top', horizontalalignment='center')



def main(args):
    if args.simulate:
        print 'here'
        #genome
        g = genome.Genome([0.25]*4,20000,'genome1')
        g.genome()
        print >> open('/tmp/genome.fa','w'), g.genome_fasta_format()
        #contigs
        ctg1, ctg2, ctg3 = g.sequence[8000:10000], g.sequence[10000:12000], g.sequence[12000:13000]
        ctgs = open('/tmp/ctgs.fa','w')
        for i,ctg in enumerate([ctg1,ctg2,ctg3]):
            ctgs.write('>ctg{0}\n{1}\n'.format(i,ctg))
        #reads
        lib = reads.DNAseq(100,100, 3000,1500)
        lib.simulate_mp_reads(g)
        reads1 = open('/tmp/reads1.fa','w')
        reads2 = open('/tmp/reads2.fa','w')
        i=0
        for read in lib.fastq_format():
            if i%2==0:
                reads1.write(read)
            else:
                reads2.write(read)

        #mapping
        align.map_paired_reads('/tmp/reads1.fa', '/tmp/reads1.fa', genome_path, output_path, args)

    #plot()
    bam_parser = bam.BamParser(args.bam_file)
    for unique_link in bam_parser.unique_reads_on_different_references('bwa'):
        print unique_link

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
    parser.add_argument('bam_file', type=str, help='Path to bam file. ')
    parser.add_argument('contig_file', type=str, help='Path to contig fasta file ')
    parser.add_argument('--simulate', action="store_true", help='mean insert size. ')
    #parser.add_argument('std_dev', type=int, help='Standard deviation of insert size ')
    #parser.add_argument('read_length', type=int, help='Length of read fragments within a pair. ')
    #parser.add_argument('coverage', type=int, help='Coverage for read library. ')
    #parser.add_argument('outpath', type=str, help='Path to output location. ')
    #parser.add_argument('experiments', type=int, help='Number of experiment for each line in sim_in.txt file to run. ')
    #parser.add_argument('genomelen', type=int, help='Length of the reference sequence. ')
    #parser.add_argument('c_len', type=int, help='Contig length. ')

    args = parser.parse_args()
    main(args)