from genomics_tools.file_format import bam

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




def main(args):

    for unique_link in bam.unique_links(args.bam_file):
        pass

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