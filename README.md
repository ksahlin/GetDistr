GetDistr
========

This is a reference implementation of the theory presented in a preprint version on bioRxiv [GetDistr](http://biorxiv.org/content/biorxiv/early/2015/08/04/023929.full.pdf) and a hevily modified preprint submitted to RECOMB 2016. This implementation contains various scripts that we used for plotting of data and inference and testing. 


### What this repository offers
The data considered here is paired reads (paired end reads, mate pair reads, fosmid/BAC ends) --- any read pair that has and "insert size", that is, an unknown distance between the reads. 


* Calculates MLE gaps between contigs given a set of observations in module (getdistr/model.py). This is an extended version of [GapEst](http://www.ncbi.nlm.nih.gov/pubmed/22923455) to also calculate better estimates for libraries that better fit a logNormal distribution (more general) rather than normal.

* Calculates an estimate of the true full library distribution given a sorted bam file of aligned reads to contigs. That is, for a read pair library aligned on a "fragmented" assembly (small contigs compared to library insert size), smaller insert size reads are more frequently observed than larger ones when looking at reads with both mates mapping within a contig. This biases the estimation of the full library. There is a way to model this (presented in [GetDistr](http://biorxiv.org/content/biorxiv/early/2015/08/04/023929.full.pdf) ) this module helps to infer the correct library distribution.

* Calculates the observed fragment length distribution over a position or region on a genome for detection of indels. As a consequence, a better null-hypothesis can be formulated for toosl such as [CLEVER](http://bioinformatics.oxfordjournals.org/content/early/2012/10/11/bioinformatics.bts566.abstract) 


Using GetDistr
-------------------------

This repostory's main purpose is to be a reference implementation (for now, but this will soon be extended). However, there are two modules that could be used for general analysis:

#### "assemblymodule"

This module is already pretty mature and consists of several modules. The module will generate lots of statistics and plost from the aligned reads of a read pair library <em>read pair s are assumed to be aligned in fr orientation!</em>. Is assumes paired reads are aligned to ether a reference genome (calling variants), or an genome assembly (calling misassemblies, sort of..).

	$ cd <path to GetDistr>/getdistr/assemblymodule/
	$ python main.py -h
	usage: main.py [-h]
	               {pipeline,get_bp_stats,cluster_pvals,filter,lib_est,get_gaps}
	               ...

	positional arguments:
	  {pipeline,get_bp_stats,cluster_pvals,filter,lib_est,get_gaps}
	                        help for subcommand
	    pipeline            Run the entire pipeline
	    filter              Filters bam file for better uniform coverage.
	    lib_est             Estimate library parameters
	    get_bp_stats        Scan bam file and calculate pvalues for each base pair
	    get_gaps            Contig assembly file for gap positions
	    cluster_pvals       Takes a pvalue file and clusters them into significan
	                        regions

	optional arguments:
	  -h, --help            show this help message and exit

For instance, the whole pipline could be run to create various library statistics, calculate p-values over every single base pair in the genome and even call varinats based on acculupation of signifivcat p-values. 

	$ python main.py pipeline -h
	usage: main.py pipeline [-h] [--n N] [--lib_min LIB_MIN] [--lib_max LIB_MAX]
	                        [--plots]
	                        bampath assembly_file outfolder

	positional arguments:
	  bampath            bam file with mapped reads.
	  assembly_file      Fasta file with assembly/genome.
	  outfolder          Outfolder.

	optional arguments:
	  -h, --help         show this help message and exit
	  --n N              Neighborhood size.
	  --lib_min LIB_MIN  Minimum insert size (if in doubt, just set lib_min =
	                     2*read_length).
	  --lib_max LIB_MAX  Maximum insert size (tight bound is not necessary, choose
	                     a larger value rather than smaller).
	  --plots            Plot pval distribution.


#### model.py 

The code in this module can be copy pasted to other code, impoted or run interactivily. Example:

	$ python
	>>>import model
	>>> library = NormalModel()
	>>> r=100
	>>> s=35
	>>> mu = 500
	>>> sigma = 100
	>>> my_library = NormalModel(mu,sigma,r,s)

This is a read pair library object containing the attrubutes of the library, read length, max allowed softclipps on eather read, mu and sigma of library. We can for instance calculate MLE of expactes mean fragment length given two contig sizes and and a gap

	>>> c1=1000
	>>> c2=2000
	>>> gap=0
	>>> ML_gap  = float(my_library.run_GapEst([245,190,237,270],c1,c2))
	>>> ML_gap
	419.5628356933594
	>>> sigma_p =float(my_library.expected_standard_deviation(gap,c1,c2))
	>>> sigma_p
	86.3915250134528

Or we can infer a gap given som observations from the library

	>>> ML_gap  = float(my_library.run_GapEst([245,190,237,270],c1,c2))
	>>> ML_gap
	419.5628356933594

We can also get a gap estimate if asuming another coverage distribution (not poisson) with the more general NegBin as:
	
	>>> precision, mean_cov = 1, 20
	>>> ML_gap = my_library.infer_mean_fast([245,190,237,270], c1, precision, b=c2, coverage = 100, coverage_model = "NegBin")
	>>> ML_gap
	404


