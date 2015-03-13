"""
	Instructions:
	1. Install python 2.6.6. Seriously, come on.. Do it! Easy peacy with pyenv, see instructions e.g. at:
		https://github.com/ksahlin/eval_kmer_choice
		then install pysam for python2.6 with pip (in the 2.6.6 environment of course)
	2. Place this script in "mrfstructvar" folder if not already there (i.e. where setup.py is)
	3. make sure to configure you pythonpath in bash_profile/bashrc as:
	 PYTHONPATH=$PYTHONPATH:[MoDIL_Home]:[MoDIL_Home]/modules 
	 where [MoDIL_Home] is the path is 
	 "/proj/b2013072/private/svest_evaluation/tools_src/structural_variations/mrfstructvar"
	 in our case.
	4. Now your ready to go! Make sure to specify complete absolute path to the argument "outfolder".

	Results are found in: outfolder/DATA/POST_PROCESSING/[chrQ]/chrQ.INDEL.txt   (where Q \in(1,22,X,Y)) (I remapped the organisms chr names to the humans)
"""
import os
import sys
import argparse
import pysam
import re
import subprocess
import shutil

# Yes, a global..
EXE_DIR_USER = os.path.dirname(os.path.abspath(__file__))

def is_proper_aligned_unique_innie(read):
    mapped_ok = not (read.is_unmapped or read.mate_is_unmapped)
    orientation_ok = (read.is_reverse and not read.mate_is_reverse and read.tlen < 0) or \
            (not read.is_reverse and read.mate_is_reverse and read.tlen > 0)
    quality_ok = read.mapq > 10 and not read.is_secondary
    same_ref = read.rname == read.mrnm

    return mapped_ok and orientation_ok and quality_ok and same_ref

def format_converter(infile,outfolder,args):
	nr_chr = len(infile.references)
	ref_names= map(lambda x: 'chr'+str(x) , range(1,nr_chr+1))
	reference_dict = dict(zip(infile.references, ref_names))
	file_dict={}
	for ref in reference_dict:
		file_dict[reference_dict[ref]] = open(os.path.join(outfolder,'INTERESTING_DB/'+reference_dict[ref]), 'w')
	print file_dict

	for i,read in enumerate(infile):
			## if proper pair
		if is_proper_aligned_unique_innie(read) and not read.is_reverse:
			assert read.tlen > 0
			ref = infile.getrname(read.rname)
			modil_stat =  args.mean - (read.mpos - (read.pos + args.readlength) - 1 + 2*args.readlength) 
			print >>file_dict[reference_dict[ref]], "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(read.pos + args.readlength, read.mpos, read.qname, modil_stat, args.mean, i  )

	print 'done converting data'
	return reference_dict

def set_data_parameters(args):

	print EXE_DIR_USER
	ROOT=EXE_DIR_USER.split("/structural_variations/mrfstructvar")[0]
	print ROOT

	lines=''
	for line in open("mrfstructvar.properties",'r').readlines():
		line = re.sub(r'DATA_DIR_USER=.+',r'DATA_DIR_USER="{0}"'.format(args.outfolder), line)
		line = re.sub(r'EXE_DIR_USER=.+',r'EXE_DIR_USER="{0}"'.format(ROOT), line)
		line = re.sub(r'RESULTS_DIR=.+',r'RESULTS_DIR=DATA_DIR_USER + "/DATA"', line)
		line = re.sub(r'INPUT_DATA_DIR=.+',r'INPUT_DATA_DIR=DATA_DIR_USER +"/INTERESTING_DB"', line)
		line = re.sub(r'MEAN_INSERT_SIZE =.+',r'MEAN_INSERT_SIZE = {0}'.format(args.mean), line)
		line = re.sub(r'STD_INSERT_SIZE =.+',r'STD_INSERT_SIZE = {0}'.format(args.sd), line)
		line = re.sub(r'READ_LENGTH =.+',r'READ_LENGTH = {0}'.format(args.readlength), line)

		lines += line
		#print line,

	print >> open("mrfstructvar.properties",'w'), lines 
	
	lines=''
	for line in open("MoDIL_simple.py",'r').readlines():
		line = re.sub(r'MODIL_HOME=.+',r'MODIL_HOME="{0}"'.format(EXE_DIR_USER), line)
		lines += line
	print >> open("MoDIL_simple.py",'w'), lines

	lines=''
	for line in open("MPProperties.py",'r').readlines():
		line = re.sub(r'MODIL_HOME=.+',r'MODIL_HOME="{0}"'.format(EXE_DIR_USER), line)
		lines += line
	print >> open("MPProperties.py",'w'), lines


def run_setup(args):
	stdout_file = open(os.path.join(args.outfolder,'setup0.stdout'),'w')
	stderr_file = open(os.path.join(args.outfolder,'setup0.stderr'), 'w')
	p = subprocess.check_call(["python", "setup.py", "0"], stderr=stderr_file,stdout=stdout_file)
	# output, err = p.communicate()
	# p.kill()
	#print output
	#print err
	#p = subprocess.check_call(["python","test/modil_demo.py"])
	#output, err = p.communicate()

def run_modil(reference_dict,args):
	print 'RUNNING setup.py 1'
	stdout_file = open(os.path.join(args.outfolder,'setup1.stdout'),'w')
	stderr_file = open(os.path.join(args.outfolder,'setup1.stderr'), 'w')

	p = subprocess.check_call(["python", "setup.py", "1"],stderr=stderr_file, stdout=stdout_file) 
	# output, err = p.communicate()

	# print output
	# print err
	print 'STARTING ACTUAL MODIL'
	steps=args.genome_length/args.stepsize
	for ref in reference_dict:
		for step in range(steps):
			stdout_file = open(os.path.join(args.outfolder,'{0}.stdout').format(ref+'_'+str(step)),'w')
			stderr_file = open(os.path.join(args.outfolder,'{0}.stderr').format(ref+'_'+str(step)), 'w')
			p = subprocess.check_call(["python", EXE_DIR_USER+"/MoDIL_simple.py", \
									"{0}".format(reference_dict[ref]), "{0}".format(step),\
									 "{0}".format(args.stepsize)], stderr=stderr_file, stdout=stdout_file) 
			# output, err = p.communicate()
			# print output
			# print err
			print 'DONE WITH MoDIL_simple.py STEP', step
			p = subprocess.Popen(["python", EXE_DIR_USER+"/qsub_MoDIL_inference.py",\
									"{0}".format(reference_dict[ref]), "{0}".format(step)],\
									stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
			out,err = p.communicate()
			print out
			print err
			# p = subprocess.check_call(["python", EXE_DIR_USER+"/qsub_MoDIL_inference.py",\
			# 						 "{0}".format(reference_dict[ref]), "{0}".format(step)],\
			# 						 stderr=stderr_file, stdout=stdout_file) 
			# output, err = p.communicate()
			# print output
			# print err
			print 'DONE WITH qsub_MoDIL_inference.py STEP', step
		p = subprocess.check_call(["python", EXE_DIR_USER+"/post_processing/modil_post_processing.py",\
									 "{0}".format(reference_dict[ref])], \
									 stderr=stderr_file, stdout=stdout_file) 
		# output, err = p.communicate()
		# print output
		# print err
		print 'DONE WITH modil_post_processing.py chr:', reference_dict[ref]

def main(args):
	set_data_parameters(args)
	run_setup(args)
	reference_dict = format_converter(pysam.Samfile(args.bam, 'rb'), args.outfolder, args)
	run_modil(reference_dict,args)

if __name__ == '__main__':

	parser = argparse.ArgumentParser()#prog="Infer variants with simple p-value test using theory of GetDistr - proof of concept.")
	parser.add_argument('bam', type=str, help='bam file. ')
	parser.add_argument('outfolder', type=str, help='outfolder path NEEDS TO BE FULL PATH with no slash at end!. e.g /Users/ksahlin/_tmp/MODIL/modil_out ')
	parser.add_argument('mean', type=int, help='lib_mean')
	parser.add_argument('sd', type=int, help='lib_sd')
	parser.add_argument('readlength', type=int, help='read length')
	parser.add_argument('ess_ratio', type=float, help='Corrected sample size')
	parser.add_argument('genome_length', type=int, help='genome length')
	parser.add_argument('step_size', type=int, default=100000, help='stepsize in modil (something related to speed/memory). set to smaller than genome size')



	#TODO: Inject ess_ratio in model code somewhere!
	args = parser.parse_args()
	if os.path.exists(args.outfolder):
		shutil.rmtree(args.outfolder)
	if not os.path.exists(args.outfolder):
		os.makedirs(args.outfolder)

	main( args )



