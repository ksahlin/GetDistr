'''
Created on Sep 18, 2013

@author: ksahlin
'''

import argparse

# class PositionStats:
# 	def get_spanning_reads(self, position, read_length, read_path):
# 	   	observations = [ ]
# 	   	for read in pysam.Samfile( read_path, "rb" ):
# 	   		if not read.is_proper_pair or read.is_read2: # or len( read.cigar ) > 1:
# 	   			continue

# 	       	if read.tlen > 0:
# 	        	if (read.positions[ 0 ] <= position) and ((read.positions[ 0 ] + read.tlen) >= position + read_length):
# 					observations.append( abs( read.tlen ) )
# 	       	else:
# 	           	if (read.positions[ -1 ] >= (position + 1)) and ((read.positions[ - 1 ] + read.tlen) < position):
# 					observations.append( abs( read.tlen ) )

# 		return observations





import sys
import pysam
from scipy.stats import norm
import math

from getdistr import model


class ReadContainer(object):
	"""docstring for ReadContainer"""
	def __init__(self, position):
		super(ReadContainer, self).__init__()
		self.position = position
		self.reads = []
	def __str__(self):
		"Prints reads on fasta format"
		pass
	def add_read(self,read):
		self.reads.append(read)

	def print_bam(self):
		"Print reads on bam format "
		pass

	def calc_observed_insert(self):
		n = len(self.reads)
		if n == 0:
			return 0

		isize_obs = 0
		modulo = 0
		for read1,read2 in zip(self.reads[:-1],self.reads[1:]):
			if modulo == 1:
				modulo = 0
				continue
			else: 
		   		assert read1.tlen == - read2.tlen
		   		assert read1.is_proper_pair == read2.is_proper_pair
		   		isize_obs += abs(read1.tlen)
				modulo = 1

		return isize_obs / float(len(self.reads)/2)

	
	def calc_pvalue(self,expected_insert,sigma):
		n = len(self.reads)
		if n == 0:
			return 0.5 # return a non significant p-value
	    #insert size
		mean_insert_obs = self.calc_observed_insert()
		# tot_insert = 0
		# for read in self.reads:
		# 	tot_insert += abs(read.tlen)
		# mean_insert_obs = tot_insert/float(n)  
		#print mean_insert_obs

		z = (mean_insert_obs - expected_insert)/(float(sigma)/math.sqrt(n)) # z-statistic
		p_value_upper_area = norm.sf(z) 
		#print z, p_value
		return p_value_upper_area


def ParseBAMfile(bamfile,read_length):
	with pysam.Samfile(bamfile, 'rb') as bam:
	# bam = pysam.Samfile(bamfile, 'rb')
	#bam2 = pysam.Samfile(bamfile, 'rb')


		ref_lengths = bam.lengths

		##
		# create a container holding all interesting reads 
		# for a given position in the reference sequence
		counter =0

	   	counter2= 0 

		container = {}
		for i in range(ref_lengths[0]): # assumes only one chromosome, i.e. a single reference strand
			container[i] = ReadContainer(i)
		print ref_lengths[0]
	   	
	   	for read1 in bam:
	   		read2 = next(bam)
			#print read1
			#print read2
			#print read1.tlen,read2.tlen
	   		assert read1.tlen == - read2.tlen
	   		assert read1.is_proper_pair == read2.is_proper_pair
		# Make sure that bamfile is sorted so that mates come after eachotern
		# Throw assertion error here otherwise.
		# use for read1,read2 in zip(bam1,bam2,2) to iterate over the samfiles simultaneosly
		# to use last aligned coordinate of mate instead of read_length (in case of softclipps)!

	   		#print read.tlen,read.is_proper_pair,read.tlen, read.aend,read.pnext
	   		counter += 1
			if counter % 20000 == 0:
				print 'Processed ', counter, 'reads.'
				#break

	   		if not read1.is_proper_pair: # or len( read.cigar ) > 1:
	   			counter2 +=1
	   			continue



			#print inner_start_pos, inner_end_pos
			if read1.tlen > 0:
				inner_start_pos = read1.aend
				inner_end_pos = read2.pos
				for pos in range(inner_start_pos,inner_end_pos):
					container[pos].add_read(read1)
					container[pos].add_read(read2)
			else:
				inner_start_pos = read2.aend
				inner_end_pos = read1.pos

				for pos in range(inner_start_pos,inner_end_pos):
					container[pos].add_read(read1)
					container[pos].add_read(read2)



	print counter2, counter
	return container


class BreakPointContainer(object):
	"""docstring for BreakPointContainer"""
	def __init__(self):
		super(BreakPointContainer, self).__init__()
		self.clusters = []


	def add_bp_to_cluster(self,pos,dist_thresh):
		new_cluster = 1
		for cluster in self.clusters:
			min_pos = min(cluster)
			max_pos = max(cluster)
			print min_pos,max_pos
			if pos <= max_pos + dist_thresh and pos >= min_pos - dist_thresh:
				new_cluster = 0
				cluster.append(pos)
				break

		if new_cluster:
			self.clusters.append([pos])


		if len(self.clusters) == 0:
			self.clusters.append([pos])

	# def __str__(self):
	# 	for cluster in self.clusters:
	# 		print cluster


def GetBasePvalue(container,mu,sigma):
	tot_mean = 0
	nr_obs =0
	est = model.NormalModel(mu,sigma,100,s_inner=0)
	# exp_stddev =  lagg till test av forvantad standard avvikelse
	exp_insert = float(est.expected_mean(0,100000))
	#exp_insert =  mu+ (sigma**2)/float(mu + 1)

	print exp_insert, type(exp_insert)


	##
	# need to select a consensus loci for a breakpoint
	# using only insert size 
	sv_container = BreakPointContainer()


	for bp in range(10000):
		#print container[bp]
		#print container[bp].calc_insert_pvalue
		p_val_upper_area = container[bp].calc_pvalue(exp_insert,sigma)
		#print p_val
		# tot_insert = 0
		# if len(container[bp].reads) > 0:
		# 	for read in container[bp].reads:
		# 		tot_insert += abs(read.tlen)
		# 	tot_mean += tot_insert/len(container[bp].reads)

		obs_mean = container[bp].calc_observed_insert()

		if len(container[bp].reads) > 10:
			nr_obs += 1
			tot_mean += obs_mean


		if p_val_upper_area < 0.001:
			sv_container.add_bp_to_cluster(bp,1)
			print 'Significant large position insert size reads: Pos: ',bp, ' p-value = ', p_val_upper_area, 'nr of reads:', len(container[bp].reads), 'obs insert: ', obs_mean,'inferred insertion in donor here'

		if p_val_upper_area > 0.999:
			sv_container.add_bp_to_cluster(bp,1)
			print 'Significant small position insert size reads: Pos: ',bp, ' p-value = ', p_val_upper_area, 'nr of reads:', len(container[bp].reads), 'obs insert: ', obs_mean,'inferred deletion in donor here'

		

	print 'avg obs_mean:', tot_mean/nr_obs
	print 'Sv clusters:'
	for sv in sv_container.clusters:
		print sv

def main(args):
    container = ParseBAMfile(args.bampath,50)
    GetBasePvalue(container,args.mean,args.stddev)




# def SVdetector(bam1,bam2):
#     def ParseBAMfile(bamfile):
#         with pysam.Samfile(bamfile, 'rb') as bam_file:
#             ref_length = bam_file.lengths[0]
#             position_stats=[]
#             for i in range(0,ref_length):
#                 position_stats.append([0,0,0,0])  # inner list is [covarage , insert size , # insert size observations , nr of split reads]
#             counter =0
#             nr_paired_reads = 0
#             insert_obs = 0
#             insert_obs_sq = 0
#             for alignedread in bam_file:
#                 counter += 1
#                 if not alignedread.is_unmapped and not alignedread.mate_is_unmapped:
#                     #coverage and split reads
#                     pos_mapped = alignedread.positions
#                     prev_pos = pos_mapped[0]
#                     for pos in pos_mapped:
#                         position_stats[pos][0] += 1 #add to coverage
#                         if pos - prev_pos > 1: #gap
#                             position_stats[pos][3] += 1 #add to split reads                            
#                         prev_pos = pos
#                     position_stats[pos_mapped[-1]][3] += 1 #don't forget last position
#                     #insert size
#                     if alignedread.is_read1:                    
#                         i_size = alignedread.tlen
#                         nr_paired_reads += 1
#                         insert_obs += alignedread.tlen
#                         insert_obs_sq += alignedread.tlen**2
#                         inner_start_pos = alignedread.aend
#                         inner_end_pos = alignedread.pnext
#                         for pos in range(inner_start_pos,inner_end_pos):
#                             position_stats[pos][1] += i_size
#                             position_stats[pos][2] += 1
#                 if counter % 5000 == 0:
#                     print 'Processed ', counter, 'positions on the genome'


#         n=float(nr_paired_reads)
#         mean_ins=insert_obs/n
#         std_dev_ins = ((insert_obs_sq - n*mean_ins**2)/(n-1))**0.5
#         tot_cov = 0
#         for pos in range(0,len(position_stats)):
#             tot_cov +=  position_stats[pos][0]
#             #print position_stats[pos][0],tot_cov
#             tot_i_size = position_stats[pos][1]
#             try:
#                 position_stats[pos][1] = tot_i_size/float(position_stats[pos][2])
#             except ZeroDivisionError:
#                 position_stats[pos][1] = 0                


#         mean_cov = tot_cov/5000.0 #replace with float(len(position_stats))
        
#         print 'Nr paired reads total: ',n,'Estimated mean insert size: ', mean_ins,'Estimated std dev insert size: ', std_dev_ins,'Estimated mean coverage: ', mean_cov
#         return(position_stats,mean_cov,std_dev_ins)


#     def GetBasePvalue(stats1,stats2,mean_cov1,std_dev_ins1,mean_cov2,std_dev_ins2,bp):
#         split_obs = stats1[3] - stats2[3]
#         mean_split1 = stats1[0]/float(50) #hardcoded read length
#         mean_split2 = stats2[0]/float(50) #hardcoded read length
#         cov_obs = stats1[0] - stats2[0]
#         try:
#             std_dev_t_test = (std_dev_ins1**2/stats1[2] + std_dev_ins2**2/stats2[2] )**0.5
#             degrees_of_frdm =   (std_dev_ins1**2/stats1[2] + std_dev_ins2**2/stats2[2] )**2 / ( (std_dev_ins1**2/stats1[2])**2 / float(stats1[2] -1 )  + (std_dev_ins2**2/stats2[2])**2  / float(stats2[2] -1 ) )
#             #print 'Std dev: ', std_dev_t_test, 'dgrs of frdm: ', degrees_of_frdm, stats1[2], stats2[2]
#         except ZeroDivisionError:
#             #just put some values not giving  a significant t-statistic
#             std_dev_t_test = 1000
#             degrees_of_frdm = 2
#         ins_obs = (stats1[1] - stats2[1]) / std_dev_t_test   #assuming equal standard deviation here, therefore I only use one!
#         #print mean_split1,mean_split2,cov_obs, ins_obs
#         pval1 = skellam.cdf(split_obs,mean_split1 ,mean_split2, loc=0)
#         pval2 = skellam.cdf(cov_obs, mean_cov1 ,mean_cov2, loc=0)
#         pval3 = (1- t.cdf(ins_obs, degrees_of_frdm, loc=0, scale = 1) ) #this only does one sided testing
#         #print pval3
#         if pval1 < 0.01:
#             print 'Significant position split reads: Pos: ',bp, ' p-value = ', pval1
#         if pval2 < 0.01:
#             print 'Significant position covarage diff: Pos: ',bp, ' p-value = ', pval2
#         if pval3 < 0.01:
#             print 'Significant position ins size diff: Pos: ',bp, ' p-value = ', pval3                        
        
#         return(pval1,pval2,pval3)
    
#     #Read in libraries from the two donor genomes
#     list1, mean_cov1,std_dev_ins1 = ParseBAMfile(bam1)
#     list2, mean_cov2,std_dev_ins2 = ParseBAMfile(bam2)
#     #do statistical testing for each bp in the reference genome
#     vector_of_pvals1=[]
#     vector_of_pvals2=[]
#     vector_of_pvals3=[]
#     for bp in range(len(list1)):
#         pval1,pval2,pval3 = GetBasePvalue(list1[bp],list2[bp],mean_cov1,std_dev_ins1,mean_cov2,std_dev_ins2,bp)
#         vector_of_pvals1.append( pval1)
#         vector_of_pvals2.append( pval2)
#         vector_of_pvals3.append( pval3)

#     return()


if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Infer variants with simple p-value test using theory of GetDistr - proof of concept.")
	parser.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	parser.add_argument('mean', type=int, help='mean insert size. ')
	parser.add_argument('stddev', type=int, help='Standard deviation of insert size ')


	args = parser.parse_args()
	main(args)



        
