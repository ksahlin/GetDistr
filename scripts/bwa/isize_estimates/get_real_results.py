'''
Created on Mar 18, 2014

@author: ksahlin
'''
import sys
import re

import matplotlib.pyplot as plt

import numpy as np

def get_bwa_results(bwa_file):
    """
    Gets the average BWA estimation of mean and standard deviation. 
    BWA estimates mean and stddev in batches, a total exact mean 
    can be derived here by linearity. However, to get the exact stddev,
    we need the original samples which we cannot get. The approximation will 
    however be OK if the means does not vary much across batches.
    """

    bwa_output = ''
    for line in open(bwa_file, 'r'):
        bwa_output += line
        #print line

    res = re.findall('[\d\.]+ \+/- [\d\.]+',bwa_output)
    #print res
    res = map(lambda x: tuple(map(lambda x: float(x,), x.split('+/-'))), res)
    #print res
    sum_mean = 0
    sum_stddev = 0

    for mean,stddev in res:
        sum_mean += mean
        sum_stddev += stddev

    mean_final = sum_mean/float(len(res))
    stddev_final = sum_stddev/float(len(res))
    #print mean_final, stddev_final     
    #mean,std_dev = map(lambda x: float(x), res[0].split('+/-'))
    #print 'BWA OUUUUT:', mean,std_dev
    return mean_final,stddev_final

def get_picard_results(picard_file):
    picard_output = ''
    for line in open(picard_file, 'r'):
        picard_output += line
        #print line

    res = re.findall("\t\d+\.\d",picard_output)
    #print res
    mean,std_dev = map(lambda x: float(x), res[:2])
    #print mean_final, stddev_final     
    #mean,std_dev = map(lambda x: float(x), res[0].split('+/-'))
    #print 'BWA OUUUUT:', mean,std_dev
    return mean,std_dev

def get_getdistr_results(getdistr_file):
	getdistr= {}
	for line in open(getdistr_file, 'r'):

		if line[0] != '#' and line[0] != 'N':
			assembly = line.strip()	
		elif line[0] != 'N':
			gd_bwa_mean, gd_bwa_stddev, gd_pic_mean,gd_pic_stddev = map(lambda x: float(x), line.split()[1:])
			# mean = line.split(',')[0][2:-1] 
			# stddev= line.split(',')[1][2:-3] 
			# getdistr[assembly] = (float(mean),float(stddev))
			getdistr[assembly] = (gd_bwa_mean,gd_bwa_stddev,gd_pic_mean,gd_pic_stddev)
	print getdistr
	return getdistr

def plot(bwa,picard,getdistr):
	#fix x-axis
	#a =['reference','ABySS', 'ABySS2', 'Allpaths-LG', 'Bambus2', 'CABOG', 'MSR-CA', 'SGA', 'SOAPdenovo', 'Velvet']
	#a =['reference', 'ABySS', 'ABySS2', 'Allpaths-LG', 'Bambus2', 'CABOG', 'MSR-CA', 'SOAPdenovo', 'Velvet']
	a =['reference', '0','10', '20', '30', '40', '50', '60', '70', '80','90']
	#b = [1,9,5,2,8,4,3,6,7] 
	b= [1,2,3,4,5,6,7,8,9,10,11] # 
	x_axis = dict(zip(b,a))
	print x_axis

	x = [1,2,3,4,5,6,7,8,9,10,11]

	##
	# means
	fig, ax = plt.subplots() 
	bwa_means = map(lambda x: bwa[x_axis[x]][0] ,x_axis)
	picard_means = map(lambda x: picard[x_axis[x]][0] ,x_axis)
	getdistr_bwa_means = map(lambda x: getdistr[x_axis[x]][0] ,x_axis)
	getdistr_picard_means = map(lambda x: getdistr[x_axis[x]][2] ,x_axis)


	# # Find and plot 1st order line of best fit 
	# coeff = np.polyfit( x, bwa_means, 1 ) 
	# p = np.poly1d( coeff ) 
	# temp = np.linspace( 0, 10, 10 ) 
	# ax.plot( temp, p(temp), '-k')# label='Best Fit Line' ) 

	# coeff = np.polyfit( x, picard_means, 1 ) 
	# p = np.poly1d( coeff ) 
	# ax.plot( temp, p(temp), ':k')# label='Best Fit Line' )

	# coeff = np.polyfit( x, getdistr_means, 1 ) 
	# p = np.poly1d( coeff ) 
	# ax.plot( temp, p(temp), '--k')# label='Best Fit Line' )


	ax.plot(x, bwa_means, '-k', marker='o', label='$BWA$')
	ax.plot(x, picard_means, '-k', marker='^', label='$Picard$')
	ax.plot(x, getdistr_bwa_means, '-k', marker='x', label='$Getdistr_{bwa}$')
	ax.plot(x, getdistr_picard_means, '-k', marker='*', label='$Getdistr_{pic}$')
	plt.xticks(x, map(lambda y: x_axis[y],x),fontsize = 8)
	ax.set_ylabel('$\mu$',fontsize=24)
	legend = ax.legend(loc='lower left', shadow=True)
	ax.grid(True)
	plt.savefig("real_data_allpaths_mean.eps", format='eps')

	##
	# stddevs
	fig, ax = plt.subplots() 
	bwa_stddevs = map(lambda x: bwa[x_axis[x]][1] ,x_axis)
	picard_stddevs = map(lambda x: picard[x_axis[x]][1] ,x_axis)
	getdistr_bwa_stddevs = map(lambda x: getdistr[x_axis[x]][1] ,x_axis)
	getdistr_picard_stddevs = map(lambda x: getdistr[x_axis[x]][3] ,x_axis)
	fig, ax = plt.subplots() 
	
	ax.plot(x, bwa_stddevs, '-k', marker='o', label='$BWA$')
	ax.plot(x, picard_stddevs, '-k', marker='^', label='$Picard$')
	ax.plot(x, getdistr_bwa_stddevs, '-k', marker='x', label='$Getdistr_{bwa}$')
	ax.plot(x, getdistr_picard_stddevs, '-k', marker='*', label='$Getdistr_{pic}$')

	plt.xticks(x, map(lambda y: x_axis[y],x),fontsize = 8)
	ax.set_ylabel('$\sigma$',fontsize=24)
	legend = ax.legend(loc='lower left', shadow=True)
	ax.grid(True)
	plt.savefig("real_data_allpaths_sigma.eps", format='eps')

def main():
	picard = {}
	bwa = {}
	# for assembly in ['reference', 'ABySS', 'ABySS2', 'Allpaths-LG', 'Bambus2', 'CABOG', 'MSR-CA', 'SOAPdenovo', 'Velvet']: #['reference','ABySS', 'ABySS2', 'Allpaths-LG', 'Bambus2', 'CABOG', 'MSR-CA', 'SGA', 'SOAPdenovo', 'Velvet']:
	#         mean,stddev = get_bwa_results("cut_N50_assemblies/bwa_results_cut/"+assembly+".bwa.1")
	#         bwa[assembly]=(mean,stddev)
	#         mean,stddev = get_picard_results("cut_N50_assemblies/picard_results_cut/"+assembly+'-picard_ref')
	#         picard[assembly] = (mean,stddev)
	# getdistr = get_getdistr_results("cut_N50_assemblies/getdistr_results_cut/getdistr_both.txt")

	for assembly in ["reference","0","10","20","30","40","50","60","70","80","90"]: # ['reference', 'ABySS', 'ABySS2', 'Allpaths-LG', 'Bambus2', 'CABOG', 'MSR-CA', 'SOAPdenovo', 'Velvet']:
	        mean,stddev = get_bwa_results("Allpaths-LG_cut/bwa/Allapths-LG."+assembly+".bwa.1")
	        bwa[assembly]=(mean,stddev)
	        mean,stddev = get_picard_results("Allpaths-LG_cut/picard/"+assembly+'-picard_ref')
	        picard[assembly] = (mean,stddev)
	getdistr = get_getdistr_results("Allpaths-LG_cut/getdistr-results.txt")

	plot(bwa,picard,getdistr)

if __name__ == '__main__':
	main()




