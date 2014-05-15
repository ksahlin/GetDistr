'''
Created on Mar 18, 2014

@author: ksahlin
'''
import sys

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter, FixedLocator
import matplotlib.pyplot as plt
import numpy as np

#from operator import itemgetter


def randrange(n, vmin, vmax):
    return (vmax-vmin)*np.random.rand(n) + vmin

def plot_sigma(a,s, sigma_bwa,sigma_gd):
	lines_1 = ["--", "_", "-.", ":"]
	col_1 = (0.6, 0.6, 0.6)
	lines_2 = ["--", "_", "-.", ":"]
	col_2 = (0.3,0.3,0.3)
	fig, ax = plt.subplots()

	ax.set_xlim(500,3000)

	bwa = map(lambda x: x[0]-x[1], zip(sigma_bwa,s))
	gd = map(lambda x: x[0]-x[1], zip(sigma_gd,s))
	#ax.plot([0] + s[:5] + [3000], [0,0,0,0,0,0,0], '-', color = (0,0,0), label='correct') # true estimation
	#ax.plot([0] + s[:5] + [3000], [0]+ s[:5]+[3000], '-', color = (0,0,0), label='correct') # true estimation

	ax.plot(s[:5], bwa[:5], '-', marker='o', color = col_1, label='$BWA$')
	ax.plot(s[:5], gd[:5], '-k', marker='x', label='$GetDistr$')

	for i in range(0,len(a), 5):
		ax.plot(s[i:i+5], bwa[i:i+5], '-', marker='o', color = col_1)
		ax.text(s[i+4], bwa[i+4], str(a[i]), fontsize=10)


		ax.plot(s[i:i+5], gd[i:i+5], '-k', marker='x')
		ax.text(s[i+4], gd[i+4], str(a[i]), fontsize=10)

	ax.set_xlabel('$\sigma$',fontsize=24)
	ax.set_ylabel('$\sigma - \hat{\sigma}$',fontsize=24)

	# Now add the legend with some customizations.
	legend = ax.legend(loc='lower left', shadow=True)
	ax.grid(True )
	plt.savefig("sigma_plot.eps", format='eps')


def plot_mean(a,s, mean_bwa,mean_gd):
	lines_1 = ["--", "_", "-.", ":"]
	col_1 = (0.6, 0.6, 0.6)
	lines_2 = ["--", "_", "-.", ":"]
	col_2 = (0.3,0.3,0.3)
	fig, ax = plt.subplots()

	bwa = mean_bwa
	gd = mean_gd
	#bwa = map(lambda x: x[0]-x[1], zip(mean_bwa,s))
	#gd = map(lambda x: x[0]-x[1], zip(mean_gd,s))

	#ax.plot([0] + s[:5] + [3000], [10000]*7, '-k', color = (0,0,0), label='correct') # true estimation
	#ax.plot([0] + s[:5] + [3000], [0]+ s[:5]+[3000], '-', color = (0,0,0), label='correct') # true estimation
	ax.set_xlim(500,3000)
	ax.plot(s[:5], bwa[:5], '-k',marker='o', color = col_1 , label='$BWA$')
	ax.plot(s[:5], gd[:5], '-k', marker='x', label='$GetDistr$')

	for i in range(0,len(a), 5):
		ax.plot(s[i:i+5], bwa[i:i+5], '-k', marker='o', color = col_1 )
		ax.text(s[i+4], bwa[i+4], str(a[i]),fontsize=10)
		#ax.scatter(s[i:i+5], bwa[i:i+5])

		ax.plot(s[i:i+5], gd[i:i+5],'-k', marker='x' )
		ax.text(s[i+4], gd[i+4], str(a[i]),fontsize=10)
		#ax.scatter(s[i:i+5], gd[i:i+5])

	ax.set_xlabel('$\sigma$',fontsize=24)
	ax.set_ylabel('$\hat{\mu}$',fontsize=24)

	# Now add the legend with some customizations.
	legend = ax.legend(loc='lower left', shadow=True)
	ax.grid(True ) #b=True, which='minor', color='k', linestyle=':')
	plt.savefig("mean_plot.eps", format='eps')

def plot_sigma3d(a,s, sigma_bwa,sigma_gd):
	print sigma_bwa
	# Create plots with pre-defined labels.
	# Alternatively, you can pass labels explicitly when calling `legend`.
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	#bwa = map(lambda x: x[0]-x[1], zip(sigma_bwa,s))
	#gd = map(lambda x: x[0]-x[1], zip(sigma_gd,s))
	c, m = ('r', 'o')
	ax.contour(a, s, sigma_bwa, c=c, marker=m)
	c, m = ('b', '^')
	ax.scatter(a, s, sigma_gd, c=c, marker=m)

	ax.set_xlabel('$a$',fontsize=24)
	ax.set_ylabel('$\sigma$',fontsize=24)
	ax.set_zlabel('Z Label',fontsize=24)


	plt.savefig("sigma_plot.eps", format='eps')

def plot_mean3d(a,s, mean_bwa,mean_gd):
	print sigma_bwa
	# Create plots with pre-defined labels.
	# Alternatively, you can pass labels explicitly when calling `legend`.
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	c, m = ('r', 'o')
	ax.scatter(a, s, mean_bwa, c=c, marker=m)
	c, m = ('b', '^')
	ax.scatter(a, s, mean_gd, c=c, marker=m)

	ax.set_xlabel('$a$',fontsize=24)
	ax.set_ylabel('$\sigma$',fontsize=24)
	ax.set_zlabel('Z Label',fontsize=24)
	plt.savefig("mean_plot.eps", format='eps')

experiment_file = open(sys.argv[1],'r')

##
# vectors indexed like:
# a = [reflen1,reflen2,...]
# s = [sigma1,sigma2,...]

a = []
s = []
sigma_bwa = []
mean_bwa = []
sigma_gd = []
mean_gd = []

##
# Parsing the values
for line in experiment_file:
	# skipping header
	if line[0] == 'l':
		continue

	a.append(int(line.split('\t')[4]))
	s.append(int(line.split('\t')[1]))
	sigma_bwa.append(int(line.split('\t')[8]))
	mean_bwa.append(int(line.split('\t')[7]))
	sigma_gd.append(int(line.split('\t')[10]))
	mean_gd.append(int(line.split('\t')[9]))


plot_sigma(a,s, sigma_bwa,sigma_gd)
plot_mean(a,s, mean_bwa,mean_gd)
