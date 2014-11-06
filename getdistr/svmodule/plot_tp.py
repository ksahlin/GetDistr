"""
Demo of the histogram (hist) function with multiple data sets.

Plot histogram with multiple sample sets and demonstrate:

    * Use of legend with multiple sample sets
    * Stacked bars
    * Step curve with a color fill
    * Data sets of different sample sizes
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt

from operator import add

def plot2(x_y_ins_true, x_y_ins_false, x_y_del_true, x_y_del_false):
	#fig, axes = plt.subplots(nrows=2, ncols=3)
	#ax0, ax1, ax2,ax3, ax4, ax5 = axes.flat	


	print 'here',len(x_y_ins_true)
	print len(x_y_del_true)
	#line_up, = plt.plot([1,2,3], label='Line 2')
	#line_down, = plt.plot([3,2,1], label='Line 1')
	fig, ax = plt.subplots()
	fig.suptitle('$H_0^{gd}$, $\mu=300$')
	plt.xlabel('indel size')
	plt.ylabel('# detected')
	sd_del_true, x_del_true, y_del_true, fp_del  = zip(*x_y_del_true)	
	sd_ins_true, x_ins_true, y_ins_true, fp_ins = zip(*x_y_ins_true)
	tot_fps = map(add, fp_del,fp_ins)
	tot_fps = map(lambda x: -x, tot_fps)

	print len(x_del_true[:6]),len( y_del_true[12:]), len( y_ins_true[12:])
	print tot_fps
	sd1, = ax.plot(x_del_true[:6], y_del_true[:6],'^-b',label='$TP_{del}$,$\sigma$=25')
	sd1, = ax.plot(x_ins_true[:6], y_ins_true[:6],'o--b',label='$TP_{ins}$,$\sigma$=25')
	sd2, = ax.plot(x_del_true[6:12], y_del_true[6:12],'^-r',label='$TP_{del}$,$\sigma$=50')
	sd2, = ax.plot(x_ins_true[6:12], y_ins_true[6:12],'o--r',label='$TP_{ins}$,$\sigma$=50')
	sd3, = ax.plot(x_del_true[12:], y_del_true[12:],'^-g',label='$TP_{del}$,$\sigma$=75')
	sd3, = ax.plot(x_ins_true[12:], y_ins_true[12:],'o--g',label='$TP_{ins}$,$\sigma$=75')
	fps20_25,fps30_25,fps40_25,fp50_25,fp75_25,fp100_25, = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[:6],color='b',width=5, label='$FP_{all}$,$\sigma=25$')
	fps20_50,fps30_50,fps40_50,fp50_50,fp75_50,fp100_50, = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[6:12],color='r',width=5, label='$FP_{all}$,$\sigma=50$',bottom=tot_fps[:6])
	fps20_75,fps30_75,fps40_75,fp50_75,fp75_75,fp100_75, = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[12:],color='g',width=5, label='$FP_{all}$,$\sigma=75$',bottom=map(add, tot_fps[:6],tot_fps[6:12]))
	legend = ax.legend(loc='center left', shadow=True)
	#ax.set_yticklabels([str(abs(x-31)) for x in ax.get_yticks()])
	ax.set_xlim([-32, 110])
	ax.set_ylim([-32, 110])
	plt.fill_between(x_del_true[:6], y_del_true[:6], y_ins_true[:6], color='b', alpha='0.3')
	plt.fill_between(x_del_true[:6], y_del_true[6:12], y_ins_true[6:12], color='r', alpha='0.3')
	plt.fill_between(x_del_true[:6], y_del_true[12:], y_ins_true[12:], color='g', alpha='0.3')
	#plt.show()
	plt.savefig('/Users/ksahlin/Dropbox/GetDist/paper-Asmblthn-contest/figures/Clever_true_300.pdf', format='pdf')

	fig, ax = plt.subplots()
	fig.suptitle('$H_0$, $\mu=300$')
	plt.xlabel('indel size')
	plt.ylabel('# detected')	
	sd_del_true, x_del_false, y_del_false, fp_del  = zip(*x_y_del_false)	
	sd_ins_true, x_ins_false, y_ins_false, fp_ins  = zip(*x_y_ins_false)
	tot_fps = map(add, fp_del,fp_ins)
	tot_fps = map(lambda x: -x, tot_fps)
	sd1, = ax.plot(x_del_false[:6], y_del_false[:6],'^-b',label='$TP_{del}$,$\sigma$=25')
	sd1, = ax.plot(x_ins_false[:6], y_ins_false[:6],'o--b',label='$TP_{ins}$,$\sigma$=25')
	sd2, = ax.plot(x_del_false[6:12], y_del_false[6:12],'^-r',label='$TP_{del}$,$\sigma$=50')
	sd2, = ax.plot(x_ins_false[6:12], y_ins_false[6:12],'o--r',label='$TP_{ins}$,$\sigma$=50')
	sd3, = ax.plot(x_del_false[12:], y_del_false[12:],'^-g',label='$TP_{del}$,$\sigma$=75')
	sd3, = ax.plot(x_ins_false[12:], y_ins_false[12:],'o--g',label='$TP_{ins}$,$\sigma$=75')
	fps20_25,fps30_25,fps40_25,fp50_25,fp75_25,fp100_25 = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[:6],color='b',width=5, label='$FP_{all}$,$\sigma=25$')
	fps20_50,fps30_50,fps40_50,fp50_50,fp75_50,fp100_50 = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[6:12],color='r',width=5, label='$FP_{all}$,$\sigma=50$',bottom=tot_fps[:6])
	fps20_75,fps30_75,fps40_75,fp50_75,fp75_75,fp100_75 = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[12:], color='g',width=5, label='$FP_{all}$,$\sigma=75$',bottom=map(add, tot_fps[:6],tot_fps[6:12]))
	legend = ax.legend(loc='center left', shadow=True)
	ax.set_xlim([-32, 110])
	ax.set_ylim([-32, 110])
	plt.fill_between(x_del_false[:6], y_del_false[:6], y_ins_false[:6], color='b', alpha='0.3')
	plt.fill_between(x_del_false[:6], y_del_false[6:12], y_ins_false[6:12], color='r', alpha='0.3')
	plt.fill_between(x_del_false[:6], y_del_false[12:], y_ins_false[12:], color='g', alpha='0.3')

	#plt.show()
	plt.savefig('/Users/ksahlin/Dropbox/GetDist/paper-Asmblthn-contest/figures/Clever_false_300.pdf', format='pdf')


def plot(ins_false,ins_true,del_false,del_true):
	n_bins = 20
	#x = np.random.randn(10, 2)
	#print x
	# inse = [20]*30 + [30]*55 + [40]*99 
	# dele =  [70]*30 + [50]*55 + [24]*99
	fig, axes = plt.subplots(nrows=2, ncols=3)
	ax0, ax1, ax2,ax3, ax4, ax5 = axes.flat

	i = 0
	for sd in (25,50,75):
		inse = ins_false[sd]
		dele = del_false[sd]
		length = max(len(inse),len(dele))
		least_points = min(len(inse),len(dele))
		if len(inse) == least_points:
			for k in range(length-least_points):
				inse.append(1000)
		elif len(dele) == least_points:
			for k in range(length-least_points):
				dele.append(1000)

		print len(inse), len(dele)
		data = zip(inse,dele)
		data = np.asarray(data)
		print data
	# for p1,p2 in data:
	# 	if p1 <=0:
	# 		print'heyyy'
	# 	if p2 <= 0:
	# 		print 'loooool'

	#print data


		colors = ['black', 'grey']
		#axes.flat[i].set_ylim(0.1,100)
		axes.flat[i].hist(data, n_bins, histtype='bar', color=colors, label=colors, range=(0,100),  bottom=0.01)
		axes.flat[i].set_title('$H_0,\sigma$:{0}'.format(sd))
		axes.flat[i].set_yscale('log')
		#axes.flat[i].axis(0,100,0,100)
		i += 1


	for sd in (25,50,75):
		inse = ins_true[sd]
		dele = del_true[sd]
		length = max(len(inse),len(dele))
		least_points = min(len(inse),len(dele))
		if len(inse) == least_points:
			for k in range(length-least_points):
				inse.append(1000)
		elif len(dele) == least_points:
			for k in range(length-least_points):
				dele.append(1000)

		#print len(inse), len(dele)
		data = zip(inse,dele)
		data = np.asarray(data)

	#print data


		colors = ['black', 'grey']
		#axes.flat[i].set_ylim(0.1,100)
		axes.flat[i].hist(data, n_bins, histtype='bar', color=colors, label=colors, range=(0,100),  bottom=0.01)
		axes.flat[i].set_title('$H_0^g,\sigma$:{0}'.format(sd))
		axes.flat[i].set_yscale('log')
		#axes.flat[i].axis(0,100,0,100)
		i += 1

		#ax0.hist(data, n_bins, histtype='bar', color=colors, label=colors)
	#ax0.legend(prop={'size': 10})
	#ax0.set_title('bars with legend')
	plt.show()
	#plt.savefig('/Users/ksahlin/Dropbox/GetDist/paper-v3/GD_results_for_ms/indel_sd3_500/insertion/figure', format='eps')

def parse_results(file_):
	false_hyp_tp = {25:[],50:[],75:[]}
	true_hyp_tp = {25:[],50:[],75:[]}
	x_y_false = []
	x_y_true = []

	for line in file_:
		cols = line.split('\t')
		#print cols
		if cols[0] == 'method' or len(cols) != 7:
			continue
		gap = int(cols[2])
		n = int(cols[3])
		sd = int(cols[1])
		fp = int(cols[4])
		if cols[0] == 'Clever':
			x_y_false.append((sd,gap,n,fp))
			false_hyp_tp[sd] += [gap]*n
		if cols[0] == 'Clever true H0':
			x_y_true.append((sd,gap,n,fp))
			true_hyp_tp[sd] += [gap]*n
	return false_hyp_tp, true_hyp_tp, x_y_true, x_y_false

def main(args):
	ins_false, ins_true, x_y_ins_true, x_y_ins_false = parse_results(open(args.inse,'r'))
	del_false, del_true, x_y_del_true, x_y_del_false = parse_results(open(args.dele,'r'))
	#print ins_false,del_false

	#plot(ins_false,ins_true,del_false,del_true)
	plot2(x_y_ins_true, x_y_ins_false, x_y_del_true, x_y_del_false)



if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
    parser.add_argument('inse', type=str, help='Path to insertion file. ')
    parser.add_argument('dele', type=str, help='Path to deletion file. ')

    args = parser.parse_args()

    main(args)