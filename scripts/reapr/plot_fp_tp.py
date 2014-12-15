import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
# import matplotlib.pyplot as plt
# import numpy as np

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# for c, z in zip(['r', 'g', 'b', 'y'], [30, 20, 10, 0]):
#     xs = np.arange(20)
#     ys = np.random.rand(20)

#     # You can provide either a single color or an array. To demonstrate this,
#     # the first bar of each set will be colored cyan.
#     cs = [c] * len(xs)
#     cs[0] = 'c'
#     ax.bar(xs, ys, zs=z, zdir='y', color=cs, alpha=0.8)

# ax.set_xlabel('Gap size')
# ax.set_ylabel('Error')
# ax.set_zlabel('Frequency')

# plt.show()

class DataBase(object):
	"""docstring for DataBase"""
	def __init__(self):
		super(DataBase, self).__init__()
		self.data = {}
	def read_in_data(self,result_file):
		for line in result_file:
			gap, error, TP, FP = line.strip().split()
			gap, error, TP, FP = map(lambda x: int(x),[gap, error, TP, FP])
			if gap not in self.data:
				self.data[gap]= {}
				self.data[gap][error] = {'TP':TP,'FP':FP}
			elif error not in self.data[gap]:
				self.data[gap][error] = {'TP':TP,'FP':FP}			

	def get_TP_FP(self, gap=None,error=None):
		if gap != None and error != None:
			return  (self.data[gap][error]['TP'], self.data[gap][error]['FP'])

		elif gap != None:
			TPs = map(lambda error: self.data[gap][error]['TP'], sorted(self.data[gap]))
			FPs = map(lambda error: self.data[gap][error]['FP'], sorted(self.data[gap]))
			return TPs, FPs

		elif error != None:
			TPs = []
			FPs = []
			for gap in sorted(self.data):
				TPs.append(self.data[gap][error]['TP'])
				FPs.append(self.data[gap][error]['FP'])
			return TPs,FPs

	def get_error_sizes(self):
		gap = self.data.keys()[0]
		return sorted(map(lambda key: key, self.data[gap]))

	def get_gap_sizes(self):
		return sorted(map(lambda key: key, self.data))

def plot(database,outfile):

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	errorsizes = database.get_error_sizes()
	gapsizes = database.get_gap_sizes()

	for z in  errorsizes:
	#for z in  (-1500, -1250, -1000, -750, -500, -250, 0, 250, 500, 750, 1000,1250,1500):
		

		#N = 7
		#TP = (20, 35, 30, 35, 27, 33, 66 )
		#FP = (25, 32, 34, 20, 25, 11, 67)

		N= len(gapsizes)
		TP,FP = database.get_TP_FP(error=z)

		gaps = np.arange(0,3*N,3)  # the x locations for the groups
		width = 0.5       # the width of the bars


		ax.bar(gaps,TP, color='r',zdir='y',zs=z,alpha=0.8) #, yerr=menStd)



		ax.bar(gaps+width, FP, color='y',zdir='y',zs=z, alpha=0.8) #, yerr=womenStd)

		# add some
	ax.set_xlabel('Gap size')
	ax.set_ylabel('Error')
	ax.set_zlabel('Number of TP/FP')
	ax.set_title('Assembly errors')



	# ax.set_xticks(gaps+width)
	# ax.set_xticklabels( ('0', '250', '500', '750', '1000','1250','1500') )
	# ax.set_yticks((-1500, -1250, -1000, -750, -500, -250, 0, 250, 500, 750, 1000,1250,1500))
	# ax.set_yticklabels( (-1500, -1250, -1000, -750, -500, -250, 0, 250, 500, 750, 1000,1250,1500) )

	ax.set_xticks(gaps+width)
	ax.set_xticklabels( gapsizes )
	ax.set_yticks(errorsizes)
	ax.set_yticklabels( errorsizes )

	plt.savefig(open(outfile+'.png','w'))
	#plt.show()


def main(args):
	database = DataBase()
	database.read_in_data(open(args.result_file,'r'))
	#data = parse_result_file()
	plot(database,args.outfile)


if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Plot results of correctly estimated variants.")
	parser.add_argument('result_file', type=str, help='Path to result file. ')
	parser.add_argument('outfile', type=str, help='Path to plot file. ')

	# parser.add_argument('mean', type=int, help='mean insert size. ')
	# parser.add_argument('stddev', type=int, help='Standard deviation of insert size ')


	args = parser.parse_args()
	main(args)

