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
			print line
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

def plot(reapr_database, getdistr_database, outfile):

	errorsizes = reapr_database.get_error_sizes()
	gapsizes = reapr_database.get_gap_sizes()

	for gap in gapsizes:
		fig, ax = plt.subplots()


		for z in  errorsizes:
			N= len(errorsizes)
			print errorsizes
			ind = np.arange(N*2,step=2)  # the x locations for the groups
			print ind
			width = 0.25       # the width of the bars
			TP_r,FP_r = reapr_database.get_TP_FP(gap=gap)
			TP_g,FP_g = getdistr_database.get_TP_FP(gap=gap)
			print TP_r

			#gaps = np.arange(0,3*N,3)  # the x locations for the groups
			t_r = ax.bar(ind, TP_r, width, color='g', alpha=1)
			t_p = ax.bar(ind+width*0.8, FP_r, width, color='r', alpha=1)
			g_r = ax.bar(ind+2*width*0.8, TP_g, width, color='b', alpha=1)
			g_p = ax.bar(ind+3*width*0.8, FP_g, width, color='y', alpha=1)

		# add some text for labels, title and axes ticks
		ax.set_ylabel('Number of TP/FP')
		ax.set_title('Assembly errors gap={0}'.format(gap))
		ax.set_xticks(ind+2*width)
		ax.set_xticklabels( errorsizes )
		ax.legend( (t_r[0], t_p[0], g_r[0], g_p[0]), ('reaper TP', 'reapr FP', 'GetDistr TP', 'Getdistr FP') )
		plt.savefig(open('{0}_{1}.png'.format(outfile,gap),'w'))


def main(args):
	reapr_database = DataBase()
	reapr_database.read_in_data(open(args.result_file1,'r'))
	getdistr_database = DataBase()
	getdistr_database.read_in_data(open(args.result_file2,'r'))

	plot(reapr_database, getdistr_database, args.outfile)


if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Plot results of correctly estimated variants.")
	parser.add_argument('result_file1', type=str, help='Path to result file for reapr. ')
	parser.add_argument('result_file2', type=str, help='Path to result file for Getdistr. ')
	parser.add_argument('outfile', type=str, help='Path to plot file. ')

	# parser.add_argument('mean', type=int, help='mean insert size. ')
	# parser.add_argument('stddev', type=int, help='Standard deviation of insert size ')


	args = parser.parse_args()
	main(args)

