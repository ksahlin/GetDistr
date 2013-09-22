'''
Created on Sep 18, 2013

@author: ksahlin
'''

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

#import matplotlib
import argparse

import model



#frag_obj =model.NormalModel(500, 20, 100, 50)
#frag_obj.expected_mean(0,10000,10000)
def plot(args):
    """
    Demo of the legend function with a few features.
    
    In addition to the basic legend, this demo shows a few optional features:
    
        * Custom legend placement.
        * A keyword argument to a drop-shadow.
        * Setting the background color.
        * Setting the font size.
        * Setting the line width.
    """


    # Example data

    std_dev = []
    for sigma in args.sigma:
        fragm_dist_object = model.NormalModel(args.mean, sigma, args.readlen, args.softclipped)
        gap_data_points = []
        for gap in args.gaps:
            print gap, args.reflen
            print fragm_dist_object.expected_mean(gap, args.reflen, args.reflen)
            gap_data_points.append(fragm_dist_object.expected_mean(gap, args.reflen, args.reflen))
        std_dev.append(gap_data_points)


    a = np.arange(0, 3, .02)
    #b = np.arange(0, 3, .02)
    c = np.exp(a)
    d = c[::-1]

    # Create plots with pre-defined labels.
    # Alternatively, you can pass labels explicitly when calling `legend`.
    fig, ax = plt.subplots()
    for sigma in args.sigma:
        #TODO: Need a way to pick 'k--' (line style automatically but different)
        ax.plot(a, c, 'k--', label='std dev = {0}'.format(sigma))
    ax.plot(a, c, 'k--', label='std dev = 50')
    ax.plot(a, d, 'k:', label='std dev = 100')
    ax.plot(a, c + d, 'k', label='std dev = 150')

    # Now add the legend with some customizations.
    legend = ax.legend(loc='upper left', shadow=True)

    # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
    frame = legend.get_frame()
    frame.set_facecolor('0.90')

    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize('large')

    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width
    plt.show()

def plot3d(args):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.arange(-5, 5, 0.25)
    Y = np.arange(-5, 5, 0.25)
    X, Y = np.meshgrid(X, Y)
    R = np.sqrt(X ** 2 + Y ** 2)
    Z = np.sin(R)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False)
    ax.set_zlim(-1.01, 1.01)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()



if __name__ == '__main__':

    # parameters
    parser = argparse.ArgumentParser(description="Generate plots of getdistr")
    parser.add_argument("-m", dest="mean", type=int, required=True,
                  help="Mean insert size of library distribution.")
    parser.add_argument("-q", dest="sigma", type=int, nargs='+', required=True,
                  help="std dev insert size of library distribution.")
    parser.add_argument("-r", dest="readlen", type=int, required=True,
                  help="Read sequence length (not insert size!).")
    parser.add_argument("-s", dest="softclipped", type=int, required=True,
                  help="Number of soft clipped bases.")
    parser.add_argument("-a", dest="reflen", type=int, required=True,
                  help="Reference sequence length.")


    # for plotting
    parser.add_argument("-z", dest="gaps", type=int, nargs='+', required=True,
                  help="sample gap sizes to fit the line to.")

    parser.add_argument("-o", dest="outfile", type=str, required=True,
                  help="outfile destination.")



    args = parser.parse_args()
    if len(args.sigma) > 1:
        plot3d(args)
    else:
        plot(args)
