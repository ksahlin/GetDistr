'''
Created on Sep 18, 2013

@author: ksahlin
'''

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, FixedLocator
import matplotlib.pyplot as plt
import numpy as np

#import matplotlib
import argparse

import model


#TODO: Split up fcns: data_readin() to get vectors plot() and plot3d() for plotting.
#frag_obj =model.NormalModel(500, 20, 100, 50)
#frag_obj.expected_mean(0,10000,10000)



def plot_mean_frag(args):
    """
    Demo of the legend function with a few features.
    
    In addition to the basic legend, this demo shows a few optional features:
    
        * Custom legend placement.
        * A keyword argument to a drop-shadow.
        * Setting the background color.
        * Setting the font size.
        * Setting the line width.
    """


    # Generating data

    std_dev = []
    for sigma in args.sigma:
        fragm_dist_object = model.NormalModel(args.mean, sigma, args.readlen, args.softclipped)
        gap_data_points = []
        for gap in args.gaps:
            #print gap, args.reflen, sigma
            #print fragm_dist_object.expected_mean(gap, args.reflen, args.reflen)
            gap_data_points.append(fragm_dist_object.expected_mean(gap, args.reflen, args.reflen))
        std_dev.append(gap_data_points)

    # Create plots with pre-defined labels.
    # Alternatively, you can pass labels explicitly when calling `legend`.
    fig, ax = plt.subplots()
    lines_ = ["k-", "k--", "k-.", "k:"]
    for i, sigma in enumerate(args.sigma):
        line = lines_[i]
        ax.plot(args.gaps, std_dev[i], line, label='$\sigma = {0}$'.format(sigma))


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

    ax.set_xlabel('$z$', fontsize=24)
    ax.set_ylabel('$x$', fontsize=24)

    plt.savefig(args.outfile)

def plot_stddev_frag(args):


    # Example data

    std_dev = []
    for sigma in args.sigma:
        fragm_dist_object = model.NormalModel(args.mean, sigma, args.readlen, args.softclipped)
        gap_data_points = []
        for gap in args.gaps:
            print gap, args.reflen, sigma
            print fragm_dist_object.expected_standard_deviation(gap, args.reflen, args.reflen)
            gap_data_points.append(fragm_dist_object.expected_standard_deviation(gap, args.reflen, args.reflen))
        std_dev.append(gap_data_points)

    # Create plots with pre-defined labels.
    # Alternatively, you can pass labels explicitly when calling `legend`.
    fig, ax = plt.subplots()
    lines_ = ["k-", "k--", "k-.", "k:"]
    for i, sigma in enumerate(args.sigma):
        line = lines_[i]
        ax.plot(args.gaps, std_dev[i], line, label='$\sigma = {0}$'.format(sigma))


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

    ax.set_xlabel('$z$', fontsize=24)
    ax.set_ylabel('$\sigma_{x|z}$', fontsize=24)

    plt.savefig(args.outfile)


def plot3d_mean_frag(args):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
#    X = np.arange(-5, 5, 0.25)
#    Y = np.arange(-5, 5, 0.25)
#    X, Y = np.meshgrid(X, Y)


    X = np.array(args.gaps)
    Y = np.array(args.sigma)
    X, Y = np.meshgrid(X, Y)
    std_dev = []
    for sigma in args.sigma:
        fragm_dist_object = model.NormalModel(args.mean, sigma, args.readlen, args.softclipped)
        given_sigma = []
        for gap in args.gaps:
            given_sigma.append(fragm_dist_object.expected_mean(gap, args.reflen, args.reflen))
        std_dev.append(given_sigma)
    Z = np.array(std_dev)

#    R = np.sqrt(X ** 2 + Y ** 2)
#    Z = np.sin(R)
#
#    for i in range(len(X)):
#        Z[i] = (X[i] + Y[i])
#    print X, len(X)
#    print Y, len(Y)
#    print Z, len(Z)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.gray, linewidth=0, antialiased=False)
    #cmap=cm.coolwarm
    ax.set_zlim(np.amin(Z) - 0.01, np.amax(Z) + 0.01)
    z_range = range(args.mean, int(np.amax(Z)) + 50, 50)
    ax.zaxis.set_major_locator(FixedLocator(z_range))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.w_xaxis.set_major_locator(FixedLocator(args.gaps))
    ax.w_yaxis.set_major_locator(FixedLocator(args.sigma))
    #ax.xticks(args.gaps)
    #ax.yticks(args.sigma)

    fig.colorbar(surf, shrink=0.5, aspect=5)


    ax.set_xlabel('$z$', fontsize=24)
    ax.set_ylabel('$\sigma$', fontsize=24)
    ax.set_zlabel('$x$', fontsize=24)

    plt.savefig(args.outfile)

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

    parser.add_argument("--3D", dest='threedim', action="store_true",
                  help="Plots a 3D plot with gaps,std_dev and function value as X,Y,Z")

    parser.add_argument("--std", dest='std', action="store_true",
                  help="Plots the expected standard deviation of fragment sizes over a gap.\
                   If --std is not specified, the mean fragment length is plotted (default = mean length)")

    args = parser.parse_args()
    if args.threedim:
        plot3d_mean_frag(args)
    elif args.std:
        plot_stddev_frag(args)
    else:
        plot_mean_frag(args)
