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
from operator import add

try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 22})
    import seaborn as sns
    sns.set_palette("husl", desat=.6)

except:
    print("Counld not import either matplotlib or seaborn!")



def plot(args):
    #NOT STACKED! We need to add up the TP/FP in each case
    # install seaborn
    # manually format the input before parsing results!
    # create a black white plot with dotted or marked bars
    # plot only sigma 50 and 75 for mu 300,400,500, remove indel size 20
    indelsizes = 6
    TP_MU300_I = (80, 70, 50, 34, 25, 25)
    TP_MU300_D = (120, 105, 62, 50, 53, 55)
    FP_MU300 = (-50, -40, -60, -15, -15, -25)

    ind = np.arange(indelsizes) + .15 # the x locations for the groups
    width = 0.2       # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, TP_MU300_D, width, color='grey',label='INS')
    rects2 = ax.bar(ind, TP_MU300_I, width, color='k',label='DEL') 
    rects3 = ax.bar(ind, FP_MU300, width, color='white', label='FP')

    TP25_MU400 = (140, 90, 78, 65, 50, 55)
    TP50_MU400 = (130, 80, 70, 60, 45, 55)
    TP75_MU400 = (120, 60, 60, 55, 44, 55)

    xtra_space = 0.05
    rects2 = ax.bar(ind + width + xtra_space , TP25_MU400, width, color='orange')
    rects2 = ax.bar(ind + width + xtra_space, TP50_MU400, width, color='cyan')
    rects2 = ax.bar(ind + width + xtra_space, TP75_MU400, width, color='purple')

    TP25_MU500 = (140, 90, 78, 65, 50, 55)
    TP50_MU500 = (130, 80, 70, 60, 45, 55)
    TP75_MU500 = (120, 60, 60, 55, 44, 55)

    rects2 = ax.bar(ind + 2*width + 2*xtra_space , TP25_MU500, width, color='r') 
    rects2 = ax.bar(ind + 2*width + 2*xtra_space, TP50_MU500, width, color='cyan') 
    rects2 = ax.bar(ind + 2*width + 2*xtra_space, TP75_MU500, width, color='purple')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Population, millions')
    ax.set_title('Population: Age Structure')

    ax.set_xticks(ind+width+xtra_space)
    ax.set_xticklabels( ('20', '30', '40', '50', '75', '100') )

    legend = ax.legend(loc='upper right', shadow=True)

    plt.savefig(args.out + '.eps', format='eps')
    print('ENDING')

# def plot(x_y_ins_true, x_y_ins_false, x_y_del_true, x_y_del_false):
#     #fig, axes = plt.subplots(nrows=2, ncols=3)
#     #ax0, ax1, ax2,ax3, ax4, ax5 = axes.flat    


#     print 'here',len(x_y_ins_true)
#     print len(x_y_del_true)
#     #line_up, = plt.plot([1,2,3], label='Line 2')
#     #line_down, = plt.plot([3,2,1], label='Line 1')
#     fig, ax = plt.subplots()
#     fig.suptitle("$H'_0$, $\mu=300$")
#     plt.xlabel('indel size')
#     plt.ylabel('# detected')
#     sd_del_true, x_del_true, y_del_true, fp_del  = zip(*x_y_del_true)   
#     sd_ins_true, x_ins_true, y_ins_true, fp_ins = zip(*x_y_ins_true)
#     tot_fps = map(add, fp_del,fp_ins)
#     tot_fps = map(lambda x: -x, tot_fps)

#     print len(x_del_true[:6]),len( y_del_true[12:]), len( y_ins_true[12:])
#     print tot_fps
#     sd1, = ax.plot(x_del_true[:6], y_del_true[:6],'^-b',label='$TP_{del}$,$\sigma$=25')
#     sd1, = ax.plot(x_ins_true[:6], y_ins_true[:6],'o--b',label='$TP_{ins}$,$\sigma$=25')
#     sd2, = ax.plot(x_del_true[6:12], y_del_true[6:12],'^-r',label='$TP_{del}$,$\sigma$=50')
#     sd2, = ax.plot(x_ins_true[6:12], y_ins_true[6:12],'o--r',label='$TP_{ins}$,$\sigma$=50')
#     sd3, = ax.plot(x_del_true[12:], y_del_true[12:],'^-g',label='$TP_{del}$,$\sigma$=75')
#     sd3, = ax.plot(x_ins_true[12:], y_ins_true[12:],'o--g',label='$TP_{ins}$,$\sigma$=75')
#     fps20_25,fps30_25,fps40_25,fp50_25,fp75_25,fp100_25, = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[:6],color='b',width=5, label='$FP_{all}$,$\sigma=25$')
#     fps20_50,fps30_50,fps40_50,fp50_50,fp75_50,fp100_50, = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[6:12],color='r',width=5, label='$FP_{all}$,$\sigma=50$',bottom=tot_fps[:6])
#     fps20_75,fps30_75,fps40_75,fp50_75,fp75_75,fp100_75, = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[12:],color='g',width=5, label='$FP_{all}$,$\sigma=75$',bottom=map(add, tot_fps[:6],tot_fps[6:12]))
#     legend = ax.legend(loc='center left', shadow=True)
#     #ax.set_yticklabels([str(abs(x-31)) for x in ax.get_yticks()])
#     ax.set_xlim([-32, 110])
#     ax.set_ylim([-32, 110])
#     plt.fill_between(x_del_true[:6], y_del_true[:6], y_ins_true[:6], color='b', alpha='0.3')
#     plt.fill_between(x_del_true[:6], y_del_true[6:12], y_ins_true[6:12], color='r', alpha='0.3')
#     plt.fill_between(x_del_true[:6], y_del_true[12:], y_ins_true[12:], color='g', alpha='0.3')
#     #plt.show()
#     ax.set_rasterized(True)
#     plt.savefig('/Users/ksahlin/Documents/workspace/getdistr_paper/data/detection_clever/plots/Clever_true_300.eps', format='eps')

#     fig, ax = plt.subplots()
#     fig.suptitle('$H_0$, $\mu=300$')
#     plt.xlabel('indel size')
#     plt.ylabel('# detected')    
#     sd_del_true, x_del_false, y_del_false, fp_del  = zip(*x_y_del_false)
#     sd_ins_true, x_ins_false, y_ins_false, fp_ins  = zip(*x_y_ins_false)
#     tot_fps = map(add, fp_del,fp_ins)
#     tot_fps = map(lambda x: -x, tot_fps)
#     sd1, = ax.plot(x_del_false[:6], y_del_false[:6],'^-b',label='$TP_{del}$,$\sigma$=25')
#     sd1, = ax.plot(x_ins_false[:6], y_ins_false[:6],'o--b',label='$TP_{ins}$,$\sigma$=25')
#     sd2, = ax.plot(x_del_false[6:12], y_del_false[6:12],'^-r',label='$TP_{del}$,$\sigma$=50')
#     sd2, = ax.plot(x_ins_false[6:12], y_ins_false[6:12],'o--r',label='$TP_{ins}$,$\sigma$=50')
#     sd3, = ax.plot(x_del_false[12:], y_del_false[12:],'^-g',label='$TP_{del}$,$\sigma$=75')
#     sd3, = ax.plot(x_ins_false[12:], y_ins_false[12:],'o--g',label='$TP_{ins}$,$\sigma$=75')
#     fps20_25,fps30_25,fps40_25,fp50_25,fp75_25,fp100_25 = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[:6],color='b',width=5, label='$FP_{all}$,$\sigma=25$')
#     fps20_50,fps30_50,fps40_50,fp50_50,fp75_50,fp100_50 = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[6:12],color='r',width=5, label='$FP_{all}$,$\sigma=50$',bottom=tot_fps[:6])
#     fps20_75,fps30_75,fps40_75,fp50_75,fp75_75,fp100_75 = ax.bar([17.5,27.5,37.5,47.5,72.5,95],tot_fps[12:], color='g',width=5, label='$FP_{all}$,$\sigma=75$',bottom=map(add, tot_fps[:6],tot_fps[6:12]))
#     legend = ax.legend(loc='center left', shadow=True)
#     ax.set_xlim([-32, 110])
#     ax.set_ylim([-32, 110])
#     plt.fill_between(x_del_false[:6], y_del_false[:6], y_ins_false[:6], color='b', alpha='0.3')
#     plt.fill_between(x_del_false[:6], y_del_false[6:12], y_ins_false[6:12], color='r', alpha='0.3')
#     plt.fill_between(x_del_false[:6], y_del_false[12:], y_ins_false[12:], color='g', alpha='0.3')

#     #plt.show()
#     ax.set_rasterized(True)
#     plt.savefig('/Users/ksahlin/Documents/workspace/getdistr_paper/data/detection_clever/plots/Clever_false_300.eps', format='eps')
#     print 'ENDING'


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
    # ins_false, ins_true, x_y_ins_true, x_y_ins_false = parse_results(open(args.inse,'r'))
    # del_false, del_true, x_y_del_true, x_y_del_false = parse_results(open(args.dele,'r'))
    #print ins_false,del_false

    plot(args)
    # plot(x_y_ins_true, x_y_ins_false, x_y_del_true, x_y_del_false)



if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
    parser.add_argument('root', type=str, help='Path to root folder. ')
    parser.add_argument('out', type=str, help='Prefix to outfile plot. ')

    args = parser.parse_args()

    main(args)