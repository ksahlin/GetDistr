'''
Created on Aug 18, 2013

@author: ksahlin
'''

from pylab import *


o = [520]
import model
est = model.NormalModel(500,100,100,s=100)
print est.expected_mean(0,10000)

likelihood_fcn  = est.get_likelihood_function(o,10000,20,b=10000)
x,y = zip(*likelihood_fcn)
avg_obs = sum(o)/len(o)
x = [i + avg_obs for i in x]
xlabel('X (mean fragment length)')
ylabel('ML function value')
plot(x,y)
show()
est.infer_mean(o,10000,1,b=10000)

#bug!?
# o = [520] -> X = 530
# o = [330] -> X = 620
# o = [310] -> X = 650
# o = [290] -> X = 690
# o = [270] -> X = 750
# o = [250] -> X = 870
# o = [230] -> X = 1130
# o = [220] -> X = 1210
# o = [210] -> X = 1200