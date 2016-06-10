# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 16:43:34 2016

@author: Andy
"""

import numpy as np
import matplotlib.pylab as plt
import ObFunc
import scipy.optimize as opt

SSA_data = np.loadtxt("SS",dtype='float')

V = SSA_data[:,0] # voltage
I = SSA_data[:,1] # current

# The parameter vector
init_params = [48.4, 49.2, 48.4, 49.2, 48.4, 
               13.6, 1, 0.023, 48.4, 13.6, 
               0.00001, 0.00001, 0.3]
step_length = 1000

dats = Activation(init_params,V,I,step_length)
model_I = dats['I_peak']
Po = dats['Po']
t = dats['t']

plt.figure()
plt.plot(t,Po)
plt.xlabel('time [ms]')
plt.ylabel('Open probability (normalized conductance)')

plt.figure()
plt.plot(t,Po)
plt.xlabel('time [ms]')
plt.ylabel('Open probability (normalized conductance)')
plt.xlim([0, 50])
plt.ylim([0, 1.1])

V = SSA_data[:,0]
I = SSA_data[:,1]
plt.figure()
plt.plot(V,I,'b-')
plt.plot(V,model_I,'r-')
plt.xlabel('Step potential [mV]')
plt.ylabel('Current (A/F)')

# Now run the optimization

[P_opt, f_opt, iters, funcalls, warnflag] = opt.fmin(cost, init_params, args=(V,I,step_length))

#outs = IKur.markov(V,params)

