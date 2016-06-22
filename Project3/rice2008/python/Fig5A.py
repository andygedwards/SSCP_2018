import rice_model_2008 as rice
from scipy.integrate import odeint
import numpy as np 
import pylab

t = np.linspace(0,1000,101)
SLvals = np.linspace(1.8,2.3,6)
p = rice.init_parameter_values(SLmin=2.5)
force_ind = rice.monitor_indices("active")     

for s in SLvals:
    init = rice.init_state_values(SL=s)
    s = odeint(rice.rhs,init,t,(p,))
    force = []
    for tn,sn in zip(t,s):
        m = rice.monitor(sn,tn,p)
        force.append(m[force_ind])
    pylab.plot(t,force)

pylab.show()

