import rice_model_2008 as rice
from scipy.integrate import odeint
import numpy as np
import pylab


t = np.linspace(0,650,101)
init = rice.init_state_values()
SL_ind = rice.state_indices("SL")
Cai = np.linspace(0,10,101)
sl_ss = []
for ca in Cai:
    p = (rice.init_parameter_values(start_time=1000,Ca_diastolic=ca,SEon=0.0),)
    s = odeint(rice.rhs,init,t,p)
    sl_ss.append(s[-1,SL_ind])
pylab.ylim([1.4,2.0])
pylab.ylabel("SL (micrometers)")
pylab.xlabel("[Ca] (microM)")
pylab.semilogx(Cai,sl_ss)
pylab.show()

"""pylab.ylim([1.4,2.0])
pylab.ylabel("SL (micrometers)")
pylab.xlim([550,1000])
pylab.xlabel("time (ms)")
pylab.legend()

#compute fitted Hill velocity:
f_hill = linspace(0,0.8,101)
a_hill = 0.16
b_hill = 2.4
v_max = 14
v_hill = (v_max*a_hill-b_hill*f_hill)/(f_hill+a_hill)

pylab.figure(2)
pylab.plot(al,v,label="Model")
pylab.plot(f_hill,v_hill,label="Hill curve")
pylab.ylabel("v (micrometers/s")
pylab.xlabel("normalized force")

pylab.legend()
"""
#pylab.ylim([1.4,2.3])
#pylab.xlim([550,1000])
pylab.show()


#pylab.plot(t,s[:,0],label='force')
#pylab.plot(t,s[:,5],label='N')
#pylab.plot(t,1-s[:,5]-s[:,8]-s[:,9],label='P')
#pylab.plot(t,s[:,8],label='XBpost')
#pylab.plot(t,s[:,9],label='XBpre')
#pylab.legend()
#pylab.show()

