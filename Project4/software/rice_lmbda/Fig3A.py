import rice_lmbda as rice
from scipy.integrate import odeint
import math
import numpy as np
import pylab

t = np.linspace(0,100,101)

SLvals = np.linspace(1.4,2.4,11)
lmbdavals = SLvals/1.85
Cai = np.linspace(0,10,101)
force_index = rice.monitor_indices("active")

for l in lmbdavals:
    init = rice.init_state_values(lmbda=l,lmbda_a=l)
    Fss = np.empty_like(Cai)
    for i in range(len(Cai)):
        p = (rice.init_parameter_values(start_time=1000,Ca_diastolic=Cai[i],SLmin=2.5),)
        s = odeint(rice.rhs,init,t,p)
        m = rice.monitor(s[-1],t[-1],p[0])
        Fss[i] = m[force_index]

    pylab.semilogx(Cai,Fss)


Fmax = 0.17
Ca50 = 3.0
n = 7.6
Fh = Fmax*np.power(Cai,n)/(math.pow(Ca50,n)+np.power(Cai,n))
pylab.plot(Cai,Fh,'--')

Fmax =0.92
Ca50= 0.87
Fh = Fmax*np.power(Cai,n)/(math.pow(Ca50,n)+np.power(Cai,n))
pylab.plot(Cai,Fh,'--')



pylab.show()


#pylab.plot(t,s[:,0],label='force')
#pylab.plot(t,s[:,5],label='N')
#pylab.plot(t,1-s[:,5]-s[:,8]-s[:,9],label='P')
#pylab.plot(t,s[:,8],label='XBpost')
#pylab.plot(t,s[:,9],label='XBpre')
#pylab.legend()
#pylab.show()

