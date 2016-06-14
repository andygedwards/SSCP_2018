# ipython --matplotlib
# run enzymes.py
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def enzymes(y, t, k1p, k1m, k2, Js, Jp):
    
    S, E, C, P = y;

    dS = k1m*C - k1p*S*E + Js;
    dE = (k1m+k2)*C - k1p*S*E;
    dC = k1p*S*E - (k2+k1m)*C;
    dP = k2*C - Jp;

    return [dS, dE, dC, dP];


#D = (k1m, k1p, k2, Js, Jp) = (10,20,30,0,0)
D = (k1m, k1p, k2, Js, Jp) = (10,20,30,1,1)

   
Y0 = [1,1,0,0];
T = np.linspace(0,1,1001);
Y = odeint(enzymes, Y0, T, D)

plt.hold(False)
plt.plot(T,Y);
plt.legend(('S','E','C','P'),0)
plt.show()
#plt.plot(T,Y[:,0]+Y[:,2]+Y[:,3]); plt.title('Sum of mass')
#plt.plot(T,Y[:,1]+Y[:,2]); plt.title('Sum of enzymes')


