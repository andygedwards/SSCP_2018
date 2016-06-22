import numpy as np
from norm_force import *
import matplotlib.pyplot as plt

t = np.linspace(0,1000,1000)
f = np.zeros_like(t)
for i, tp in enumerate(t):
    f[i] = force_transient(tp)

plt.plot(t,f)
plt.show()


