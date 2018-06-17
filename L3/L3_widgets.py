"""
Implements widgets that are used in the L2 notebook. Each widget is 
implemented as a class that can be imported. To use a widget, create
an object of the class in question and call its display method.

Example:
========
from L2_widgets import ReactionWidget
ReactionWidget().display()
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from ipywidgets import interact, IntSlider, FloatSlider


class VoltageClampWidget():
    """A widget solving the simple voltage clamp circuit with a step change."""
    Cm = 0.05 # nF
    Rs = 10 # MOhm

    def V_target(self, t):
        return (t > 2)*(t < 6)*40 - 80

    def I_app(self, V, t):
        return (V - self.V_target(t))/self.Rs

    def I_cap(self, V, t):
        return -self.I_app(V, t)

    def dV_dt(self, V, t):
        return -self.I_app(V, t)/self.Cm

    def solve_and_plot(self, Cm, Rs):
        self.Cm = Cm
        self.Rs = Rs
        time = np.arange(0, 10.1, 0.01)
        V0 = -80
        V = odeint(self.dV_dt, V0, time)
        V = V[:, 0]

        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

        # Potential Plot
        ax1.plot(time, V)
        ax1.plot(time, self.V_target(time))
        ax1.set_ylabel('Mem. Potential [mV]')

        # I_cap plot
        ax2.plot(time, self.I_cap(V, time))
        ax2.set_ylabel('Cap. current [?]')
        ax2.set_xlabel('Time [ms]')
        ax2.axis((0, 10, -5, 5))
        plt.show()
        
    def display(self):
        widget = interact(self.solve_and_plot,
                          Cm = FloatSlider(value=0.05, min=0.005, max=0.1, step=0.005),  
                          Rs = IntSlider(value=10, min=5, max=20, step=1))


class MembraneWidget():
    Cm = 0.05

    g_Na = 0.005
    g_Ca = 0.002
    g_K = 0.02

    E_Na = 70
    E_K = -86
    E_Ca = 114

    def dV_dt(self,V, t):
        g_Na, g_Ca, g_K = self.g_Na, self.g_Ca, self.g_K
        E_Na, E_Ca, E_K = self.E_Na, self.E_Ca, self.E_K
        return -(g_Na*(V-E_Na) + g_K*(V - E_K) + g_Ca*(V - E_Ca))/self.Cm

    def solve_and_plot(self, g_Na, g_Ca, g_K):
        self.g_Na = g_Na
        self.g_Ca = g_Ca
        self.g_K = g_K

        time = np.linspace(0, 200, 2001)
        V0 = 0
        V = odeint(self.dV_dt, V0, time)
        V_eq = V[-1, 0]

        plt.axhline(V_eq, linewidth=2.0)
        plt.axhline(114, alpha=0.5, color='black', linestyle='--')
        plt.axhline(70, alpha=0.5, color='black', linestyle='--')
        plt.axhline(-86, alpha=0.5, color='black', linestyle='--')
        plt.xlabel('Time [ms]')
        plt.ylabel('Mem. Potential [mV]')
        plt.axis((0., 200, -90, 120))
        plt.show()

    def display(self):
        widget = interact(self.solve_and_plot,
                          g_Na = IntSlider(value=5, min=0, max=30, step=1),
                          g_Ca = IntSlider(value=5, min=0, max=30, step=1),
                          g_K  = IntSlider(value=5, min=0, max=30, step=1))


if __name__ == '__main__':
   #VoltageClampWidget().display()
   MembraneWidget().display()