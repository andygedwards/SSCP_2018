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


class ReactionWidget():
    """A widget solving the Law of Mass action for a combustion reaction."""
    def rhs(self, y, t, k):
        H2, O2, H2O = y
        dH2_dt = -2*k*H2*H2*O2
        dO2_dt = -k*H2*H2*O2
        dH2O_dt = k*H2*H2*O2       
        return [dH2_dt, dO2_dt, dH2O_dt]

    def solve_and_plot(self, k, H2_0, O2_0):
        time = np.arange(0, 5.1, 0.01)
        initial_condition = (H2_0, O2_0, 0.0)
        concentrations = odeint(self.rhs, initial_condition, time, (k,))
        H2, O2, H2O = np.hsplit(concentrations, 3)

        plt.plot(time, H2, label=r'H$_2$')
        plt.plot(time, O2, label=r'O$_2$')
        plt.plot(time, H2O, label=r'H$_2$O')
        plt.axis((0, 5, 0, 5))
        plt.xlabel('Time')
        plt.ylabel('Concentrations')
        plt.legend()
        plt.show()       

    def display(self):
        widget = interact(self.solve_and_plot,
                          k = FloatSlider(value=1, min=0, max=5, step=0.1),
                          H2_0 = FloatSlider(value=3, min=0, max=5, step=0.1),
                          O2_0 = FloatSlider(value=3, min=0, max=5, step=0.1))


if __name__ == '__main__':
    ReactionWidget().display()
