import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.integrate

from array import array
from pylab import *
from scipy.integrate.odepack import odeint

from scipy.integrate import odeint
from math import exp, log, sqrt, pi
from ipywidgets import interact, IntSlider, FloatSlider


class CICR_Widget ():

  def solve_and_plot(self, units, tau_SR, tau_cyt, BC):
    #initial conditions
    #units = 10
    self.units = units
    self.tau_SR = tau_SR
    self.tau_cyt = tau_cyt
    self.BC = BC

    y0 = [0.0795,4.1725];
    y0MultiUnit = np.tile(y0,self.units)
    time = np.linspace(0.0,1000.0, num = 1000)

    #ode function call
    y = scipy.integrate.odeint(self.CICR, y0MultiUnit, time);

    plt.subplot(2,1,1); plt.plot(time,y[:,0::2]); plt.title('c'); plt.axis([0, 1000, 0, 1])
    plt.subplot(2,1,2); plt.plot(time,y[:,1::2]); plt.title('cSR'); plt.axis([0, 1000, 0, 10])
    plt.show()


  def display(self):
    widget = interact(self.solve_and_plot,
            units = IntSlider(value = 10, min = 1, max = 30, step= 1),
            tau_SR = IntSlider(value = 25, min = 10, max = 40, step= 1),
            tau_cyt = IntSlider(value = 25, min = 10, max = 40, step= 1),
            BC = IntSlider(value = 1, min = 1, max = 2, step= 1))


  def CICR(self, y, t):
    #constants
    k_1 = 2*10**(-5);
    k_2 = 0.13;
    k_4 = 0.9;
    kappa_1 = 0.013;
    kappa_2 = 0.58;
    K_d = 0.5;
    n = 3;
    gamma = 4.17;
    c0 = 1000;

    #input
    #input
    c = y[0::2];
    cSR = y[1::2];

    #calculation of k_3
    k_3 = kappa_1 + (kappa_2*c**n)/(K_d**n+c**n)

    #calcium entry
    J_L1 = k_1 * (c0 - c)

    #calcium extrusion
    J_P1 = k_2 * c;

    #calcium release
    J_L2  = k_3 *(cSR - c);

    #calcium uptake
    J_P2 = k_4 *c;

    # intra-unit fluxes
    tau_cyt_flux = ones(c.shape)*self.tau_cyt
    tau_SR_flux = ones(cSR.shape)*self.tau_SR

    if (self.units != 1):
        #Boundary Conditions
        Cai_BC = pad(c, (1, 1), mode='constant')
        Ca_sr_BC = pad(cSR, (1, 1), mode='constant')


        if self.BC == 2:
            Cai_BC[0] = 0.0795
            Cai_BC[self.units+1] = 0.0795
            Ca_sr_BC[0] = 4.1725
            Ca_sr_BC[self.units+1] = 4.1725

        #fluxes between same compartments in different units
        J_diff_cyt = (Cai_BC[0:self.units] + Cai_BC[2:self.units+2] - 2*c)/tau_cyt_flux
        J_diff_SR = (Ca_sr_BC[0:self.units] + Ca_sr_BC[2:self.units+2] - 2*cSR)/tau_SR_flux
    else:
        J_diff_cyt = 0
        J_diff_SR = 0

    #calculate time dependent functions
    dc = (J_L1-J_P1+J_L2-J_P2) + J_diff_cyt#dc
    dcSR = gamma*(-J_L2+J_P2) + J_diff_SR  #dcSR

    ydotMultiUnit = zeros(y.shape)
    ydotMultiUnit[0::2] = dc
    ydotMultiUnit[1::2] = dcSR


    return ydotMultiUnit
