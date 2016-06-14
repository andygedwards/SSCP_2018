%***********************************************************************%
%   Beeler-Reuter Ventricular Myocyte AP Model                          %
%   File:   BR.m                                                        %
%   Author: Stuart Campbell                                             %
%   Description: This script simulates the action potential of a        %
%   ventricular myocyte using the Beeler-Reuter Model (J Physiol, 1977. %
%   268(1): p.177-210).  Simulation time, time step size, and stimululs %
%   parameters can be set by the user.  The ODE solver used may also be %
%   changed by the user to note differences in the numberical solution. %
%***********************************************************************%

clear all, close all, clc

global STIMULUS;

tic
%-----------------------------------------------%
%-------------Simulation Parameters-------------%
%-----------------------------------------------%

period = 500.0;                     %[ms] - Total time of simulation
delta  = 0.03;                      %[ms] - Time step
steps  = floor(period / delta + 1); %     - Total number of time steps
t      = [0:delta:period];          %[ms] - Vector containing time points
method = 'stdRK4';                  %Solver option: 'stdRK4' or 'vtsRK4'

%---Stimulation Current Parameters---%
event   = [1:2];                    %[ms] - Stimulation event times
event_t = floor(event / delta);     %     - Convert to time step #'s
stim    = 500.0;                    %[uA] - Magnitude of current stimulation
i_stim  = zeros(1,steps);           %Initialize stim. current vector
i_stim(event_t) = stim;             %Set stim. event(s)


%-----------------------------------------------%
%-----------Initialize State Variables----------%
%-----------------------------------------------%

Vm   = zeros(1,steps);
Ca_i = zeros(1,steps);
x1   = zeros(1,steps);
m    = zeros(1,steps);
h    = zeros(1,steps);
j    = zeros(1,steps);
d    = zeros(1,steps);
f    = zeros(1,steps);
i_K1 = zeros(1,steps);
i_Na = zeros(1,steps);
i_Ca = zeros(1,steps);
i_x1 = zeros(1,steps);

Vm(1)   = -84.5732;          %Membrane Voltage
Ca_i(1) = 1.782e-7;          %Intracellular Calcium
x1(1)   = 0.0057;            %Gating parameters...
m(1)    = 0.011;
h(1)    = 0.9877;
j(1)    = 0.9748;
d(1)    = 0.003;
f(1)    = 1;
i_K1(1) = 0.4678;            %Ionic Currents...
i_Na(1) = -0.4044;
i_Ca(1) = -0.0547;
i_x1(1) = -0.0086;

%Set initial condition vector
p    = [Vm(1) Ca_i(1) x1(1) m(1) h(1) j(1) d(1) f(1) i_K1(1) i_Na(1) i_Ca(1) i_x1(1)]';

%-----------------------------------------------%
%----------------Call ODE Solvers---------------%
%-----------------------------------------------%

switch method
    case {'stdRK4'}   %Standard fixed time-step fourth order RK scheme
        
        %----------------Integration Loop---------------%
        for i = 1:(steps-1)
            new_p     = RK4('BR_deriv_std', t(i), delta, p, i_stim(i));
            Vm(i+1)   = new_p(1);
            Ca_i(i+1) = new_p(2);
            x1(i+1)   = new_p(3);
            m(i+1)    = new_p(4);
            h(i+1)    = new_p(5);
            j(i+1)    = new_p(6);
            d(i+1)    = new_p(7);
            f(i+1)    = new_p(8);
            i_K1(i+1) = new_p(9);
            i_Na(i+1) = new_p(10);
            i_Ca(i+1) = new_p(11);
            i_x1(i+1) = new_p(12);

            p = new_p;  %Update state variable vector
        end
        
    case {'vtsRK4'}
            a = 0;
            b = 0.2;
            tspan = [a b];
            STIMULUS = stim;
            [t1, new_p] = ode15s(@BR_deriv_vts, tspan, p);
            Vm   = new_p(:,1);
            Ca_i = new_p(:,2);
            x1   = new_p(:,3);
            m    = new_p(:,4);
            h    = new_p(:,5);
            j    = new_p(:,6);
            d    = new_p(:,7);
            f    = new_p(:,8);
            i_K1 = new_p(:,9);
            i_Na = new_p(:,10);
            i_Ca = new_p(:,11);
            i_x1 = new_p(:,12);
            
            tspan = [t1(end) period];
            STIMULUS = 0;
            p    = new_p(end,:);
            [t2, new_p] = ode45(@BR_deriv_vts, tspan, p);
            Vm_2   = new_p(:,1);
            Ca_i_2 = new_p(:,2);
            x1_2   = new_p(:,3);
            m_2    = new_p(:,4);
            h_2    = new_p(:,5);
            j_2    = new_p(:,6);
            d_2    = new_p(:,7);
            f_2    = new_p(:,8);
            i_K1_2 = new_p(:,9);
            i_Na_2 = new_p(:,10);
            i_Ca_2 = new_p(:,11);
            i_x1_2 = new_p(:,12);
end

toc
