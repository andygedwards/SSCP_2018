%***********************************************************************%
%   Fourth-Order Runge Kutta ODE Solver                                 %
%   File:   RK4.m                                                       %
%   Author: Stuart Campbell                                             %
%   Description: This function approximates the value of a state        %
%   variable (or vector or state variables) p at time = t + delta by    %
%   repeatedly calling a separate function of the form dp/dt = f(t,p).  %
%   Inputs:  deriv_func --> A string containing the name of the         %
%                           derivative function                         %
%            t          --> Time
%            delta      --> Size of the time step                       %
%            p          --> Value of state var's at time t              %
%            Param      --> Array containing model parameters to be     %
%                           passed to the derivative function           %
%   Output:  new_p      --> Value of state var's at time t + delta      %
%***********************************************************************%

function new_p= RK4(deriv_func, t, delta, p, Param)
p1= p;
   k1= feval(deriv_func,t,p1,Param);
p2= p + k1*delta/2;
   k2= feval(deriv_func,t,p2,Param);
p3= p + k2*delta/2;
   k3= feval(deriv_func,t,p3,Param);
p4= p + k3*delta;
   k4= feval(deriv_func,t,p4,Param);
new_p= p + (k1 + 2*k2 + 2*k3 + k4)*delta/6;
return