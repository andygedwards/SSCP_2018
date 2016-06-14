%*************************************************************************%
%   Beeler-Reuter Ventricular Myocyte AP Model                            %
%   File:   BR_deriv.m                                                    %
%   Author: Stuart Campbell                                               %
%   Description: This function returns the values of the time derivatives %
%   of state variables for use in numerical integration of the solution.  %
%   Most fixed model parameters are contained within this function.       %
%*************************************************************************%

function pdot = BR_deriv(t,p,i_stim)

%-----------------------------------------------%
%----------------Extract St. Vars---------------%
%-----------------------------------------------%

Vm   = p(1);     %[mV] - Transmembrane Voltage
Ca_i = p(2);     %[mol] - Intracellular Calcium

%These six are unitless channel gating parameters
x1 = p(3);     
m  = p(4);
h  = p(5);
j  = p(6);
d  = p(7);
f  = p(8);

%These four are ionic currents - not state var's; for output purposes only
i_K1old = p(9);
i_Naold = p(10);
i_Caold = p(11);
i_x1old = p(12);

%-----------------------------------------------%
%----------------Model Parameters---------------%
%-----------------------------------------------%

E_Na  = 50.0;        %[mV]        - Nernst potential from Na
g_Na  = 4.0;         %[mmho/cm^2] - Membrane conductance parameter
g_NaC = 0.003;       %[mmho/cm^2] - Membrane conductance parameter
g_s   = 0.09;        %[mmho/cm^2] - Membrane conductance parameter
C_m   = 1.0;         %[uF/cm^2]   - Membrane capacitance

%Rate constant parameters for determining alpha and beta
C(:,:,1) = [5e-4   8.3e-2  50  0   0   5.7e-2  1; ...
            0      0       47  -1  47  -0.1   -1; ...
            0.126  -0.25   77  0   0   0       0; ...
            5.5e-2 -0.25   78  0   0   -0.2    1; ...
            9.5e-2 -0.01   -5  0   0   -0.072  1; ...
            1.2e-2 -8.0e-3 28  0   0   0.15    1]; 
     
C(:,:,2) = [1.3e-3 -6.0e-2 20  0   0   -4.0e-2 1; ...
            40     -5.6e-2 72  0   0   0       0; ...
            1.7    0      22.5 0   0   -8.2e-2 1; ...
            0.3    0       32  0   0   -0.1    1; ...
            7.0e-2 -1.7e-2 44  0   0   5.0e-2  1; ...
            6.5e-3 -2.0e-2 30  0   0   -0.2    1];

%-----------------------------------------------%
%------------Preliminary Calculations-----------%
%-----------------------------------------------%

%Calculate alpha and beta values
a_b = (C(:,1,:) .* exp(C(:,2,:) .* (Vm + C(:,3,:))) + C(:,4,:) .* (Vm + C(:,5,:)))...
            ./(exp(C(:,6,:) .* (Vm + C(:,3,:))) + C(:,7,:));

%Calculate y_inf and tau
y_inf = a_b(:,:,1)./(a_b(:,:,1) + a_b(:,:,2));
tau   = (a_b(:,:,1) + a_b(:,:,2)).^-1;

E_Ca  = -82.3 - 13.0287 * log(Ca_i);

%-----------------------------------------------%
%----------Ionic Current Calculations-----------%
%-----------------------------------------------%

% %i_K1 = 0.1 * (4 * (exp(0.04 * (Vm + 85)) - 1) / (exp(0.08 * (Vm + 53))...
%      + exp(0.04 * (Vm + 53))) + 0.2 * ((Vm + 23) / (1 - exp(-0.04 * (Vm + 23)))));
i_K1 = 0.35 * (4 * (exp(0.04 * (Vm + 85)) - 1) / (exp(0.08 * (Vm + 53))...
     + exp(0.04 * (Vm + 53))) + 0.2 * ((Vm + 23) / (1 - exp(-0.04 * (Vm + 23)))));
i_x1 = (x1) * 0.8 * (exp(0.04 * (Vm + 77)) - 1)/exp(0.04 * (Vm + 35));
i_Na = (g_Na * m^3 * h * j + g_NaC) * (Vm - E_Na);
i_Ca = g_s * d * f * (Vm - E_Ca);


%-----------------------------------------------%
%------------Derivative Calculations------------%
%-----------------------------------------------%

pdot(1,1)   = -(1/C_m) * (i_K1 + i_x1 + i_Na + i_Ca - i_stim);  %dVm/dt
pdot(2,1)   = -1e-7 * i_Ca + 0.07 * (1e-7 - Ca_i);              %d[Ca]/dt
pdot(3:8,1) = (y_inf - p(3:8)) ./ tau;
pdot(9,1)   = i_K1 - i_K1old;
pdot(10,1)  = i_Na - i_Naold;
pdot(11,1)  = i_Ca - i_Caold;
pdot(12,1)  = i_x1 - i_x1old;