function [q0, tCycle ] = Exercise_COHR_Relation( CO )
% Returns physiological values of flow (q0) and cycle duration (tCycle) for
% a given cardiac output and HR data that can serve as input for simulation
% of exercise
% 
% 07 August 2017, Joost Lumens and John Walmsley

% Linear fit on following data
% source: W.F. Boron, E.L. Boulpaep, Medical Physiology, 2003:

CO_data = [5,  7.5,  10,   12.5,  15,   17.5,  20,   22.5,  25,   27.5,  30];
HR_data = [71, 89,   106,  124,   141,  159,   176,  194,   211,  229,   246];

[fit_COHR]= fit(CO_data',HR_data','poly1');
slope = fit_COHR.p1;
interc = fit_COHR.p2;

HR = CO.*slope + interc; % Heart Rate [bpm]
q0 = CO./60000; % Cardiac Output [m3/s]
tCycle = 60./HR; % Cycle time [s]
