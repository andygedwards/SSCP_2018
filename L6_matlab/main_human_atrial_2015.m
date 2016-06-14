function [ output_args ] = main_human_atrial(initialcond, outputfile, tend, basic_cycle_length, pacinginit, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration)
%function [ output_args ] = main_human_atrial(stimulus, initialcond, outputfile, tend, basic_cycle_length, pacinginit, modelvariant, drug_compound, drug_concentration)
%function [ output_args ] = main_human_atrial(stimulus, initialcond, outputfile, tend, basic_cycle_length, pacinginit, modelvariant)

% -------------------------------------------------------------------------
% Inputs 
%
% stimulus = matrix, which has time in the first column and stimuluscurrent
% in the second column. If needed time point isn't found in the matrix,
% program interpolates to the previous timepoint. Stimcreator.m can be used
% to create stimulus timecourse for pacing experiments.
%
% initialcond = filename for initialconditions
%
% outputfile = name of the file where data is stored
%
% tend = tendtime
%
% BCL = pacing BCL used in simulation to notify the stimulus current
% during simulation (shorter max-step used during stimulation).
%
% pacinginit = initial values 0: use initialcond-file in pacing simulation or
%                             1: use quiecent-ss in pacing simulation 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OUTPUTS
% 
% Calculated data to files (outputfiless = steady-state-calculations, 
% outputfile = stimulation simulation) and calculation time displayed to workspace
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% CODE BEGINS
% -------------------------------------------------------------------------

clear BCL;
global BCL;
BCL = basic_cycle_length;

% Start timer
tic

% Create global matrix from stimulus
%clear stim;
%global stim;
%stim = stimulus;

% Output filename
global ofile;

% Load initial conditions
clear y0;
y0 = load(initialcond);

% Solve the problem

% No-stimulus steady-state calculations
if pacinginit
ofile = ['quiescent_' outputfile];
fid = fopen(ofile,'w');
fclose(fid);
options = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'OutputFcn',@tofile); 
[a, b] = ode15s(@dy_human_atrial, [0 4000], y0, options, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration);
end

% Stimulus calculations

% Solver options
% MaxStep has to be smaller than the length of stimulus during stimulation. Otherwise the solver doesn't notice the stimulus. 
if BCL > 0
    stimpoint = BCL/1000; % THIS ASSUMES THAT YOU ARE USING THIS APPROXIMATION IN STIMULUSCURRENT (e.g. stimcreator.m) TOO!!!
end
%options1 = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'MaxStep', 0.4e-3,'OutputFcn', @tofile); 
%options2 = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'MaxStep', 100,'OutputFcn', @tofile); 
options1 = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'MaxStep', 0.9e-3,'OutputFcn', @tofile); 
options2 = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'MaxStep', 100,'OutputFcn', @tofile); 

% Choose initial conditions
if pacinginit
    y0 = b(end, :); % Take values from the end of quiecent simulation
end

% Data reduction
global counter redu;
redu = 1; % every "redu"th point is saved 
counter = 1;

% Solve problem & store data to file
ofile = ['pacing_' outputfile];
fid = fopen(ofile,'w');
fclose(fid);

if BCL > 0
    if (round(tend/(BCL/1000))-1) >= 1
        for n = 1:((round(tend/(BCL/1000))-1))
            [a, b] = ode15s(@dy_human_atrial, [((n-1)*stimpoint) (((n-1)*stimpoint) + 1e-3)], y0, options1, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration);
            y0 = b(end,:);
            [a, b] = ode15s(@dy_human_atrial, [(((n-1)*stimpoint) + 1e-3) (n*stimpoint)], y0, options2, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration);
            y0 = b(end,:);
        end
        redu = 1;
        [a, b] = ode15s(@dy_human_atrial, [(n*stimpoint) tend], y0, options1, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration);
    else
        redu = 1;
        [a, b] = ode15s(@dy_human_atrial, [0 tend], y0, options1, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration);
    end
else
[a, b] = ode15s(@dy_human_atrial, [0 tend], y0, options2, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration);
end
    

% Stop timer
toc

% -------------------------------------------------------------------------
% SUBFUNCTIONS BEGIN
% -------------------------------------------------------------------------

function status = tofile(t,y,flag,varargin)
% Stores data to file, no headers & doesn't erase old data in file
% Stores some differential variables and some ion currents & fluxes

% Modified from odeprint by Topi Korhonen
% Original function by
%   Mark W. Reichelt and Lawrence F. Shampine, 3-24-94
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.17 $  $Date: 2002/04/08 20:26:49 $

global fid ofile; % Pointer to file is available for fprintf when flag is 'empty'

% Data reduction
global counter redu

% % Functions which values will be saved to file
global Istim INa ICaL It Isus IKr IKs IK1 IKACh IKCa If INab ICab ICaP INaK INaCa Jrelss Jrel1 Jrel2 Jrel3 J_bulkSERCA1 J_bulkSERCA2 J_bulkSERCA3 J_bulkSERCAss JSRCaleak1 JSRCaleak2 JSRCaleak3 JSRCaleakss

i_SIZE = length(y);

% Create outputstyle
style = '';

%for n = 1:(1 + 27 + 4 + 3 + 1 + 4*2 + 17 + 4*2) % time, differential variables and functions
%for n = 1:(1 + 27 + 4 + 3 + 1 + 4*2 + 17 + 4*2 + 1*2) % time, differential variables and functions + gKCa&IKCa
for n = 1:(1 + 28+4+3+1+2*4 + 16 + 3*4) % time, differential variables, currents and fluxes
    style = [style ' %6.6e'];
end

style = [style '\n'];

% Outputs to file
if nargin < 3 | isempty(flag) % This is done after every succesful integration step
    if mod(counter,redu) == 0
        fprintf(fid, style, [ t y' Istim INa ICaL It Isus IKr IKs IK1 IKACh IKCa If INab ICab ICaP INaK INaCa Jrelss Jrel1 Jrel2 Jrel3 J_bulkSERCA1 J_bulkSERCA2 J_bulkSERCA3 J_bulkSERCAss JSRCaleak1 JSRCaleak2 JSRCaleak3 JSRCaleakss]);
        counter = counter + 1;
    else
        counter = counter + 1;
    end
else
    switch(flag)
        case 'init'               % This is done in the initial point
            % Create file for outputdata
             fid = fopen(ofile,'a');
        case 'done'               % This is done when solving ends
             fclose(fid);
    end
end

status = 0;
