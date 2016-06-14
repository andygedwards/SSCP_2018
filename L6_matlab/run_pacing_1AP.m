%run_pacing_1AP
%hold on
%(change model)
%run_pacing_1AP

y0 = load('y0_ss_BCL_1000.dat');
global BCL
BCL = 1000;
%dy_human_atrial(t, y, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration)
%[t,y] = ode15s(@dy_human_atrial, [0,1], y0, [], 'stim_I', 10 ,'cAF_all', 'control', 0);
[t,y] = ode15s(@dy_human_atrial, [0,1], y0, [], 'stim_I', 10 ,'nSR', 'control', 0);
plot(t,y(:,(1:1)))





