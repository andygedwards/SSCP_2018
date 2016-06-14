function [] = run_pacing_1AP(BCL, modelvariant, drug_compound, drug_concentration)
%% Run a single AP. For example,
%  >> run_pacing_1AP_Ca(1000, 'nSR', 'control', 0)
%run_pacing_1AP_Ca 
%(change model)
%run_pacing_1AP)Ca



stimulus_type = 'stim_I';
stimulus_value = 0;

sim_length = BCL/1000;

if BCL < 1000
    main_human_atrial_2015(['y0_ss_BCL_0' num2str(BCL) '.dat'], ['1AP_BCL_0' num2str(BCL) '.dat'], sim_length, BCL, 0, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration);
    results_ss_human_atrial(['Results_BCL_0' num2str(BCL)], ['pacing_1AP_BCL_0' num2str(BCL) '.dat']);
else
    main_human_atrial_2015(['y0_ss_BCL_' num2str(BCL) '.dat'], ['1AP_BCL_' num2str(BCL) '.dat'], sim_length, BCL, 0, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration);
    results_ss_human_atrial(['Results_BCL_' num2str(BCL)], ['pacing_1AP_BCL_' num2str(BCL) '.dat']);
end
