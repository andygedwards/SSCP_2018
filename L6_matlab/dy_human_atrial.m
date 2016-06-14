function [dy] = dy_human_atrial(t, y, stimulus_type, stimulus_value, modelvariant, drug_compound, drug_concentration)
%**************************************************************************
% Human atrial myocyte model
%
% Citation: N/A
%
% -------------------------------------------------------------------------
% - Modifications to 2014 version:
%   * reformulated sodium current, I_Na (including late current)
%   * small conductance calcium activated potassium current, I_KCa
%   * included, I_KACh as in Voigt et al. (2013)
%   * updated I_Kr and I_Ks formulations as in Grandi et al. (2011)
%   * updated NKA formulation as in Grandi et al. (2011)
%   * increased I_to conductance and reformulated I_CaL to get a more
%     spike-and-dome like AP shape
%   * dependence of ECa_app on Ca_ext as in Campbell et al. (1988)
%   * temperature adjusted to 37C (from 33C)
%   * cAF model variant extended to include:
%     > reduced I_Na
%     > reduced I_KCNa
%     > increased I_Ks
% 
% -------------------------------------------------------------------------
% Model pedigree:
%
% Koivumaeki JT, Seemann G, Maleckar MM, Tavi P (2014)
% In Silico Screening of the Key Cellular Remodeling Targets in Chronic
% Atrial Fibrillation.
% PLoS Comput Biol 10(5): e1003620.
% http://dx.doi.org/10.1371/journal.pcbi.1003620
%
% - Modifications to 2011 version:
%   * reformulated L-type calcium current
%   * SERCA pump with explicit PLB and SLN effects
%   * corrected Ito and Isus conductances
%
% Koivumaeki JT, Korhonen T, Tavi P (2011) 
% Impact of Sarcoplasmic Reticulum Calcium Release on Calcium Dynamics and
% Action Potential Morphology in Human Atrial Myocytes: A Computational Study.
% PLoS Comput Biol 7(1): e1001067.
% http://dx.doi.org/10.1371/journal.pcbi.1001067
%
% - Model developed based on:
%   * Nygren et al. (1998) model, including
%   * updated K+ current formulations as in Maleckar et al. (2009)
%
%**************************************************************************

%**********************************************
% Globals
%**********************************************

global BCL
% LIS???? T??H??N LISTAAN INaL
global Istim INa ICaL It Isus IKr IKs IK1 IKACh IKCa If INab ICab ICaP INaK INaCa Jrelss Jrel1 Jrel2 Jrel3 J_bulkSERCA1 J_bulkSERCA2 J_bulkSERCA3 J_bulkSERCAss JSRCaleak1 JSRCaleak2 JSRCaleak3 JSRCaleakss

%**********************************************
% Definition of differential variables
%**********************************************

i_V = 1;
i_lastend = i_V;

% INa
i_start = i_lastend + 1;
i_INam = i_start; i_INah1 = i_start + 1; i_INah2 = i_start + 2;
i_lastend = i_INah2;

% ICaL 
i_start = i_lastend + 1;
i_ICaLd = i_start; i_ICaLf1 = i_start + 1; i_ICaLf2 = i_start + 2; i_ICaLfca = i_start + 3;
i_lastend = i_ICaLfca;

% It
i_start = i_lastend + 1;
i_Itr = i_start; i_Its = i_start + 1;
i_lastend = i_Its;

% Isus (Ikur)
i_start = i_lastend + 1;
i_Isusr = i_start; i_Isuss = i_start + 1;
i_lastend = i_Isuss;

% IKs
i_start = i_lastend + 1;
i_IKsn = i_start;
i_lastend = i_IKsn;

% IKr
i_start = i_lastend + 1;
i_IKrpa = i_start;
i_lastend = i_IKrpa;

% If
i_start = i_lastend + 1;
i_Ify = i_start;
i_lastend = i_Ify;

% gKCa
i_start = i_lastend + 1;
i_IKCa_O = i_start;
i_lastend = i_IKCa_O;

% RyR
i_start = i_lastend + 1;
i_RyRoss = i_start; i_RyRcss = i_start + 1; i_RyRass = i_start + 2;
i_RyRo1 = i_start + 3; i_RyRc1 = i_start + 4; i_RyRa1 =  i_start + 5;
i_RyRo2 = i_start + 6; i_RyRc2 = i_start + 7; i_RyRa2 =  i_start + 8;
i_RyRo3 = i_start + 9; i_RyRc3 = i_start + 10; i_RyRa3 =  i_start + 11;
i_lastend = i_RyRa3;

% SERCA
i_start = i_lastend + 1;
i_SERCACa1 = i_start; % first compartment from center
i_SERCACa2 = i_start + 1; % 2nd compartment from center
i_SERCACa3 = i_start + 2; % 3rd compartment from center
i_SERCACass = i_start + 3; % ss-SERCA
i_lastend = i_SERCACass;

% Nai ja Ki
i_start = i_lastend + 1;
i_Nass = i_start; i_Nai = i_start + 1;
i_Ki = i_start + 2;
i_lastend = i_Ki;

% Cai and CaSR 
i_start = i_lastend + 1;
i_Cass = i_start;
i_Cacenter = i_start + 1;
%i_lastend = i_Cacenter;

%**********************************************
%% Select model variant
%**********************************************

switch modelvariant

    case 'cAF_all'
        cAF_lcell = 1.10; cAF_gNa = 0.82; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gKs = 2.70; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2.00; cAF_RyR = 2.00; cAF_gKCa = 0.50;

    case 'cAF_dilation'
        cAF_lcell = 1.10; cAF_gNa = 1.00; cAF_gCaL = 1.00; cAF_gt = 1.00; cAF_gsus = 1.00; cAF_gKs = 1.00; cAF_gK1 = 1.00; cAF_kNaCa = 1.00; cAF_cpumps = 1.00; cAF_PLB = 1.00; cAF_SLN = 1.00; cAF_phos = 1.00; cAF_RyR = 1.00; cAF_gKCa = 1.00;
    case 'cAF_ICaL'
        cAF_lcell = 1.00; cAF_gNa = 1.00; cAF_gCaL = 0.41; cAF_gt = 1.00; cAF_gsus = 1.00; cAF_gKs = 1.00; cAF_gK1 = 1.00; cAF_kNaCa = 1.00; cAF_cpumps = 1.00; cAF_PLB = 1.00; cAF_SLN = 1.00; cAF_phos = 1.00; cAF_RyR = 1.00; cAF_gKCa = 1.00;
    case 'cAF_IK1'
        cAF_lcell = 1.00; cAF_gNa = 1.00; cAF_gCaL = 1.00; cAF_gt = 1.00; cAF_gsus = 1.00; cAF_gKs = 1.00; cAF_gK1 = 1.62; cAF_kNaCa = 1.00; cAF_cpumps = 1.00; cAF_PLB = 1.00; cAF_SLN = 1.00; cAF_phos = 1.00; cAF_RyR = 1.00; cAF_gKCa = 1.00;
    case 'cAF_IKs'
        cAF_lcell = 1.00; cAF_gNa = 1.00; cAF_gCaL = 1.00; cAF_gt = 1.00; cAF_gsus = 1.00; cAF_gKs = 2.70; cAF_gK1 = 1.00; cAF_kNaCa = 1.00; cAF_cpumps = 1.00; cAF_PLB = 1.00; cAF_SLN = 1.00; cAF_phos = 1.00; cAF_RyR = 1.00; cAF_gKCa = 1.00;
    case 'cAF_INa'
        cAF_lcell = 1.00; cAF_gNa = 0.82; cAF_gCaL = 1.00; cAF_gt = 1.00; cAF_gsus = 1.00; cAF_gKs = 1.00; cAF_gK1 = 1.00; cAF_kNaCa = 1.00; cAF_cpumps = 1.00; cAF_PLB = 1.00; cAF_SLN = 1.00; cAF_phos = 1.00; cAF_RyR = 1.00; cAF_gKCa = 1.00;
    case 'cAF_Isus'
        cAF_lcell = 1.00; cAF_gNa = 1.00; cAF_gCaL = 1.00; cAF_gt = 1.00; cAF_gsus = 0.62; cAF_gKs = 1.00; cAF_gK1 = 1.00; cAF_kNaCa = 1.00; cAF_cpumps = 1.00; cAF_PLB = 1.00; cAF_SLN = 1.00; cAF_phos = 1.00; cAF_RyR = 1.00; cAF_gKCa = 1.00;
    case 'cAF_Ito'
        cAF_lcell = 1.00; cAF_gNa = 1.00; cAF_gCaL = 1.00; cAF_gt = 0.38; cAF_gsus = 1.00; cAF_gKs = 1.00; cAF_gK1 = 1.00; cAF_kNaCa = 1.00; cAF_cpumps = 1.00; cAF_PLB = 1.00; cAF_SLN = 1.00; cAF_phos = 1.00; cAF_RyR = 1.00; cAF_gKCa = 1.00;
    case 'cAF_NCX'
        cAF_lcell = 1.00; cAF_gNa = 1.00; cAF_gCaL = 1.00; cAF_gt = 1.00; cAF_gsus = 1.00; cAF_gKs = 1.00; cAF_gK1 = 1.00; cAF_kNaCa = 1.50; cAF_cpumps = 1.00; cAF_PLB = 1.00; cAF_SLN = 1.00; cAF_phos = 1.00; cAF_RyR = 1.00; cAF_gKCa = 1.00;
    case 'cAF_RyR'
        cAF_lcell = 1.00; cAF_gNa = 1.00; cAF_gCaL = 1.00; cAF_gt = 1.00; cAF_gsus = 1.00; cAF_gKs = 1.00; cAF_gK1 = 1.00; cAF_kNaCa = 1.00; cAF_cpumps = 1.00; cAF_PLB = 1.00; cAF_SLN = 1.00; cAF_phos = 1.00; cAF_RyR = 2; cAF_gKCa = 1.00;
    case 'cAF_SERCA'
        cAF_lcell = 1.00; cAF_gNa = 1.00; cAF_gCaL = 1.00; cAF_gt = 1.00; cAF_gsus = 1.00; cAF_gKs = 1.00; cAF_gK1 = 1.00; cAF_kNaCa = 1.00; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2; cAF_RyR = 1.00; cAF_gKCa = 1.00;

    case 'cAF_no_dilation'
        cAF_lcell = 1.00; cAF_gNa = 0.82; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gKs = 2.70; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2.00; cAF_RyR = 2.00; cAF_gKCa = 0.50;
    case 'cAF_no_ICaL'
        cAF_lcell = 1.10; cAF_gNa = 0.82; cAF_gCaL = 1.00; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gKs = 2.70; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2.00; cAF_RyR = 2.00; cAF_gKCa = 0.50;
    case 'cAF_no_IK1'
        cAF_lcell = 1.10; cAF_gNa = 0.82; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gKs = 2.70; cAF_gK1 = 1.00; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2.00; cAF_RyR = 2.00; cAF_gKCa = 0.50;
    case 'cAF_no_IKs'
        cAF_lcell = 1.10; cAF_gNa = 0.82; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gKs = 1.00; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2.00; cAF_RyR = 2.00; cAF_gKCa = 0.50;
    case 'cAF_no_INa'
        cAF_lcell = 1.10; cAF_gNa = 1.00; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gKs = 2.70; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2.00; cAF_RyR = 2.00; cAF_gKCa = 0.50;
    case 'cAF_no_Isus'
        cAF_lcell = 1.10; cAF_gNa = 0.82; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 1.00; cAF_gKs = 2.70; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2.00; cAF_RyR = 2.00; cAF_gKCa = 0.50;
    case 'cAF_no_Ito'
        cAF_lcell = 1.10; cAF_gNa = 0.82; cAF_gCaL = 0.41; cAF_gt = 1.00; cAF_gsus = 0.62; cAF_gKs = 2.70; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2.00; cAF_RyR = 2.00; cAF_gKCa = 0.50;
    case 'cAF_no_NCX'
        cAF_lcell = 1.10; cAF_gNa = 0.82; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gKs = 2.70; cAF_gK1 = 1.62; cAF_kNaCa = 1.00; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2.00; cAF_RyR = 2.00; cAF_gKCa = 0.50;
    case 'cAF_no_RyR'
        cAF_lcell = 1.10; cAF_gNa = 0.82; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gKs = 2.70; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 0.84; cAF_PLB = 1.18; cAF_SLN = 0.60; cAF_phos = 2.00; cAF_RyR = 1.00; cAF_gKCa = 0.50;
    case 'cAF_no_SERCA'
        cAF_lcell = 1.10; cAF_gNa = 0.82; cAF_gCaL = 0.41; cAF_gt = 0.38; cAF_gsus = 0.62; cAF_gKs = 2.70; cAF_gK1 = 1.62; cAF_kNaCa = 1.50; cAF_cpumps = 1.00; cAF_PLB = 1.00; cAF_SLN = 1.00; cAF_phos = 1.00; cAF_RyR = 2.00; cAF_gKCa = 0.50;
        
    case 'nSR'
        cAF_lcell = 1; cAF_gNa = 1; cAF_gCaL = 1; cAF_gt = 1; cAF_gsus = 1; cAF_gKs = 1; cAF_gK1 = 1; cAF_kNaCa = 1; cAF_cpumps = 1; cAF_PLB = 1; cAF_SLN = 1; cAF_phos = 1; cAF_RyR = 1; cAF_gKCa = 1;
        
    otherwise
        disp('Unknown model variant.')
end

%**********************************************
%% Select stimulus
%**********************************************

switch stimulus_type

    case 'stim_AP_clamp'
        % Set the two-column stim vector global also in the script.
        global stim;
        st = stim( find(stim(:,1) <= t & stim(:,1) >= 0) , 2 ); % Find all (time, current)-data which t is smaller or equal to t
        if length(st) == 0
            Istim = 0;
        else
            Istim = -(st(end) - y(i_V)) / 0.0001; % for voltage clamp simulation, according to rate = 0.001
        end
        
    case 'stim_I_old'
        % This uses previous datarow in stim-matrix for the value of stimuluscurrent at time t. In other words: interpolating to previous time-point.
        st = stim( find(stim(:,1) <= t & stim(:,1) >= 0) , 2 ); % Find all (time, current)-data which t is smaller or equal to t
        if length(st) == 0
            Istim = 0;
        else
            Istim = st(end); % use the last value (highest t) % for normal pacing
        end

    case 'stim_I'
        stim_offset = 0.0; % has to be smaller than the max step in solver
        stim_duration = 0.002;
        stim_amplitude_multi = 1.25; % standard
        stim_amplitude = -460*stim_amplitude_multi * cAF_gK1; % new model
        %stim_amplitude = -950*stim_amplitude_multi * cAF_gK1; % LuoRudy
        %stim_amplitude = -675*stim_amplitude_multi * cAF_gK1; % Nygren
        %stim_amplitude = -240*stim_amplitude_multi * cAF_gK1; % TenTusscher
        stim_period = BCL/1000;
        if ((t-floor(t/stim_period)*stim_period >= stim_offset) && (t-floor(t/stim_period)*stim_period <= stim_offset + stim_duration))
           Istim = stim_amplitude;
        else
           Istim = 0.0;
        end
        
    case 'stim_I_adjustable_duration_and_amplitude'
        stim_offset = 0.0; % has to be smaller than the max step in solver
        stim_duration = stimulus_value(1);
        stim_amplitude = stimulus_value(2);
        stim_period = BCL/1000;
        if ((t-floor(t/stim_period)*stim_period >= stim_offset) && (t-floor(t/stim_period)*stim_period <= stim_offset + stim_duration))
           Istim = stim_amplitude;
        else
           Istim = 0.0;
        end
        
    case 'stim_I_adjustable_baseline'
        stim_offset = 0.0; % has to be smaller than the max step in solver
        stim_duration = 0.002; stim_amplitude = -460*2; % 2-times the threshold
        stim_baseline = stimulus_value(1);
        stim_period = BCL/1000;
        if ((t-floor(t/stim_period)*stim_period >= stim_offset) && (t-floor(t/stim_period)*stim_period <= stim_offset + stim_duration))
           Istim = stim_amplitude + stim_baseline;
        else
           Istim = stim_baseline;
        end
        
    case 'stim_I_induce_to_pace'
        stim_offset = 0.100; stim_duration = 1.800;
        stim_amplitude = -30/4; % outward hump of IK1 is about 30 pA.
        stim_period = BCL/1000;
        if ((t-floor(t/stim_period)*stim_period >= stim_offset) && (t-floor(t/stim_period)*stim_period <= stim_offset + stim_duration))
           Istim = stim_amplitude;
        else
           Istim = 0.0;
        end
        
    case 'stim_VC'
        stim_offset = 0.5; stim_duration = 0.5;
        stim_Vhold = stimulus_value(1); stim_Vtest = stimulus_value(2);
        if ((t >= stim_offset) && (t <= stim_offset + stim_duration))
            Istim = -(stim_Vtest - y(i_V)) / 0.0001;
        else
            Istim = -(stim_Vhold - y(i_V)) / 0.0001;
        end

    case 'stim_VC_INaL'
        % As in Olesen et al. (2012) for measuring INa late vs. peak
        stim_offset = 0.5; stim_duration = 0.5;
        stim_Vhold = stimulus_value(1); stim_Vtest = stimulus_value(2);
        if ((t >= stim_offset) && (t <= stim_offset + stim_duration))
            Istim = -(stim_Vtest - y(i_V)) / 0.0001;
        else
            Istim = -(stim_Vhold - y(i_V)) / 0.0001;
        end

    case 'stim_VC_train'
        stim_offset = 0;
        stim_duration = 0.050;
        stim_Vhold = stimulus_value(1);
        stim_V_increment = stimulus_value(2);
        stim_period = BCL/1000;
        
        if ((t-floor(t/stim_period)*stim_period >= stim_offset) && (t-floor(t/stim_period)*stim_period <= stim_offset + stim_duration))
            Istim = -(stim_Vhold + stim_V_increment*(t - mod(t, stim_period))/stim_period - y(i_V)) / 0.0001;
        else
            Istim = -(stim_Vhold - y(i_V)) / 0.0001;
        end

    case 'stim_VC_INa_recovery'
        % as in Schneider et al. (1994)
        stim_offset = 0.5;
        stim_duration = 0.2;
        stim_Vtest = -20;
        stim_Vhold = stimulus_value(1);
        stim_int = stimulus_value(2);
        if ((t >= stim_offset) && (t <= stim_offset + stim_duration))
            Istim = -(stim_Vtest - y(i_V)) / 0.0001;
        elseif ((t >= stim_offset + stim_duration + stim_int) && (t <= stim_offset + stim_duration + stim_int + stim_duration))
            Istim = -(stim_Vtest - y(i_V)) / 0.0001;
        else
            Istim = -(stim_Vhold - y(i_V)) / 0.0001;
        end

    case 'stim_VC_ramp'
        stim_offset = 1.000;
        stim_duration = 1.000;
        stim_V_increment = 190.0;
        stim_period = BCL;
        stim_Vhold = -80;
        stim_Vmin = -140;
        
        if ((t-floor(t/stim_period)*stim_period >= stim_offset) && (t-floor(t/stim_period)*stim_period <= stim_offset + stim_duration))
            Istim = -(stim_Vmin + (t - stim_offset) * stim_V_increment - y(i_V)) / 0.003;
        else
            Istim = -(stim_Vhold - y(i_V)) / 0.003;
        end

    otherwise
        disp('Unknown stimulus type.')
end

%**********************************************
%% Select drug or control
%**********************************************
% INa INaL ICaL It Isus IKr IKs IK1 IKACh IKCa If INaK INaCa
drug_vector = ones(13,1);

switch drug_compound
    
    case 'control'
        ACh = 0.0;
        
    case 'Isoprenaline'
        ISO = 0;
    case 'Acetylcholine'
        ACh = drug_concentration;
        
    case 'E4031'
        drug_vector(6) = 1/(1 + drug_concentration/20e-9);
        ACh = 0.0;
    case 'Barium'
        drug_vector(8) = 0.85/(1 + drug_concentration/2.412) + 0.15; %uM
        ACh = 0.0;
        
    case 'IKr_block'
        drug_vector(6) = 0;
        ACh = 0.0;
    case 'IKs_block'
        drug_vector(7) = 0;
        ACh = 0.0;
    case 'IKCa_block'
        drug_vector(10) = 0;
        ACh = 0.0;
    case 'INKA_block'
        drug_vector(11) = 0;
        ACh = 0.0;
    case 'INCX_block'
        drug_vector(12) = 0;
        ACh = 0.0;
                
    otherwise
        disp('Unknown compound.')
end

%**********************************************
%% Numerical parameters
%**********************************************

% Physical & environmental constants
F = 96487;
R = 8314;
T = 37 + 273.15; % 37C
%T = 24 + 273.15; % Schneider et al. INa recovery, Torsten's INa I-V curve
%T = 21 + 273.15; % Olesen et al. INaL

%Nao = 5; Nass = 5; % Torsten's INa solutions
Nao = 130; Nass = y(i_Nass); % normal solutions

Cao = 1.8;
%global Caext; Cao = Caext;

%Ko = 4.5;
Ko = 5.4;
Ki = y(i_Ki);

Cm = 0.05; %nF

% Temperature dependencies
q10exp = (T - 310.15)/10;
q10_DCaNa = 1.18;
q10_DCaBmSR = 1.425; % Sidell & Hazel (1987)
%q10_ICaL = 2.6;
q10_INa = 3; % Benndorf (1994)
INa_V_shift_act = -4.65 * q10exp; % linear temperature-dependent shift
INa_V_shift_inact = -7.85 * q10exp; % linear temperature-dependent shift, not including Yuan et al.
q10_Ito = 2.6; % Radicke et al. (2013): 2.4/2.8 for act/inact
q10_Isus = 2.2;
q10_SERCA = 2;
q10_RyR = 1.5;
q10_NKA = 1.63; %q10_KmNai = 1.49;
q10_CaP = 2.35;
q10_NCX = 1.57;

% Cell dilation in cAF
Ddcell = (cAF_lcell - 1)*(20/10) + 1;
Dvcell = cAF_lcell*Ddcell^2;

% Geometry
Vss = 4.99232e-5 * Dvcell; %nL
rjunct = 6.5 * Ddcell; % mum
lcell = 122.051 * cAF_lcell; % mum

% Ca diffusion grid
dx = 1.625 * Ddcell; % mum
rstart = 0 + 0.5*dx;
rend = rjunct - 0.5*dx;
j = round(rstart/dx):1:round(rjunct/dx); % Spatial index of Cai diffusion
j = j';

Aj_nj = pi*rjunct*2*lcell*0.5; % Area between junct and nonjunct
xj_nj = 0.02/2 * Ddcell + dx/2; % diffusion distance from center to center of junct to first njunct
xj_nj_Nai = 0.02/2 * Ddcell + 2*dx; % diffusion distance from center of junct to center of njunct (between 2nd and 3rd njunct)

% Diffusion compartment volumes
Vnonjunct = zeros(length(j),1);
Vnonjunct = (pi.*(j.*dx).^2.*lcell-pi.*((j-1).*dx).^2.*lcell).*1e-6.*0.5; %nL

Vcytosol = sum(Vnonjunct) + Vss;

VSR = 0.05.*Vnonjunct./2*0.9 / Dvcell;

Vnonjunct_Nai = sum(Vnonjunct);

% Non-junct Cai data & CaSR data
Cai = zeros(length(j),1);
Cai = y(i_Cacenter:length(j)+i_Cacenter-1);
CaSR = y(length(j)+i_Cacenter:length(j)*2+i_Cacenter-1);

% Cytosol Ca Buffers
Begta = 0;
%Begta = 5; % 5 mM EGTA
%Begta = 10; % 10 mM EGTA
Bcmdn = 24e-3;
BCa = Bcmdn + Begta;
%BCa = 24e-3;
SLlow = 165;
SLhigh = 13;

KdBegta = 0.12e-3;
KdBcmdn = 2.38e-3;
KdBCa = (Bcmdn*KdBcmdn + Begta*KdBegta) / (Bcmdn + Begta);
KdSLlow = 1.1;
KdSLhigh = 13e-3;

% SR Ca buffers
CSQN =  6.7;
KdCSQN = 0.8;

% Sarcolemmal Na burrering
% Area relation to Grandi et al. model 6.5^2 * pi() * 122 / (10.25^2 * pi() * 100) = 0.49
% Bmax = 7.651*0.11 + 1.65*0.89 = 2.31
BNa = 0.49 * 2.31;
KdBNa = 10;

%% Ion channel conductances & permeabilities & other parameters
gNaL = 0.47*0.9; % accounting for the 10% loss of INa in enzymatic cell isolation

%ECa_app = 60; kCa = 0.6e-3; gCaL = 15 * cAF_gCaL; % PLoS 2014
ECa_app = 60 + 29.2*log10(Cao/1.8); % Campbell et al. (1988)
%gCaL = 6.75; % Nygren
%gCaL = 6.2; % Courtemanche
kCa = 0.6e-3; gCaL = 7.0 * cAF_gCaL;

%gt = 8.25 * cAF_gt; % Maleckar
gt = 11 * cAF_gt; % increased ~33%
gsus = 2.25 * cAF_gsus; % Maleckar
%gKs = 1; % Nygren
% Grandi: 0.0035*50 = 0.175
% Courtemanche: 0.0294*50 = 1.47
gKs = 0.175 * drug_vector(7) * cAF_gKs; % Grandi
%gKr = 0.50; % PLoS 2014
gKr = 3.4 * drug_vector(6);
gK1 = 2.9 * cAF_gK1 * drug_vector(8); % further reduced with IKACh and IKCa, and increased IKr

EC50_ACh = 0.125e-3; % Voigt et al. (2013)
hill_ACh = 1; % Voigt et al. (2013)
gKACh = 6.375 * 1 / (1 + (EC50_ACh/ACh)^hill_ACh); % Voigt et al. (2013) 0.1275*50

gKCa = 3.7 * drug_vector(10) * cAF_gKCa;

gNab = 0.060599;
%gCab = 0.0952;
gCab = 0.084; % to reduce SR Ca load

% = 70.8253; kNaKK = 1; kNaKNa = 11; % Nygren
IbarNaK = q10_NKA^q10exp * 113; KmNaip = 11; KmKo = 1.5; % Grandi

ICaPmax = q10_CaP^q10exp * 2.0; kCaP = 0.0005;

kNaCa = q10_NCX^q10exp * 0.0088 * cAF_kNaCa; gam = 0.45; dNaCa = 0.0003;

gIf = 1;

% Ca and Na diffusion
DCa   = q10_DCaNa^q10exp * 833;
DCaSR = q10_DCaBmSR^q10exp * 50.7;
DCaBm = q10_DCaBmSR^q10exp * 28.8;
DNa   = q10_DCaNa^q10exp * 0.146;

% SERCA
base_phos = 0.1 * cAF_phos; % Baseline phosphorylation
PLB_SERCA_ratio = cAF_PLB;
SLN_SERCA_ratio = cAF_SLN;
Kmf_PLBKO = 0.15e-3; Kmf_PLB = 0.12e-3; Kmf_SLN = 0.07e-3;
Kmr_PLBKO = 2.5; Kmr_PLB = 0.88; Kmr_SLN = 0.5;
SERCAKmf = Kmf_PLBKO + Kmf_PLB * PLB_SERCA_ratio * (1 - base_phos) + Kmf_SLN * SLN_SERCA_ratio * (1 - base_phos);
SERCAKmr = Kmr_PLBKO - Kmr_PLB * PLB_SERCA_ratio * (1 - base_phos) - Kmr_SLN * SLN_SERCA_ratio * (1 - base_phos);
k4 = q10_SERCA^q10exp * 17; % pump rate
k3 = k4 / SERCAKmr^2;
k1 = 1000^2 * k4;
k2 = k1 * SERCAKmf^2;
cpumps = 30e-3 / Dvcell * cAF_cpumps; % pump concentration in cytosol volume, in cAF SR does not dilate

% RyR
%RyRtauadapt   = 1; % Pit??isk?? t??t?? nopeuttaa ???
RyRtauadapt   = (1/q10_RyR^q10exp) * 0.85;
RyRtauactss   = (1/q10_RyR^q10exp) * 4.3e-3;
RyRtauinactss = (1/q10_RyR^q10exp) * 13e-3;
RyRtauact     = (1/q10_RyR^q10exp) * 16e-3;
RyRtauinact   = (1/q10_RyR^q10exp) * 74e-3;

% SR Ca leak
kSRleak = 6e-3; 


%**********************************************
%% Analytical equations
%**********************************************

% Reversal potentials
ENa = R*T/F * log ( Nao / Nass );
EK = R*T/F * log ( Ko / Ki );
ECa = R*T/F/2 * log ( Cao / y(i_Cass) );

% INa
INa_model = 'new_model';
%INa_model = 'LuoRudy_INa';
%INa_model = 'Nygren_INa';
%INa_model = 'TenTusscher_INa';

switch INa_model

    case 'new_model'
        % threshold stimulus amplitude = 460
        gNa = 0.9*620 * cAF_gNa; % Vmax/s ~235 V/s, accounting for the 10% loss of INa in enzymatic cell isolation
        INaF = gNa * y(i_INam)^3 * y(i_INah1) * y(i_INah2) * (y(i_V) - ENa);
        INaminf = 1 / (1 + exp((y(i_V) + 39 + INa_V_shift_act)/-7.2));
        INahinf = 1 / (1 + exp((y(i_V) + 67 + INa_V_shift_inact)/6.0));

        INah1inf = INahinf;
        INah2inf = INahinf;
        
        INamtau = (1/q10_INa^q10exp) * (0.00001 + 0.00013*exp(-((y(i_V) + 48 + INa_V_shift_act)/15)^2) + 0.000045 / (1 + exp((y(i_V) + 42 + INa_V_shift_act)/-5)));
        INah1tau = (1/q10_INa^q10exp) * (0.00007 + 0.034 / (1 + exp((y(i_V) + 41 + INa_V_shift_inact)/5.5) + exp(-(y(i_V) + 41 + INa_V_shift_inact)/14)) + 0.0002 / (1 + exp(-(y(i_V) + 79 + INa_V_shift_act)/14))); % q10 = 2.8
        INah2tau = (1/q10_INa^q10exp) * (0.0007 + 0.15 / (1 + exp((y(i_V) + 41 + INa_V_shift_inact)/5.5) + exp(-(y(i_V) + 41 + INa_V_shift_inact)/14)) + 0.002 / (1 + exp(-(y(i_V) + 79 + INa_V_shift_act)/14))); % q10 = 2.8

    case 'LuoRudy_INa'
        % threshold stimulus amplitude = 950
        gNa = 7.8*50 * cAF_gNa; % Courtemanche: 7.8 nS/pF
        INaF = gNa * y(i_INam)^3 * y(i_INah1) * y(i_INah2) * (y(i_V) - ENa);
        alpha_m = 0.32 * (y(i_V) + 47.13) ./ (1.0 - exp(-0.1 * (y(i_V) + 47.13)));
        beta_m = 0.8 * exp(-y(i_V)/11);
        INamtau = 0.001 ./ (alpha_m + beta_m); % in seconds
        INaminf = alpha_m ./ (alpha_m + beta_m);

        if y(i_V) >= -40
            alpha_h = 0.0;
            beta_h = 1 ./ (0.13*(1.0+exp((y(i_V)+10.66)/-11.1)));
            INah1tau = 0.001 ./ (alpha_h+beta_h); % in seconds
            INah1inf = alpha_h ./ (alpha_h + beta_h);

            alpha_j = 0.0;
            beta_j = 0.3*exp(-2.535e-7*y(i_V)) ./ (1.0+exp(-0.1*(y(i_V)+32)));
            INah2tau = 0.001 ./ (alpha_j+beta_j); % in seconds
            INah2inf = alpha_j ./ (alpha_j + beta_j);
        else
            alpha_rec_h = 0.135*exp(-(y(i_V)+80.0)/6.8);
            beta_rec_h = 3.56*exp(0.079*y(i_V)) + 3.1e5*exp(0.35*y(i_V));
            INah1tau = 0.001 ./ (alpha_rec_h + beta_rec_h); % in seconds
            INah1inf = alpha_rec_h ./ (alpha_rec_h + beta_rec_h);

            alpha_rec_j = (-1.2714e5*exp(0.2444*y(i_V)) - 3.474e-5 .* exp(-0.04391*y(i_V))) .* (y(i_V)+37.78)/1.0 ./ (1.0+exp(0.311*(y(i_V)+79.23)));
            beta_rec_j = 0.1212*exp(-0.01052*y(i_V)) ./ (1.0+exp(-0.1378*(y(i_V)+40.14)));
            INah2tau = 0.001 ./ (alpha_rec_j + beta_rec_j); % in seconds
            INah2inf = alpha_rec_j ./ (alpha_rec_j + beta_rec_j);
        end

    case 'Nygren_INa'
        % threshold stimulus amplitude = 650
        PNa = 0.0016; % Nygren : 0.0016nL/s
        INaF = PNa .* y(i_INam).^3 .* ( 0.9.*y(i_INah1) + 0.1.*y(i_INah2) ) * Nao * y(i_V) * F^2/(R*T) * ( exp( (y(i_V)-ENa)*F/R/T ) - 1) / ( exp( y(i_V)*F/R/T ) - 1);
        INaminf = 1/(1+exp((y(i_V)+27.12)/-8.21));
        INahinf = 1/(1+exp((y(i_V)+63.6)/5.3));
        INah1inf = INahinf;
        INah2inf = INahinf;
        INamtau = 0.000042*exp( -((y(i_V)+25.57)/28.8).^2 ) + 0.000024;
        INah1tau = 0.03/(1+exp((y(i_V)+35.1)/3.2)) + 0.0003;
        INah2tau = 0.12/(1+exp((y(i_V)+35.1)/3.2)) + 0.003;

    case 'TenTusscher_INa'
        % threshold stimulus amplitude = 220
        gNa = 23*50 * cAF_gNa; % Grandi: 23 mS/?F = 23 nS/pF
        INaF = gNa * y(i_INam)^3 * y(i_INah1) * y(i_INah2) * (y(i_V) - ENa);
        INaminf = 1.0 ./ (1.0+exp((-56.86-y(i_V))/9.03)).^2.0;
        INahinf = 1.0 ./ (1.0+exp((y(i_V)+71.55)/7.43)) .^2.0;
        INah1inf = INahinf;
        INah2inf = INahinf;
        alpha_m = 1.0 ./ (1.0+exp((-60.0-y(i_V))/5.0));
        beta_m = 0.1 ./ (1.0+exp((y(i_V)+35.0)/5.0))+0.1 ./ (1.0+exp((y(i_V)-50.0)/200.0));
        INamtau = 0.001*alpha_m .* beta_m; % in seconds
        if y(i_V) >= -40
            alpha_h = 0.0;
            beta_h = 0.77 ./ (0.13*(1.0+exp((y(i_V)+10.66)/-11.1)));
            INah1tau = 0.001 ./ (alpha_h+beta_h); % in seconds
            alpha_j = 0.0;
            beta_j = 0.6*exp(0.057*y(i_V)) ./ (1.0+exp(-0.1*(y(i_V)+32.0)));
            INah2tau = 0.001 ./ (alpha_j+beta_j); % in seconds

        else
            alpha_h = 0.057*exp(-(y(i_V)+80.0)/6.8);
            beta_h = 2.7*exp(0.079*y(i_V))+310000.0*exp(0.3485*y(i_V));
            INah1tau = 0.001 ./ (alpha_h + beta_h); % in seconds
            alpha_j = (-25428.0*exp(0.2444*y(i_V))-6.948e-6*exp(-0.04391*y(i_V))) .* (y(i_V)+37.78)/1.0 ./ (1.0+exp(0.311*(y(i_V)+79.23)));
            beta_j = 0.02424*exp(-0.01052*y(i_V)) ./ (1.0+exp(-0.1378*(y(i_V)+40.14)));
            INah2tau = 0.001 ./ (alpha_j+beta_j); % in seconds
        end
        
end

INaL_hinf = 1 / (1 + exp((y(i_V) + 72 + INa_V_shift_inact)/5.1));
INaL_tauh = (1/q10_INa^q10exp) *  0.200; % O'Hara et al.
INaL = gNaL * y(i_INam)^3 * y(i_ICaLf1) * (y(i_V) - ENa);

INa = INaF + INaL;

% ICaL
%ICaLfcainf = 1 / ( 1 + (y(i_Cass)/kCa)^2); % PLoS 2014
ICaLfcainf = 1 / ( 1 + (y(i_Cass)/kCa)); % Courtemanche with larger kCa
ICaL = gCaL * y(i_ICaLd) * y(i_ICaLf2) * y(i_ICaLfca) * (y(i_V) - ECa_app); % Nygren without f1
%ICaLdinf = 1/(1+exp((y(i_V)+9)/-5.8)); % Nygren
%ICaLdinf = (1 + exp((-10 - y(i_V))/8))^-1; % Courtemanche
ICaLdinf = 1/(1+exp((y(i_V)+9.5)/-6.9)); % "ave"
ICaLfinf = 0.04 + 0.96 / (1 + exp((y(i_V) + 25.5)/8.4)) + 1 / (1 + exp(-(y(i_V) - 60)/8.0)); % fit Li et al. data
ICaLdtau = 0.00065 * exp(-((y(i_V) + 35)/30).^2) + 0.0005; % Nygren four-fold faster, according Cavalie
ICaLf2tau = 1.34*exp( -((y(i_V)+40)/14.2).^2 ) + 0.04; % average from Li et al. and Christ et al. (55 + 51.8 + 20.6 + 33.9)/4
ICaLfcatau = 2e-3;

% Ito
It = gt * y(i_Itr) * y(i_Its) * (y(i_V) - EK);
Itrinf = 1/(1+exp((y(i_V)-1)/-11)); % Nygren
Itsinf = 1/(1+exp((y(i_V)+40.5)/11.5)); % Nygren
Itrtau = (1/q10_Ito^q10exp) * (0.0024*exp( -((y(i_V) + 0)/30)^2 ) + 0.0010); % Nygren at 37C
%Itstau = 0.025635*exp( -((y(i_V)+52.45)/15.8827).^2 ) + 0.01414; % Maleckar et al.
Itstau = (1/q10_Ito^q10exp) * (0.018*exp( -((y(i_V) + 52.45)/15.88)^2 ) + 0.0096); % Maleckar at 37C

% Isus
Isus = gsus * y(i_Isusr) * y(i_Isuss) * (y(i_V) - EK);
Isusrinf = 1/(1 + exp((y(i_V) + 6)/-8.6)); % Maleckar et al.
Isussinf = 1/(1 + exp((y(i_V) + 7.5)/10)); % Maleckar et al.
Isusrtau = (1/q10_Isus^q10exp) * (0.0066/(1 + exp((y(i_V) + 5)/12)) + 0.00036); % Maleckar at 37C
Isusstau = (1/q10_Isus^q10exp) * (0.43/(1 + exp((y(i_V) + 60)/10)) + 2.2); % Maleckar at 37C

% IKr
IKrpi = 1/(1 + exp((y(i_V) + 74)/24)); % Grandi
IKrpainf = 1/(1 + exp(-(y(i_V) + 10)/5)); % xrss
IKrpatau = 0.550/(1 + exp((-22-y(i_V))/9)) * 6/(1 + exp((y(i_V)-(-11))/9)) + 0.230/(1 + exp((y(i_V)-(-40))/20)); % tauxr
IKr = gKr * sqrt(Ko/5.4) * y(i_IKrpa) * IKrpi * (y(i_V) - EK);

% IKs
IKsninf = 1 / (1+exp(-(y(i_V) + 3.8)/14.25)); % Grandi, xsss
IKsntau = 0.9901 / (1+exp(-(y(i_V) + 2.436)/14.12)); % tauxs
IKs = gKs * y(i_IKsn)^2 * (y(i_V) - EK);

% IK1
IK1 = gK1 * Ko.^0.4457 * (y(i_V) - EK) / (1 + exp(1.5*(y(i_V)-EK+3.6)*F/R/T));

% IKACh
IKACh = gKACh * sqrt(Ko/5.4) * (y(i_V) - EK) * (0.055 + 0.400 / (1 + exp((y(i_V) - EK + 9.53)/17.18)));

% Background leaks
INab = gNab * (y(i_V) - ENa);
ICab = gCab * (y(i_V) - ECa);

% INaK
%INaK = INaKmax * Ko./(Ko + kNaKK) * Nass.^1.5/(Nass.^1.5 + kNaKNa.^1.5) * (y(i_V) + 150) / (y(i_V) + 200); % Nygren
% Grandi
sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(i_V)*F/R/T)+0.0365*sigma*exp(-y(i_V)*F/R/T));
INaK = IbarNaK*fnak*Ko /(1+(KmNaip/Nass)^4) /(Ko+KmKo);

% INaCa
fCaNCX = 1;
INaCa = kNaCa * ( (exp( gam.*y(i_V)*F/R/T ) .* Nass.^3 .* Cao - exp( (gam-1).*y(i_V)*F/R/T ) .* Nao^3 .* y(i_Cass)*fCaNCX ) / ( 1 + dNaCa*(Nao^3 .* y(i_Cass)*fCaNCX + Nass.^3 .* Cao) ) );

% ICaP
ICaP = ICaPmax * y(i_Cass) / (kCaP + y(i_Cass));

% If, Zorn-Pauly LAW fit
Ifyinf = 1 / (1 + exp((y(i_V)+97.82874)/12.48025));
Ifytau = 1 ./ (0.00332.*exp(-y(i_V)./16.54103)+23.71839.*exp(y(i_V)./16.54103));
IfNa = gIf * y(i_Ify)*((0.2677)*(y(i_V)-ENa));
IfK = gIf * y(i_Ify)*((1-0.2677)*(y(i_V)-EK)); 
If = IfK + IfNa;

% IKCa, rectification based on Tao et al. (2015)
IKCa = gKCa * y(i_IKCa_O) * (1 / (1 + exp((y(i_V) - EK + 120)/45))) * (y(i_V) - EK);

% Ca buffers
betass = ( 1 + SLlow*KdSLlow./(y(i_Cass) + KdSLlow).^2 + SLhigh*KdSLhigh./(y(i_Cass) + KdSLhigh).^2 + BCa*KdBCa./(y(i_Cass) + KdBCa).^2  ).^(-1);
betai = ( 1 + BCa.*KdBCa./(Cai + KdBCa).^2  ).^(-1);
gammai = BCa.*KdBCa./(Cai + KdBCa).^2;

betaSR = ( 1 + CSQN.*KdCSQN./(CaSR + KdCSQN).^2 ).^(-1); 

betaNass = ( 1 + BNa*KdBNa./(y(i_Nass) + KdBNa).^2 ).^(-1);

%  Diffusion from junct to non-junct
Jj_nj = DCa * Aj_nj / xj_nj * (y(i_Cass)-Cai(end)).*1e-6;

% SERCA fluxes 
J_SERCASR1 = (-k3*CaSR(1).^2*(cpumps-y(i_SERCACa1))+k4*y(i_SERCACa1))*Vnonjunct(1)*2; % in 1 nl volume
J_bulkSERCA1 = (k1*Cai(1).^2*(cpumps-y(i_SERCACa1))-k2*y(i_SERCACa1))*Vnonjunct(1)*2; % in 1 nl volume

J_SERCASR2 = (-k3*CaSR(2).^2*(cpumps-y(i_SERCACa2))+k4*y(i_SERCACa2))*Vnonjunct(2)*2; % in 1 nl volume
J_bulkSERCA2 = (k1*Cai(2).^2*(cpumps-y(i_SERCACa2))-k2*y(i_SERCACa2))*Vnonjunct(2)*2; % in 1 nl volume

J_SERCASR3 = (-k3*CaSR(3).^2*(cpumps-y(i_SERCACa3))+k4*y(i_SERCACa3))*Vnonjunct(3)*2; % in 1 nl volume
J_bulkSERCA3 = (k1*Cai(3).^2*(cpumps-y(i_SERCACa3))-k2*y(i_SERCACa3))*Vnonjunct(3)*2; % in 1 nl volume

J_SERCASRss = (-k3*CaSR(4).^2*(cpumps-y(i_SERCACass))+k4*y(i_SERCACass))*Vss*2; % in 1 nl volume
J_bulkSERCAss = (k1*y(i_Cass).^2*(cpumps-y(i_SERCACass))-k2*y(i_SERCACass))*Vss*2; % in 1 nl volume

% RyR
nuss = 900*Vss / Dvcell; % in cAF SR does not dilate
RyRSRCass = (1 - 1./(1 +  exp((CaSR(4)-0.3/cAF_RyR)./0.1)));
RyRainfss = 0.505-0.427./(1 + exp((y(i_Cass).*1000-0.29)./0.082));
RyRoinfss = (1 - 1./(1 +  exp((y(i_Cass).*1000-(y(i_RyRass) + 0.22/cAF_RyR))./0.03)));
RyRcinfss = (1./(1 + exp((y(i_Cass).*1000-(y(i_RyRass)+0.02))./0.01)));
Jrelss = nuss * ( y(i_RyRoss) ) * y(i_RyRcss) * RyRSRCass * ( CaSR(4) -  y(i_Cass) ); 

nu1 = 1.6*Vnonjunct(1) / Dvcell; % in cAF SR does not dilate
RyRSRCa1 = (1 - 1./(1 +  exp((CaSR(1)-0.3/cAF_RyR)./0.1)));
RyRainf1 = 0.505-0.427./(1 + exp((Cai(1).*1000-0.29)./0.082));
RyRoinf1 = (1 - 1./(1 +  exp(( Cai(1).*1000-(y(i_RyRa1) + 0.22/cAF_RyR))./0.03)));
RyRcinf1 = (1./(1 +  exp(( Cai(1).*1000-(y(i_RyRa1)+0.02))./0.01)));
Jrel1 = nu1 * ( y(i_RyRo1) ) * y(i_RyRc1) * RyRSRCa1 * ( CaSR(1) -  Cai(1) ); 

nu2 = 1.6*Vnonjunct(2) / Dvcell; % in cAF SR does not dilate
RyRSRCa2 = (1 - 1./(1 +  exp((CaSR(2)-0.3/cAF_RyR)./0.1)));
RyRainf2 =  0.505-0.427./(1 + exp((Cai(2).*1000-0.29)./0.082));
RyRoinf2 = (1 - 1./(1 +  exp(( Cai(2).*1000-(y(i_RyRa2) + 0.22/cAF_RyR))./0.03)));
RyRcinf2 = (1./(1 +  exp(( Cai(2).*1000-(y(i_RyRa2)+0.02))./0.01)));
Jrel2 = nu2 * ( y(i_RyRo2) ) * y(i_RyRc2) * RyRSRCa2 * ( CaSR(2) -  Cai(2) ); 

nu3 = 1.6*Vnonjunct(3) / Dvcell; % in cAF SR does not dilate
RyRSRCa3 = (1 - 1./(1 +  exp((CaSR(3)-0.3/cAF_RyR)./0.1)));
RyRainf3 =  0.505-0.427./(1 + exp((Cai(3).*1000-0.29)./0.082));
RyRoinf3 = (1 - 1./(1 +  exp(( Cai(3).*1000-(y(i_RyRa3) + 0.22/cAF_RyR))./0.03)));
RyRcinf3 = (1./(1 +  exp(( Cai(3).*1000-(y(i_RyRa3)+0.02))./0.01)));
Jrel3 = nu3 * ( y(i_RyRo3) ) * y(i_RyRc3) * RyRSRCa3 * ( CaSR(3) -  Cai(3) ); 

% SR leak fluxes
JSRCaleak1 = kSRleak * ( CaSR(1) - Cai(1) ) * Vnonjunct(1) / Dvcell; % in cAF SR does not dilate
JSRCaleak2 = kSRleak * ( CaSR(2) - Cai(2) ) * Vnonjunct(2) / Dvcell; % in cAF SR does not dilate
JSRCaleak3 = kSRleak * ( CaSR(3) - Cai(3) ) * Vnonjunct(3) / Dvcell; % in cAF SR does not dilate
JSRCaleakss = kSRleak * ( CaSR(4) - y(i_Cass) ) * Vss / Dvcell; % in cAF SR does not dilate

% Cafluxes in 1 nl volume
JCa = zeros(length(j),1);
JCa(1) = -J_bulkSERCA1 + JSRCaleak1 + Jrel1;
JCa(2) = -J_bulkSERCA2 + JSRCaleak2 + Jrel2;
JCa(3) = -J_bulkSERCA3 + JSRCaleak3 + Jrel3;
JCa(4) = Jj_nj;
JCass = -Jj_nj + JSRCaleakss - J_bulkSERCAss + Jrelss;

JSRCa = zeros(length(j),1);
JSRCa(1) = J_SERCASR1 - JSRCaleak1 - Jrel1;
JSRCa(2) = J_SERCASR2 - JSRCaleak2 - Jrel2;
JSRCa(3) = J_SERCASR3 - JSRCaleak3 - Jrel3;
JSRCa(4) = J_SERCASRss - JSRCaleakss - Jrelss;

% Naflux in 1 nl volume
JNa = DNa * Aj_nj / xj_nj_Nai * (y(i_Nass) - y(i_Nai))* 1e-6;


%**********************************************
% Differential equations
%**********************************************

dy = zeros(length(j)*2+i_Cacenter-1,1);

% V
dy(i_V) = (INa + ICaL + It + Isus + IK1 + IKACh + IKr + IKs + INab + ICab + INaK + ICaP + INaCa + If + IKCa + Istim)/(-Cm); % currents are in (pA)

% INa
dy(i_INam) = (INaminf - y(i_INam))/INamtau;
dy(i_INah1) = (INah1inf - y(i_INah1))/INah1tau;
dy(i_INah2) = (INah2inf - y(i_INah2))/INah2tau;

% ICaL
dy(i_ICaLd) = (ICaLdinf - y(i_ICaLd))/ICaLdtau;
%dy(i_ICaLf1) = (ICaLfinf - y(i_ICaLf1))/ICaLf1tau;
dy(i_ICaLf1) = (INaL_hinf - y(i_ICaLf1))/INaL_tauh; % for INaL
dy(i_ICaLf2) = (ICaLfinf - y(i_ICaLf2))/ICaLf2tau;
dy(i_ICaLfca) = (ICaLfcainf - y(i_ICaLfca))/ICaLfcatau;

% It
dy(i_Itr) = (Itrinf - y(i_Itr))/Itrtau;
dy(i_Its) = (Itsinf - y(i_Its))/Itstau;

% Isus
dy(i_Isusr) = (Isusrinf - y(i_Isusr))/Isusrtau;
dy(i_Isuss) = (Isussinf - y(i_Isuss))/Isusstau;

% IKs
dy(i_IKsn) = (IKsninf - y(i_IKsn))/IKsntau;

% IKr
dy(i_IKrpa) = (IKrpainf - y(i_IKrpa))/IKrpatau;

% If
dy(i_Ify) = (Ifyinf - y(i_Ify))/Ifytau;

% IKCa, rates from Hirschberg et al. (1998)
dy(i_IKCa_O) = (1 - y(i_IKCa_O)) * 47e6 * y(i_Cass)^2 - y(i_IKCa_O) * 13;

% SERCACa
dy(i_SERCACa1) = 0.5*(-J_SERCASR1 + J_bulkSERCA1)/Vnonjunct(1); 
dy(i_SERCACa2) = 0.5*(-J_SERCASR2 + J_bulkSERCA2)/Vnonjunct(2); 
dy(i_SERCACa3) = 0.5*(-J_SERCASR3 + J_bulkSERCA3)/Vnonjunct(3); 
dy(i_SERCACass) = 0.5*(-J_SERCASRss + J_bulkSERCAss)/Vss; 

% RyR
dy(i_RyRoss) = (RyRoinfss-y(i_RyRoss))./RyRtauactss;
dy(i_RyRcss) = (RyRcinfss-y(i_RyRcss))./RyRtauinactss;
dy(i_RyRass) = (RyRainfss-y(i_RyRass))./RyRtauadapt;
dy(i_RyRo1) = (RyRoinf1-y(i_RyRo1))./RyRtauact;
dy(i_RyRc1) = (RyRcinf1-y(i_RyRc1))./RyRtauinact;
dy(i_RyRa1) = (RyRainf1-y(i_RyRa1))./RyRtauadapt;
dy(i_RyRo2) = (RyRoinf2-y(i_RyRo2))./RyRtauact;
dy(i_RyRc2) = (RyRcinf2-y(i_RyRc2))./RyRtauinact;
dy(i_RyRa2) = (RyRainf2-y(i_RyRa2))./RyRtauadapt;
dy(i_RyRo3) = (RyRoinf3-y(i_RyRo3))./RyRtauact;
dy(i_RyRc3) = (RyRcinf3-y(i_RyRc3))./RyRtauinact;
dy(i_RyRa3) = (RyRainf3-y(i_RyRa3))./RyRtauadapt;

% Nai & Ki
dy(i_Nass) = betaNass * (-JNa/Vss -(INa + INab + 3*INaK + 3*INaCa + IfNa) / (Vss*F));
dy(i_Nai) = JNa/Vnonjunct_Nai;
dy(i_Ki) = -(It + Isus + IK1 + IKACh + IKr + IKs - 2*INaK + IfK + IKCa + Istim) / (Vcytosol*F);

% Ca
dy(i_Cass) = betass * ( JCass/Vss + (-ICaL - ICab - ICaP + 2*INaCa) / (2*Vss*F) );

dCaidt = zeros(length(Cai),1); 
dCaidt(1) = betai(1) .* (DCa + gammai(1).*DCaBm) .* ( (Cai(2)-2.*Cai(1)+Cai(1))./dx.^2 + (Cai(2)-Cai(1))./(2.*j(1).*dx.^2) ) - 2.*betai(1).*gammai(1).*DCaBm./(KdBCa + Cai(1)) .* ((Cai(2)-Cai(1))./(2.*dx)).^2 + JCa(1)./Vnonjunct(1).*betai(1); 
dCaidt(2:end-1) = betai(2:end-1) .* (DCa + gammai(2:end-1).*DCaBm) .* ( (Cai(3:end)-2.*Cai(2:end-1)+Cai(1:end-2))./dx.^2 + (Cai(3:end)-Cai(1:end-2))./(2.*j(2:end-1).*dx.^2) ) - 2.*betai(2:end-1).*gammai(2:end-1).*DCaBm./(KdBCa + Cai(2:end-1)) .* ((Cai(3:end)-Cai(1:end-2))./(2.*dx)).^2 + JCa(2:end-1)./Vnonjunct(2:end-1).*betai(2:end-1); 
dCaidt(end) = betai(end) .* (DCa + gammai(end).*DCaBm) .* ( (Cai(end)-2.*Cai(end)+Cai(end-1))./dx.^2 + (Cai(end)-Cai(end-1))./(2.*j(end).*dx.^2) ) - 2.*betai(end).*gammai(end).*DCaBm./(KdBCa + Cai(end)) .* ((Cai(end)-Cai(end-1))./(2.*dx)).^2 + JCa(end)./Vnonjunct(end).*betai(end); 

dy(i_Cacenter:length(Cai)+i_Cacenter-1) = dCaidt;

dCaSRdt = zeros(length(CaSR),1); 
dCaSRdt(1) = betaSR(1) .* (DCaSR) .* ( (CaSR(2)-2.*CaSR(1)+CaSR(1))./dx.^2 + (CaSR(2)-CaSR(1))./(2.*j(1).*dx.^2) ) + JSRCa(1)./VSR(1).*betaSR(1); 
dCaSRdt(2:end-1) = betaSR(2:end-1) .* (DCaSR) .* ( (CaSR(3:end)-2.*CaSR(2:end-1)+CaSR(1:end-2))./dx.^2 + (CaSR(3:end)-CaSR(1:end-2))./(2.*j(2:end-1).*dx.^2) ) + JSRCa(2:end-1)./VSR(2:end-1).*betaSR(2:end-1); 
dCaSRdt(end) = betaSR(end) .* (DCaSR) .* ( (CaSR(end)-2.*CaSR(end)+CaSR(end-1))./dx.^2 + (CaSR(end)-CaSR(end-1))./(2.*j(end).*dx.^2) ) + JSRCa(end)./VSR(end).*betaSR(end); 

dy(i_Cacenter+length(j):length(j).*2+i_Cacenter-1) = dCaSRdt;
