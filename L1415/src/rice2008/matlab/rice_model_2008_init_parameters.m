function [parameters, varargout] = rice_model_2008_init_parameters()
  % % Default parameter values for ODE model: rice_model_2008
  % % -------------------------------------------------------
  % %
  % % parameters = rice_model_2008_init_parameters();
  % % [parameters, parameters_names] = rice_model_2008_init_parameter();

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(56, 1);

  % --- Thin filament regulation and crossbridge cycling rates ---
  parameters(1) = 6.25; % Qfapp;
  parameters(2) = 2.5; % Qgapp;
  parameters(3) = 6.25; % Qgxb;
  parameters(4) = 6.25; % Qhb;
  parameters(5) = 6.25; % Qhf;
  parameters(6) = 0.5; % fapp;
  parameters(7) = 0.07; % gapp;
  parameters(8) = 6; % gslmod;
  parameters(9) = 0.07; % gxb;
  parameters(10) = 0.4; % hb;
  parameters(11) = 0; % hbmdc;
  parameters(12) = 2; % hf;
  parameters(13) = 5; % hfmdc;
  parameters(14) = 1; % sigman;
  parameters(15) = 8; % sigmap;
  parameters(16) = 1; % xbmodsp;

  % --- Normalised active and passive force ---
  parameters(17) = 1; % KSE;
  parameters(18) = 0.02; % PCon_c;
  parameters(19) = 0.002; % PCon_t;
  parameters(20) = 70; % PExp_c;
  parameters(21) = 10; % PExp_t;
  parameters(22) = 1; % SEon;
  parameters(23) = 2.25; % SL_c;
  parameters(24) = 2.4; % SLmax;
  parameters(25) = 1.4; % SLmin;
  parameters(26) = 1.85; % SLrest;
  parameters(27) = 1.9; % SLset;
  parameters(28) = 0; % fixed_afterload;
  parameters(29) = 120; % kxb_normalised;
  parameters(30) = 50; % massf;
  parameters(31) = 3; % visc;

  % --- Equation for simulated calcium transient ---
  parameters(32) = 1.45; % Ca_amplitude;
  parameters(33) = 0.09; % Ca_diastolic;
  parameters(34) = 5; % start_time;
  parameters(35) = 20; % tau1;
  parameters(36) = 110; % tau2;

  % --- Model parameters ---
  parameters(37) = 24; % TmpC;
  parameters(38) = 0.1; % len_hbare;
  parameters(39) = 1.65; % len_thick;
  parameters(40) = 1.2; % len_thin;
  parameters(41) = 0.007; % x_0;

  % --- Ca binding to troponin to thin filament regulation ---
  parameters(42) = 1.6; % Qkn_p;
  parameters(43) = 1.3; % Qkoff;
  parameters(44) = 1.5; % Qkon;
  parameters(45) = 1.6; % Qkp_n;
  parameters(46) = 0.5; % kn_p;
  parameters(47) = 0.025; % koffH;
  parameters(48) = 0.25; % koffL;
  parameters(49) = 1; % koffmod;
  parameters(50) = 0.05; % kon;
  parameters(51) = 0.05; % kp_n;
  parameters(52) = 15; % nperm;
  parameters(53) = 0.5; % perm50;

  % --- Mean strain of strongly bound states ---
  parameters(54) = 2; % xPsi;

  % --- Calculation of micromolar per millisecondes of Ca for apparent Ca
  % binding ---
  parameters(55) = 70; % Trop_conc;
  parameters(56) = 120; % kxb;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(56, 1);

    % --- Thin filament regulation and crossbridge cycling rates ---
    parameter_names{1} = 'Qfapp';
    parameter_names{2} = 'Qgapp';
    parameter_names{3} = 'Qgxb';
    parameter_names{4} = 'Qhb';
    parameter_names{5} = 'Qhf';
    parameter_names{6} = 'fapp';
    parameter_names{7} = 'gapp';
    parameter_names{8} = 'gslmod';
    parameter_names{9} = 'gxb';
    parameter_names{10} = 'hb';
    parameter_names{11} = 'hbmdc';
    parameter_names{12} = 'hf';
    parameter_names{13} = 'hfmdc';
    parameter_names{14} = 'sigman';
    parameter_names{15} = 'sigmap';
    parameter_names{16} = 'xbmodsp';

    % --- Normalised active and passive force ---
    parameter_names{17} = 'KSE';
    parameter_names{18} = 'PCon_c';
    parameter_names{19} = 'PCon_t';
    parameter_names{20} = 'PExp_c';
    parameter_names{21} = 'PExp_t';
    parameter_names{22} = 'SEon';
    parameter_names{23} = 'SL_c';
    parameter_names{24} = 'SLmax';
    parameter_names{25} = 'SLmin';
    parameter_names{26} = 'SLrest';
    parameter_names{27} = 'SLset';
    parameter_names{28} = 'fixed_afterload';
    parameter_names{29} = 'kxb_normalised';
    parameter_names{30} = 'massf';
    parameter_names{31} = 'visc';

    % --- Equation for simulated calcium transient ---
    parameter_names{32} = 'Ca_amplitude';
    parameter_names{33} = 'Ca_diastolic';
    parameter_names{34} = 'start_time';
    parameter_names{35} = 'tau1';
    parameter_names{36} = 'tau2';

    % --- Model parameters ---
    parameter_names{37} = 'TmpC';
    parameter_names{38} = 'len_hbare';
    parameter_names{39} = 'len_thick';
    parameter_names{40} = 'len_thin';
    parameter_names{41} = 'x_0';

    % --- Ca binding to troponin to thin filament regulation ---
    parameter_names{42} = 'Qkn_p';
    parameter_names{43} = 'Qkoff';
    parameter_names{44} = 'Qkon';
    parameter_names{45} = 'Qkp_n';
    parameter_names{46} = 'kn_p';
    parameter_names{47} = 'koffH';
    parameter_names{48} = 'koffL';
    parameter_names{49} = 'koffmod';
    parameter_names{50} = 'kon';
    parameter_names{51} = 'kp_n';
    parameter_names{52} = 'nperm';
    parameter_names{53} = 'perm50';

    % --- Mean strain of strongly bound states ---
    parameter_names{54} = 'xPsi';

    % --- Calculation of micromolar per millisecondes of Ca for apparent Ca
    % binding ---
    parameter_names{55} = 'Trop_conc';
    parameter_names{56} = 'kxb';
    varargout(1) = {parameter_names};
  end
end