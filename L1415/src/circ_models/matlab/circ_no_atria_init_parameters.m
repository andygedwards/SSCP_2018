function [parameters, varargout] = circ_no_atria_init_parameters()
  % % Default parameter values for ODE model: circ_no_atria
  % % -----------------------------------------------------
  % %
  % % parameters = circ_no_atria_init_parameters();
  % % [parameters, parameters_names] = circ_no_atria_init_parameter();

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(36, 1);

  % --- Unloaded ventricular volumes ---
  parameters(1) = 28.0; % restVlvd;
  parameters(2) = 19.0; % restVlvs;
  parameters(3) = 23.3; % restVrvd;
  parameters(4) = 17.0; % restVrvs;

  % --- Max and min ventricular elastances ---
  parameters(5) = 0.333; % Emaxlv;
  parameters(6) = 0.0667; % Emaxrv;
  parameters(7) = 0.0133; % Eminlv;
  parameters(8) = 0.0027; % Eminrv;

  % --- Unloaded atrial volumes ---
  parameters(9) = 14.0; % restVlad;
  parameters(10) = 13.0; % restVlas;
  parameters(11) = 14.0; % restVrad;
  parameters(12) = 13.0; % restVras;

  % --- Mmax and min atrial elastances ---
  parameters(13) = 0.07823464; % Emaxla;
  parameters(14) = 0.03000767; % Emaxra;
  parameters(15) = 0.0711224; % Eminla;
  parameters(16) = 0.0272797; % Eminra;

  % --- Systemic circulation ---
  parameters(17) = 246.9382; % R1_s;
  parameters(18) = 0.5; % Rao_s;
  parameters(19) = 0.5; % Rmit;
  parameters(20) = 32.5461712573; % c1s;
  parameters(21) = 432.567206469; % c2s;
  parameters(22) = 0.0; % pext_s;
  parameters(23) = 0.0; % rest_v1s;
  parameters(24) = 0.0; % rest_v2s;

  % --- Pulmonic circulation ---
  parameters(25) = 5.0496; % R1_p;
  parameters(26) = 0.5; % Rao_p;
  parameters(27) = 0.5; % Rtric;
  parameters(28) = 41.7248684095; % c1p;
  parameters(29) = 50.0351139313; % c2p;
  parameters(30) = 0.0; % p_epi;
  parameters(31) = 0.0; % pext_p;
  parameters(32) = 0.0; % rest_v1p;
  parameters(33) = 0.0; % rest_v2p;

  % --- Time ---
  parameters(34) = 600.0; % bcl;
  parameters(35) = 80; % t_atria;
  parameters(36) = 300; % twitchperiod;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(36, 1);

    % --- Unloaded ventricular volumes ---
    parameter_names{1} = 'restVlvd';
    parameter_names{2} = 'restVlvs';
    parameter_names{3} = 'restVrvd';
    parameter_names{4} = 'restVrvs';

    % --- Max and min ventricular elastances ---
    parameter_names{5} = 'Emaxlv';
    parameter_names{6} = 'Emaxrv';
    parameter_names{7} = 'Eminlv';
    parameter_names{8} = 'Eminrv';

    % --- Unloaded atrial volumes ---
    parameter_names{9} = 'restVlad';
    parameter_names{10} = 'restVlas';
    parameter_names{11} = 'restVrad';
    parameter_names{12} = 'restVras';

    % --- Mmax and min atrial elastances ---
    parameter_names{13} = 'Emaxla';
    parameter_names{14} = 'Emaxra';
    parameter_names{15} = 'Eminla';
    parameter_names{16} = 'Eminra';

    % --- Systemic circulation ---
    parameter_names{17} = 'R1_s';
    parameter_names{18} = 'Rao_s';
    parameter_names{19} = 'Rmit';
    parameter_names{20} = 'c1s';
    parameter_names{21} = 'c2s';
    parameter_names{22} = 'pext_s';
    parameter_names{23} = 'rest_v1s';
    parameter_names{24} = 'rest_v2s';

    % --- Pulmonic circulation ---
    parameter_names{25} = 'R1_p';
    parameter_names{26} = 'Rao_p';
    parameter_names{27} = 'Rtric';
    parameter_names{28} = 'c1p';
    parameter_names{29} = 'c2p';
    parameter_names{30} = 'p_epi';
    parameter_names{31} = 'pext_p';
    parameter_names{32} = 'rest_v1p';
    parameter_names{33} = 'rest_v2p';

    % --- Time ---
    parameter_names{34} = 'bcl';
    parameter_names{35} = 't_atria';
    parameter_names{36} = 'twitchperiod';
    varargout(1) = {parameter_names};
  end
end