function [parameters, varargout] = systemic_init_parameters()
  % % Default parameter values for ODE model: systemic
  % % ------------------------------------------------
  % %
  % % parameters = systemic_init_parameters();
  % % [parameters, parameters_names] = systemic_init_parameter();

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(14, 1);

  % --- Unloaded ventricular volumes ---
  parameters(1) = 28.0; % restVlvd;
  parameters(2) = 19.0; % restVlvs;

  % --- Max and min ventricular elastances ---
  parameters(3) = 0.333; % Emaxlv;
  parameters(4) = 0.0133; % Eminlv;

  % --- Systemic circulation ---
  parameters(5) = 250; % R1_s;
  parameters(6) = 0.01; % Rao_s;
  parameters(7) = 0.01; % Rmit;
  parameters(8) = 40; % c1s;
  parameters(9) = 400; % c2s;
  parameters(10) = 0.0; % rest_v1s;
  parameters(11) = 0.0; % rest_v2s;

  % --- Time varying elastance ---
  parameters(12) = 600.0; % bcl;
  parameters(13) = 80; % t_atria;
  parameters(14) = 300; % twitchperiod;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(14, 1);

    % --- Unloaded ventricular volumes ---
    parameter_names{1} = 'restVlvd';
    parameter_names{2} = 'restVlvs';

    % --- Max and min ventricular elastances ---
    parameter_names{3} = 'Emaxlv';
    parameter_names{4} = 'Eminlv';

    % --- Systemic circulation ---
    parameter_names{5} = 'R1_s';
    parameter_names{6} = 'Rao_s';
    parameter_names{7} = 'Rmit';
    parameter_names{8} = 'c1s';
    parameter_names{9} = 'c2s';
    parameter_names{10} = 'rest_v1s';
    parameter_names{11} = 'rest_v2s';

    % --- Time varying elastance ---
    parameter_names{12} = 'bcl';
    parameter_names{13} = 't_atria';
    parameter_names{14} = 'twitchperiod';
    varargout(1) = {parameter_names};
  end
end