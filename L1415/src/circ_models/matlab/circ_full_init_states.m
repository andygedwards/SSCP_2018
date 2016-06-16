function [states, varargout] = circ_full_init_states()
  % % Default state values for ODE model: circ_full
  % % ---------------------------------------------
  % %
  % % states = circ_full_init_states();
  % % [states, states_names] = circ_full_init_states();

  % --- Default initial state values --- 
  states = zeros(8, 1);

  % --- Ventricular volumes ---
  states(1) = 93.42; % vlv;
  states(2) = 70.136; % vrv;

  % --- Atrial volumes ---
  states(3) = 39.7717; % vla;
  states(4) = 37.2625; % vra;

  % --- Systemic volumes ---
  states(5) = 764; % v1s;
  states(6) = 1860; % v2s;

  % --- Pulmonary volumes ---
  states(7) = 112; % v1p;
  states(8) = 101; % v2p;

  if nargout == 2

    % --- State names --- 
    state_names = cell(8, 1);

    % --- Ventricular volumes ---
    state_names{1} = 'vlv';
    state_names{2} = 'vrv';

    % --- Atrial volumes ---
    state_names{3} = 'vla';
    state_names{4} = 'vra';

    % --- Systemic volumes ---
    state_names{5} = 'v1s';
    state_names{6} = 'v2s';

    % --- Pulmonary volumes ---
    state_names{7} = 'v1p';
    state_names{8} = 'v2p';
    varargout(1) = {state_names};
  end
end