function [states, varargout] = circ_no_atria_init_states()
  % % Default state values for ODE model: circ_no_atria
  % % -------------------------------------------------
  % %
  % % states = circ_no_atria_init_states();
  % % [states, states_names] = circ_no_atria_init_states();

  % --- Default initial state values --- 
  states = zeros(6, 1);

  % --- Ventricular volumes ---
  states(1) = 332.42; % vlv;
  states(2) = 226.5; % vrv;

  % --- Systemic volumes ---
  states(3) = 1656; % v1s;
  states(4) = 333; % v2s;

  % --- Pulmonary volumes ---
  states(5) = 243; % v1p;
  states(6) = 208; % v2p;

  if nargout == 2

    % --- State names --- 
    state_names = cell(6, 1);

    % --- Ventricular volumes ---
    state_names{1} = 'vlv';
    state_names{2} = 'vrv';

    % --- Systemic volumes ---
    state_names{3} = 'v1s';
    state_names{4} = 'v2s';

    % --- Pulmonary volumes ---
    state_names{5} = 'v1p';
    state_names{6} = 'v2p';
    varargout(1) = {state_names};
  end
end