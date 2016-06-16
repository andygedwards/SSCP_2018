function [states, varargout] = systemic_init_states()
  % % Default state values for ODE model: systemic
  % % --------------------------------------------
  % %
  % % states = systemic_init_states();
  % % [states, states_names] = systemic_init_states();

  % --- Default initial state values --- 
  states = zeros(3, 1);

  % --- Ventricular volume ---
  states(1) = 229; % vlv;

  % --- Systemic volumes ---
  states(2) = 1439; % v1s;
  states(3) = 1071; % v2s;

  if nargout == 2

    % --- State names --- 
    state_names = cell(3, 1);

    % --- Ventricular volume ---
    state_names{1} = 'vlv';

    % --- Systemic volumes ---
    state_names{2} = 'v1s';
    state_names{3} = 'v2s';
    varargout(1) = {state_names};
  end
end