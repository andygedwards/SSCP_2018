function [states, varargout] = hybrid_init_states()
  % % Default state values for ODE model: hybrid
  % % ------------------------------------------
  % %
  % % states = hybrid_init_states();
  % % [states, states_names] = hybrid_init_states();

  % --- Default initial state values --- 
  states = zeros(11, 1);

  % --- Troponin ---
  states(1) = 0.0; % TCa;

  % --- Tropomyosin ---
  states(2) = 0.0; % TMon;

  % --- Crossbridge ---
  states(3) = 0.0; % AMATP;
  states(4) = 0.0385; % MATP;
  states(5) = 0.3846; % MADPP;
  states(6) = 0.5769; % AwMADPP;
  states(7) = 0.0; % AsMADPP;
  states(8) = 0.0; % AsMADP;
  states(9) = 0.0; % AMADP;
  states(10) = 0.0; % MADP;
  states(11) = 0.0; % M;

  if nargout == 2

    % --- State names --- 
    state_names = cell(11, 1);

    % --- Troponin ---
    state_names{1} = 'TCa';

    % --- Tropomyosin ---
    state_names{2} = 'TMon';

    % --- Crossbridge ---
    state_names{3} = 'AMATP';
    state_names{4} = 'MATP';
    state_names{5} = 'MADPP';
    state_names{6} = 'AwMADPP';
    state_names{7} = 'AsMADPP';
    state_names{8} = 'AsMADP';
    state_names{9} = 'AMADP';
    state_names{10} = 'MADP';
    state_names{11} = 'M';
    varargout(1) = {state_names};
  end
end