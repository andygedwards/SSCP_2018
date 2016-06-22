function [states, varargout] = rice_model_2008_init_states()
  % % Default state values for ODE model: rice_model_2008
  % % ---------------------------------------------------
  % %
  % % states = rice_model_2008_init_states();
  % % [states, states_names] = rice_model_2008_init_states();

  % --- Default initial state values --- 
  states = zeros(11, 1);

  % --- Normalised active and passive force ---
  states(1) = -4.51134525104e-06; % intf;
  states(2) = 1.89999811516; % SL;

  % --- Ca binding to troponin to thin filament regulation ---
  states(3) = 0.0147730085064; % TRPNCaL;
  states(4) = 0.130660965615; % TRPNCaH;

  % --- Regulation and crossbridge cycling state equations ---
  states(5) = 0.999999959256; % N_NoXB;
  states(6) = 4.07437173989e-08; % P_NoXB;
  states(7) = 0.99999783454; % N;
  states(8) = 1.81017564384e-06; % XBpostr;
  states(9) = 3.049496488e-07; % XBprer;

  % --- Mean strain of strongly bound states ---
  states(10) = 3.41212828972e-08; % xXBprer;
  states(11) = 0.00700005394874; % xXBpostr;

  if nargout == 2

    % --- State names --- 
    state_names = cell(11, 1);

    % --- Normalised active and passive force ---
    state_names{1} = 'intf';
    state_names{2} = 'SL';

    % --- Ca binding to troponin to thin filament regulation ---
    state_names{3} = 'TRPNCaL';
    state_names{4} = 'TRPNCaH';

    % --- Regulation and crossbridge cycling state equations ---
    state_names{5} = 'N_NoXB';
    state_names{6} = 'P_NoXB';
    state_names{7} = 'N';
    state_names{8} = 'XBpostr';
    state_names{9} = 'XBprer';

    % --- Mean strain of strongly bound states ---
    state_names{10} = 'xXBprer';
    state_names{11} = 'xXBpostr';
    varargout(1) = {state_names};
  end
end