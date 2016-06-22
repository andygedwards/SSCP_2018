function [parameters, varargout] = hybrid_init_parameters()
  % % Default parameter values for ODE model: hybrid
  % % ----------------------------------------------
  % %
  % % parameters = hybrid_init_parameters();
  % % [parameters, parameters_names] = hybrid_init_parameter();

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(47, 1);

  % --- hybrid ---
  parameters(1) = 1.0; % stretch;
  parameters(2) = 0.0; % velocity;

  % --- Equation for simulated calcium transient ---
  parameters(3) = 1.45; % Ca_amplitude;
  parameters(4) = 0.09; % Ca_diastolic;
  parameters(5) = 5; % start_time;
  parameters(6) = 20; % tau1;
  parameters(7) = 110; % tau2;

  % --- Crossbridge ---
  parameters(8) = 4.0; % ATP;
  parameters(9) = 1.0; % F_physiol;
  parameters(10) = 0.5; % Fmax;
  parameters(11) = 10.0; % N_v;
  parameters(12) = 50.0; % TCaMax;
  parameters(13) = 0.0; % TCaMin;
  parameters(14) = 1.0; % TCa_stretch;
  parameters(15) = 2.0; % TMon_coop;
  parameters(16) = 2.0; % TMon_pow;
  parameters(17) = 10.0; % detachVel;
  parameters(18) = 1.0; % k5_stretch;
  parameters(19) = 1.5; % k5_xb;
  parameters(20) = 2.25; % k7_base;
  parameters(21) = 1.0; % k7_force;
  parameters(22) = 1.25; % k7_stretch;
  parameters(23) = 1.0; % k_1;
  parameters(24) = 1.0; % k_10;
  parameters(25) = 1.0; % k_11;
  parameters(26) = 0.05; % k_12;
  parameters(27) = 1.0; % k_13;
  parameters(28) = 1.0; % k_14;
  parameters(29) = 1.0; % k_2;
  parameters(30) = 0.15; % k_3;
  parameters(31) = 1.5; % k_4;
  parameters(32) = 0.025; % k_5;
  parameters(33) = 0.05; % k_6;
  parameters(34) = 0.03; % k_7;
  parameters(35) = 0.2; % k_8;
  parameters(36) = 1.0; % k_9;
  parameters(37) = 0.01; % k_m1;
  parameters(38) = 0.015; % k_m3;
  parameters(39) = 1.0; % k_m4;
  parameters(40) = 0.008; % k_m5;
  parameters(41) = 0.02; % k_m6;
  parameters(42) = 0.005; % k_m8;
  parameters(43) = 0.04; % k_off;
  parameters(44) = 0.04; % k_on;
  parameters(45) = 0.035; % tm_off;
  parameters(46) = 0.012; % tm_on;
  parameters(47) = 3.0; % v50;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(47, 1);

    % --- hybrid ---
    parameter_names{1} = 'stretch';
    parameter_names{2} = 'velocity';

    % --- Equation for simulated calcium transient ---
    parameter_names{3} = 'Ca_amplitude';
    parameter_names{4} = 'Ca_diastolic';
    parameter_names{5} = 'start_time';
    parameter_names{6} = 'tau1';
    parameter_names{7} = 'tau2';

    % --- Crossbridge ---
    parameter_names{8} = 'ATP';
    parameter_names{9} = 'F_physiol';
    parameter_names{10} = 'Fmax';
    parameter_names{11} = 'N_v';
    parameter_names{12} = 'TCaMax';
    parameter_names{13} = 'TCaMin';
    parameter_names{14} = 'TCa_stretch';
    parameter_names{15} = 'TMon_coop';
    parameter_names{16} = 'TMon_pow';
    parameter_names{17} = 'detachVel';
    parameter_names{18} = 'k5_stretch';
    parameter_names{19} = 'k5_xb';
    parameter_names{20} = 'k7_base';
    parameter_names{21} = 'k7_force';
    parameter_names{22} = 'k7_stretch';
    parameter_names{23} = 'k_1';
    parameter_names{24} = 'k_10';
    parameter_names{25} = 'k_11';
    parameter_names{26} = 'k_12';
    parameter_names{27} = 'k_13';
    parameter_names{28} = 'k_14';
    parameter_names{29} = 'k_2';
    parameter_names{30} = 'k_3';
    parameter_names{31} = 'k_4';
    parameter_names{32} = 'k_5';
    parameter_names{33} = 'k_6';
    parameter_names{34} = 'k_7';
    parameter_names{35} = 'k_8';
    parameter_names{36} = 'k_9';
    parameter_names{37} = 'k_m1';
    parameter_names{38} = 'k_m3';
    parameter_names{39} = 'k_m4';
    parameter_names{40} = 'k_m5';
    parameter_names{41} = 'k_m6';
    parameter_names{42} = 'k_m8';
    parameter_names{43} = 'k_off';
    parameter_names{44} = 'k_on';
    parameter_names{45} = 'tm_off';
    parameter_names{46} = 'tm_on';
    parameter_names{47} = 'v50';
    varargout(1) = {parameter_names};
  end
end