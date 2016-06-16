function [monitored_names] = systemic_monitored_names()
  % % Monitored value names for ODE model: systemic
  % % ---------------------- ----------------------
  % %
  % % monitored_names = systemic_monitored_names();

  % --- Monitored names --- 
  monitored_names = cell(13, 1);

  % --- Time varying elastance ---
  monitored_names{1} = 't_lv';
  monitored_names{2} = 'yv';
  monitored_names{3} = 'Elv';
  monitored_names{4} = 'restVlv';
  monitored_names{5} = 'plv';

  % --- Systemic pressures and flows ---
  monitored_names{6} = 'p1s';
  monitored_names{7} = 'p2s';
  monitored_names{8} = 'qart_s';
  monitored_names{9} = 'q_mit';
  monitored_names{10} = 'q1_s';

  % --- Ventricular volume ---
  monitored_names{11} = 'dvlv_dt';

  % --- Systemic volumes ---
  monitored_names{12} = 'dv1s_dt';
  monitored_names{13} = 'dv2s_dt';
end