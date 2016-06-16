function [monitored_names] = circ_no_atria_monitored_names()
  % % Monitored value names for ODE model: circ_no_atria
  % % ---------------------- ---------------------------
  % %
  % % monitored_names = circ_no_atria_monitored_names();

  % --- Monitored names --- 
  monitored_names = cell(24, 1);

  % --- Time varying elastance ---
  monitored_names{1} = 't_lv';
  monitored_names{2} = 'yv';
  monitored_names{3} = 'Elv';
  monitored_names{4} = 'restVlv';
  monitored_names{5} = 'plv';
  monitored_names{6} = 'Erv';
  monitored_names{7} = 'restVrv';
  monitored_names{8} = 'prv';

  % --- Pressures and flows ---
  monitored_names{9} = 'p1s';
  monitored_names{10} = 'p2s';
  monitored_names{11} = 'p1p';
  monitored_names{12} = 'p2p';
  monitored_names{13} = 'qart_s';
  monitored_names{14} = 'q_mit';
  monitored_names{15} = 'q1_s';
  monitored_names{16} = 'qart_p';
  monitored_names{17} = 'q_tric';
  monitored_names{18} = 'q1_p';

  % --- Ventricular volumes ---
  monitored_names{19} = 'dvlv_dt';
  monitored_names{20} = 'dvrv_dt';

  % --- Systemic volumes ---
  monitored_names{21} = 'dv1s_dt';
  monitored_names{22} = 'dv2s_dt';

  % --- Pulmonary volumes ---
  monitored_names{23} = 'dv1p_dt';
  monitored_names{24} = 'dv2p_dt';
end