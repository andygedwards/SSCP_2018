function [monitored_names] = circ_full_monitored_names()
  % % Monitored value names for ODE model: circ_full
  % % ---------------------- -----------------------
  % %
  % % monitored_names = circ_full_monitored_names();

  % --- Monitored names --- 
  monitored_names = cell(36, 1);

  % --- Time varying elastance ---
  monitored_names{1} = 't_lv';
  monitored_names{2} = 'yv';
  monitored_names{3} = 'Elv';
  monitored_names{4} = 'restVlv';
  monitored_names{5} = 'plv';
  monitored_names{6} = 'Erv';
  monitored_names{7} = 'restVrv';
  monitored_names{8} = 'prv';
  monitored_names{9} = 't_la';
  monitored_names{10} = 'ya';
  monitored_names{11} = 'Ela';
  monitored_names{12} = 'restVla';
  monitored_names{13} = 'pla';
  monitored_names{14} = 'Era';
  monitored_names{15} = 'restVra';
  monitored_names{16} = 'pra';

  % --- Systemic pressures and flows ---
  monitored_names{17} = 'p1s';
  monitored_names{18} = 'p2s';
  monitored_names{19} = 'qart_s';
  monitored_names{20} = 'q_mit';
  monitored_names{21} = 'q2_s';
  monitored_names{22} = 'q1_s';

  % --- Pulmonary pressures and flows ---
  monitored_names{23} = 'p1p';
  monitored_names{24} = 'p2p';
  monitored_names{25} = 'qart_p';
  monitored_names{26} = 'q_tric';
  monitored_names{27} = 'q2_p';
  monitored_names{28} = 'q1_p';

  % --- Ventricular volumes ---
  monitored_names{29} = 'dvlv_dt';
  monitored_names{30} = 'dvrv_dt';

  % --- Atrial volumes ---
  monitored_names{31} = 'dvla_dt';
  monitored_names{32} = 'dvra_dt';

  % --- Systemic volumes ---
  monitored_names{33} = 'dv1s_dt';
  monitored_names{34} = 'dv2s_dt';

  % --- Pulmonary volumes ---
  monitored_names{35} = 'dv1p_dt';
  monitored_names{36} = 'dv2p_dt';
end