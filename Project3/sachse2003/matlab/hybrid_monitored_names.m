function [monitored_names] = hybrid_monitored_names()
  % % Monitored value names for ODE model: hybrid
  % % ---------------------- --------------------
  % %
  % % monitored_names = hybrid_monitored_names();

  % --- Monitored names --- 
  monitored_names = cell(27, 1);

  % --- Equation for simulated calcium transient ---
  monitored_names{1} = 'beta';
  monitored_names{2} = 'Cai';

  % --- Crossbridge ---
  monitored_names{3} = 'overlap';
  monitored_names{4} = 'k5';
  monitored_names{5} = 'k7';
  monitored_names{6} = 'velFactor';
  monitored_names{7} = 't1';
  monitored_names{8} = 't2';
  monitored_names{9} = 't3';
  monitored_names{10} = 't4';
  monitored_names{11} = 't5';
  monitored_names{12} = 't6';
  monitored_names{13} = 't7';
  monitored_names{14} = 't8';
  monitored_names{15} = 'active';

  % --- Troponin ---
  monitored_names{16} = 'tCaTCa';
  monitored_names{17} = 'dTCa_dt';

  % --- Tropomyosin ---
  monitored_names{18} = 'dTMon_dt';

  % --- Crossbridge ---
  monitored_names{19} = 'dAMATP_dt';
  monitored_names{20} = 'dMATP_dt';
  monitored_names{21} = 'dMADPP_dt';
  monitored_names{22} = 'dAwMADPP_dt';
  monitored_names{23} = 'dAsMADPP_dt';
  monitored_names{24} = 'dAsMADP_dt';
  monitored_names{25} = 'dAMADP_dt';
  monitored_names{26} = 'dMADP_dt';
  monitored_names{27} = 'dM_dt';
end