function [monitored_names] = rice_model_2008_monitored_names()
  % % Monitored value names for ODE model: rice_model_2008
  % % ---------------------- -----------------------------
  % %
  % % monitored_names = rice_model_2008_monitored_names();

  % --- Monitored names --- 
  monitored_names = cell(65, 1);

  % --- Thin filament regulation and crossbridge cycling rates ---
  monitored_names{1} = 'fappT';
  monitored_names{2} = 'gapslmd';
  monitored_names{3} = 'gappT';
  monitored_names{4} = 'hfmd';
  monitored_names{5} = 'hbmd';
  monitored_names{6} = 'hfT';
  monitored_names{7} = 'hbT';
  monitored_names{8} = 'gxbmd';
  monitored_names{9} = 'gxbT';

  % --- Normalised active and passive force ---
  monitored_names{10} = 'SSXBprer';
  monitored_names{11} = 'SSXBpostr';
  monitored_names{12} = 'Fnordv';
  monitored_names{13} = 'force';
  monitored_names{14} = 'active';
  monitored_names{15} = 'ppforce_t';
  monitored_names{16} = 'ppforce_c';
  monitored_names{17} = 'ppforce';
  monitored_names{18} = 'preload';
  monitored_names{19} = 'afterload';
  monitored_names{20} = 'dSL';

  % --- Equation for simulated calcium transient ---
  monitored_names{21} = 'beta';
  monitored_names{22} = 'Cai';

  % --- Ca binding to troponin to thin filament regulation ---
  monitored_names{23} = 'konT';
  monitored_names{24} = 'koffLT';
  monitored_names{25} = 'koffHT';
  monitored_names{26} = 'dTRPNCaL';
  monitored_names{27} = 'dTRPNCaH';
  monitored_names{28} = 'Tropreg';
  monitored_names{29} = 'permtot';
  monitored_names{30} = 'inprmt';
  monitored_names{31} = 'kn_pT';
  monitored_names{32} = 'kp_nT';

  % --- Regulation and crossbridge cycling state equations ---
  monitored_names{33} = 'dXBpostr';
  monitored_names{34} = 'P';
  monitored_names{35} = 'dXBprer';

  % --- Mean strain of strongly bound states ---
  monitored_names{36} = 'dutyprer';
  monitored_names{37} = 'dutypostr';
  monitored_names{38} = 'dxXBprer';
  monitored_names{39} = 'dxXBpostr';

  % --- Calculation of micromolar per millisecondes of Ca for apparent Ca
  % binding ---
  monitored_names{40} = 'FrSBXB';
  monitored_names{41} = 'dFrSBXB';
  monitored_names{42} = 'dsovr_ze';
  monitored_names{43} = 'dsovr_cle';
  monitored_names{44} = 'dlen_sovr';
  monitored_names{45} = 'dSOVFThin';
  monitored_names{46} = 'dSOVFThick';
  monitored_names{47} = 'TropTot';
  monitored_names{48} = 'dTropTot';
  monitored_names{49} = 'dforce';

  % --- Sarcomere geometry ---
  monitored_names{50} = 'sovr_ze';
  monitored_names{51} = 'sovr_cle';
  monitored_names{52} = 'len_sovr';
  monitored_names{53} = 'SOVFThick';
  monitored_names{54} = 'SOVFThin';

  % --- Normalised active and passive force ---
  monitored_names{55} = 'dintf_dt';
  monitored_names{56} = 'dSL_dt';

  % --- Ca binding to troponin to thin filament regulation ---
  monitored_names{57} = 'dTRPNCaL_dt';
  monitored_names{58} = 'dTRPNCaH_dt';

  % --- Regulation and crossbridge cycling state equations ---
  monitored_names{59} = 'dN_NoXB_dt';
  monitored_names{60} = 'dP_NoXB_dt';
  monitored_names{61} = 'dN_dt';
  monitored_names{62} = 'dXBpostr_dt';
  monitored_names{63} = 'dXBprer_dt';

  % --- Mean strain of strongly bound states ---
  monitored_names{64} = 'dxXBprer_dt';
  monitored_names{65} = 'dxXBpostr_dt';
end