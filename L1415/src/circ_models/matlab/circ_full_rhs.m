function [values] = circ_full_rhs(t, states,  parameters)
  % Compute the right hand side of the circ_full ODE

  % Assign states
  if length(states)~=8
    error('Expected the states array to be of size 8.');
  end
  vlv=states(1); vrv=states(2); vla=states(3); vra=states(4); v1s=states(5);...
    v2s=states(6); v1p=states(7); v2p=states(8);

  % Assign parameters
  if length(parameters)~=38
    error('Expected the parameters array to be of size 38.');
  end
  restVlvd=parameters(1); restVlvs=parameters(2); restVrvd=parameters(3);...
    restVrvs=parameters(4); Emaxlv=parameters(5); Emaxrv=parameters(6);...
    Eminlv=parameters(7); Eminrv=parameters(8); restVlad=parameters(9);...
    restVlas=parameters(10); restVrad=parameters(11);...
    restVras=parameters(12); Emaxla=parameters(13); Emaxra=parameters(14);...
    Eminla=parameters(15); Eminra=parameters(16); R1_s=parameters(17);...
    R2_s=parameters(18); Rao_s=parameters(19); Rmit=parameters(20);...
    c1s=parameters(21); c2s=parameters(22); pext_s=parameters(23);...
    rest_v1s=parameters(24); rest_v2s=parameters(25); R1_p=parameters(26);...
    R2_p=parameters(27); Rao_p=parameters(28); Rtric=parameters(29);...
    c1p=parameters(30); c2p=parameters(31); p_epi=parameters(32);...
    pext_p=parameters(33); rest_v1p=parameters(34); rest_v2p=parameters(35);...
    bcl=parameters(36); t_atria=parameters(37); twitchperiod=parameters(38);

  % Init return args
  values = zeros(8, 1);

  % Expressions for the Time varying elastance component
  t_lv = mod(t, bcl);
  yv = ((t_lv < twitchperiod)*(0.5 - 0.5*cos(2.0*pi*t_lv/twitchperiod)) +...
    ~(t_lv < twitchperiod)*(0));
  Elv = Eminlv + (Emaxlv - Eminlv)*yv;
  restVlv = restVlvs + (1 - yv)*(restVlvd - restVlvs);
  plv = p_epi + (-restVlv + vlv)*Elv;
  Erv = Eminrv + (Emaxrv - Eminrv)*yv;
  restVrv = restVrvs + (1 - yv)*(restVrvd - restVrvs);
  prv = p_epi + (-restVrv + vrv)*Erv;
  t_la = ((t > t_atria)*(mod(t - t_atria, bcl)) + ~(t > t_atria)*(0));
  ya = ((t_atria < twitchperiod)*(0.5 - 0.5*cos(2.0*pi*t_la/twitchperiod)) +...
    ~(t_atria < twitchperiod)*(0));
  Ela = Eminla + (Emaxla - Eminla)*ya;
  restVla = restVlas + (1 - ya)*(restVlad - restVlas);
  pla = p_epi + (-restVla + vla)*Ela;
  Era = Eminra + (Emaxra - Eminra)*ya;
  restVra = restVras + (1 - ya)*(restVrad - restVras);
  pra = p_epi + (-restVra + vra)*Era;

  % Expressions for the Systemic pressures and flows component
  p1s = ((v1s > 0)*(-pext_s + (-rest_v1s + v1s)/c1s) + ~(v1s > 0)*(-pext_s -...
    rest_v1s/c1s));
  p2s = ((v2s > 0)*(-pext_s + (-rest_v2s + v2s)/c2s) + ~(v2s > 0)*(-pext_s -...
    rest_v2s/c2s));
  qart_s = ((plv > p1s)*((-p1s + plv)/Rao_s) + ~(plv > p1s)*(0));
  q_mit = ((pla > plv)*((-plv + pla)/Rmit) + ~(pla > plv)*(0));
  q2_s = ((p2s > pra)*((-pra + p2s)/R2_s) + ~(p2s > pra)*(0));
  q1_s = (-p2s + p1s)/R1_s;

  % Expressions for the Pulmonary pressures and flows component
  p1p = ((v1p > 0)*(-pext_p + (-rest_v1p + v1p)/c1p) + ~(v1p > 0)*(-pext_p -...
    rest_v1p/c1p));
  p2p = ((v2p > 0)*(-pext_p + (-rest_v2p + v2p)/c2p) + ~(v2p > 0)*(-pext_p -...
    rest_v2p/c2p));
  qart_p = ((prv > p1p)*((-p1p + prv)/Rao_p) + ~(prv > p1p)*(0));
  q_tric = ((pra > prv)*((-prv + pra)/Rtric) + ~(pra > prv)*(0));
  q2_p = ((p2p > pla)*((-pla + p2p)/R2_p) + ~(p2p > pla)*(0));
  q1_p = (-p2p + p1p)/R1_p;

  % Expressions for the Ventricular volumes component
  values(1) = -qart_s + q_mit;
  values(2) = -qart_p + q_tric;

  % Expressions for the Atrial volumes component
  values(3) = -q_mit + q2_p;
  values(4) = -q_tric + q2_s;

  % Expressions for the Systemic volumes component
  values(5) = -q1_s + qart_s;
  values(6) = -q2_s + q1_s;

  % Expressions for the Pulmonary volumes component
  values(7) = -q1_p + qart_p;
  values(8) = -q2_p + q1_p;
end