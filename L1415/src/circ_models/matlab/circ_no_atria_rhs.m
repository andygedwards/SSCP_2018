function [values] = circ_no_atria_rhs (t, states, parameters)
  % Compute the right hand side of the circ_no_atria ODE

  % Assign states
  if length(states)~=6
    error('Expected the states array to be of size 6.');
  end
  vlv=states(1); vrv=states(2); v1s=states(3); v2s=states(4); v1p=states(5);...
    v2p=states(6);

  % Assign parameters
  if length(parameters)~=36
    error('Expected the parameters array to be of size 36.');
  end
  restVlvd=parameters(1); restVlvs=parameters(2); restVrvd=parameters(3);...
    restVrvs=parameters(4); Emaxlv=parameters(5); Emaxrv=parameters(6);...
    Eminlv=parameters(7); Eminrv=parameters(8); R1_s=parameters(17);...
    Rao_s=parameters(18); Rmit=parameters(19); c1s=parameters(20);...
    c2s=parameters(21); pext_s=parameters(22); rest_v1s=parameters(23);...
    rest_v2s=parameters(24); R1_p=parameters(25); Rao_p=parameters(26);...
    Rtric=parameters(27); c1p=parameters(28); c2p=parameters(29);...
    p_epi=parameters(30); pext_p=parameters(31); rest_v1p=parameters(32);...
    rest_v2p=parameters(33); bcl=parameters(34); twitchperiod=parameters(36);

  % Init return args
  values = zeros(6, 1);

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

  % Expressions for the Pressures and flows component
  p1s = ((v1s > 0)*(-pext_s + (-rest_v1s + v1s)/c1s) + ~(v1s > 0)*(-pext_s -...
    rest_v1s/c1s));
  p2s = ((v2s > 0)*(-pext_s + (-rest_v2s + v2s)/c2s) + ~(v2s > 0)*(-pext_s -...
    rest_v2s/c2s));
  p1p = ((v1p > 0)*(-pext_p + (-rest_v1p + v1p)/c1p) + ~(v1p > 0)*(-pext_p -...
    rest_v1p/c1p));
  p2p = ((v2p > 0)*(-pext_p + (-rest_v2p + v2p)/c2p) + ~(v2p > 0)*(-pext_p -...
    rest_v2p/c2p));
  qart_s = ((plv > p1s)*((-p1s + plv)/Rao_s) + ~(plv > p1s)*(0));
  q_mit = ((p2p > plv)*((-plv + p2p)/Rmit) + ~(p2p > plv)*(0));
  q1_s = (-p2s + p1s)/R1_s;
  qart_p = ((prv > p1p)*((-p1p + prv)/Rao_p) + ~(prv > p1p)*(0));
  q_tric = ((p2s > prv)*((-prv + p2s)/Rtric) + ~(p2s > prv)*(0));
  q1_p = (-p2p + p1p)/R1_p;

  % Expressions for the Ventricular volumes component
  values(1) = -qart_s + q_mit;
  values(2) = -qart_p + q_tric;

  % Expressions for the Systemic volumes component
  values(3) = -q1_s + qart_s;
  values(4) = -q_tric + q1_s;

  % Expressions for the Pulmonary volumes component
  values(5) = -q1_p + qart_p;
  values(6) = -q_mit + q1_p;
end