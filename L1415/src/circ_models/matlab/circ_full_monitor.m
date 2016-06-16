function [monitored] = circ_full_monitor(t, states,  parameters)
  % Computes monitored expressions of the circ_full ODE

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
  monitored = zeros(36, 1);

  % Expressions for the Time varying elastance component
  monitored(1) = mod(t, bcl);
  monitored(2) = ((monitored(1) < twitchperiod)*(0.5 -...
    0.5*cos(2.0*pi*monitored(1)/twitchperiod)) + ~(monitored(1) <...
    twitchperiod)*(0));
  monitored(3) = Eminlv + (Emaxlv - Eminlv)*monitored(2);
  monitored(4) = restVlvs + (1 - monitored(2))*(restVlvd - restVlvs);
  monitored(5) = p_epi + (-monitored(4) + vlv)*monitored(3);
  monitored(6) = Eminrv + (Emaxrv - Eminrv)*monitored(2);
  monitored(7) = restVrvs + (1 - monitored(2))*(restVrvd - restVrvs);
  monitored(8) = p_epi + (-monitored(7) + vrv)*monitored(6);
  monitored(9) = ((t > t_atria)*(mod(t - t_atria, bcl)) + ~(t > t_atria)*(0));
  monitored(10) = ((t_atria < twitchperiod)*(0.5 -...
    0.5*cos(2.0*pi*monitored(9)/twitchperiod)) + ~(t_atria <...
    twitchperiod)*(0));
  monitored(11) = Eminla + (Emaxla - Eminla)*monitored(10);
  monitored(12) = restVlas + (1 - monitored(10))*(restVlad - restVlas);
  monitored(13) = p_epi + (-monitored(12) + vla)*monitored(11);
  monitored(14) = Eminra + (Emaxra - Eminra)*monitored(10);
  monitored(15) = restVras + (1 - monitored(10))*(restVrad - restVras);
  monitored(16) = p_epi + (-monitored(15) + vra)*monitored(14);

  % Expressions for the Systemic pressures and flows component
  monitored(17) = ((v1s > 0)*(-pext_s + (-rest_v1s + v1s)/c1s) + ~(v1s >...
    0)*(-pext_s - rest_v1s/c1s));
  monitored(18) = ((v2s > 0)*(-pext_s + (-rest_v2s + v2s)/c2s) + ~(v2s >...
    0)*(-pext_s - rest_v2s/c2s));
  monitored(19) = ((monitored(5) > monitored(17))*((-monitored(17) +...
    monitored(5))/Rao_s) + ~(monitored(5) > monitored(17))*(0));
  monitored(20) = ((monitored(13) > monitored(5))*((-monitored(5) +...
    monitored(13))/Rmit) + ~(monitored(13) > monitored(5))*(0));
  monitored(21) = ((monitored(18) > monitored(16))*((-monitored(16) +...
    monitored(18))/R2_s) + ~(monitored(18) > monitored(16))*(0));
  monitored(22) = (-monitored(18) + monitored(17))/R1_s;

  % Expressions for the Pulmonary pressures and flows component
  monitored(23) = ((v1p > 0)*(-pext_p + (-rest_v1p + v1p)/c1p) + ~(v1p >...
    0)*(-pext_p - rest_v1p/c1p));
  monitored(24) = ((v2p > 0)*(-pext_p + (-rest_v2p + v2p)/c2p) + ~(v2p >...
    0)*(-pext_p - rest_v2p/c2p));
  monitored(25) = ((monitored(8) > monitored(23))*((-monitored(23) +...
    monitored(8))/Rao_p) + ~(monitored(8) > monitored(23))*(0));
  monitored(26) = ((monitored(16) > monitored(8))*((-monitored(8) +...
    monitored(16))/Rtric) + ~(monitored(16) > monitored(8))*(0));
  monitored(27) = ((monitored(24) > monitored(13))*((-monitored(13) +...
    monitored(24))/R2_p) + ~(monitored(24) > monitored(13))*(0));
  monitored(28) = (-monitored(24) + monitored(23))/R1_p;

  % Expressions for the Ventricular volumes component
  monitored(29) = -monitored(19) + monitored(20);
  monitored(30) = -monitored(25) + monitored(26);

  % Expressions for the Atrial volumes component
  monitored(31) = -monitored(20) + monitored(27);
  monitored(32) = -monitored(26) + monitored(21);

  % Expressions for the Systemic volumes component
  monitored(33) = -monitored(22) + monitored(19);
  monitored(34) = -monitored(21) + monitored(22);

  % Expressions for the Pulmonary volumes component
  monitored(35) = -monitored(28) + monitored(25);
  monitored(36) = -monitored(27) + monitored(28);
end