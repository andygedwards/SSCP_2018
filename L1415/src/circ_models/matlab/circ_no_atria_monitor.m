function [monitored] = circ_no_atria_monitor(t, states,  parameters)
  % Computes monitored expressions of the circ_no_atria ODE

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
  monitored = zeros(24, 1);

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

  % Expressions for the Pressures and flows component
  monitored(9) = ((v1s > 0)*(-pext_s + (-rest_v1s + v1s)/c1s) + ~(v1s >...
    0)*(-pext_s - rest_v1s/c1s));
  monitored(10) = ((v2s > 0)*(-pext_s + (-rest_v2s + v2s)/c2s) + ~(v2s >...
    0)*(-pext_s - rest_v2s/c2s));
  monitored(11) = ((v1p > 0)*(-pext_p + (-rest_v1p + v1p)/c1p) + ~(v1p >...
    0)*(-pext_p - rest_v1p/c1p));
  monitored(12) = ((v2p > 0)*(-pext_p + (-rest_v2p + v2p)/c2p) + ~(v2p >...
    0)*(-pext_p - rest_v2p/c2p));
  monitored(13) = ((monitored(5) > monitored(9))*((-monitored(9) +...
    monitored(5))/Rao_s) + ~(monitored(5) > monitored(9))*(0));
  monitored(14) = ((monitored(12) > monitored(5))*((-monitored(5) +...
    monitored(12))/Rmit) + ~(monitored(12) > monitored(5))*(0));
  monitored(15) = (-monitored(10) + monitored(9))/R1_s;
  monitored(16) = ((monitored(8) > monitored(11))*((-monitored(11) +...
    monitored(8))/Rao_p) + ~(monitored(8) > monitored(11))*(0));
  monitored(17) = ((monitored(10) > monitored(8))*((-monitored(8) +...
    monitored(10))/Rtric) + ~(monitored(10) > monitored(8))*(0));
  monitored(18) = (-monitored(12) + monitored(11))/R1_p;

  % Expressions for the Ventricular volumes component
  monitored(19) = -monitored(13) + monitored(14);
  monitored(20) = -monitored(16) + monitored(17);

  % Expressions for the Systemic volumes component
  monitored(21) = -monitored(15) + monitored(13);
  monitored(22) = -monitored(17) + monitored(15);

  % Expressions for the Pulmonary volumes component
  monitored(23) = -monitored(18) + monitored(16);
  monitored(24) = -monitored(14) + monitored(18);
end