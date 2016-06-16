function [monitored] = systemic_monitor(t, states, parameters)
  % Computes monitored expressions of the systemic ODE

  % Assign states
  if length(states)~=3
    error('Expected the states array to be of size 3.');
  end
  vlv=states(1); v1s=states(2); v2s=states(3);

  % Assign parameters
  if length(parameters)~=14
    error('Expected the parameters array to be of size 14.');
  end
  restVlvd=parameters(1); restVlvs=parameters(2); Emaxlv=parameters(3);...
    Eminlv=parameters(4); R1_s=parameters(5); Rao_s=parameters(6);...
    Rmit=parameters(7); c1s=parameters(8); c2s=parameters(9);...
    rest_v1s=parameters(10); rest_v2s=parameters(11); bcl=parameters(12);...
    twitchperiod=parameters(14);

  % Init return args
  monitored = zeros(13, 1);

  % Expressions for the Time varying elastance component
  monitored(1) = mod(t, bcl);
  monitored(2) = ((monitored(1) < twitchperiod)*(0.5 -...
    0.5*cos(2.0*pi*monitored(1)/twitchperiod)) + ~(monitored(1) <...
    twitchperiod)*(0));
  monitored(3) = Eminlv + (Emaxlv - Eminlv)*monitored(2);
  monitored(4) = restVlvs + (1 - monitored(2))*(restVlvd - restVlvs);
  monitored(5) = (-monitored(4) + vlv)*monitored(3);

  % Expressions for the Systemic pressures and flows component
  monitored(6) = ((v1s > 0)*((-rest_v1s + v1s)/c1s) + ~(v1s >...
    0)*(-rest_v1s/c1s));
  monitored(7) = ((v2s > 0)*((-rest_v2s + v2s)/c2s) + ~(v2s >...
    0)*(-rest_v2s/c2s));
  monitored(8) = ((monitored(5) > monitored(6))*((-monitored(6) +...
    monitored(5))/Rao_s) + ~(monitored(5) > monitored(6))*(0));
  monitored(9) = ((monitored(7) > monitored(5))*((-monitored(5) +...
    monitored(7))/Rmit) + ~(monitored(7) > monitored(5))*(0));
  monitored(10) = (-monitored(7) + monitored(6))/R1_s;

  % Expressions for the Ventricular volume component
  monitored(11) = -monitored(8) + monitored(9);

  % Expressions for the Systemic volumes component
  monitored(12) = -monitored(10) + monitored(8);
  monitored(13) = -monitored(9) + monitored(10);
end