function [values] = systemic_rhs(t, states, parameters)
  % Compute the right hand side of the systemic ODE

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
  values = zeros(3, 1);

  % Expressions for the Time varying elastance component
  t_lv = mod(t, bcl);
  yv = ((t_lv < twitchperiod)*(0.5 - 0.5*cos(2.0*pi*t_lv/twitchperiod)) +...
    ~(t_lv < twitchperiod)*(0));
  Elv = Eminlv + (Emaxlv - Eminlv)*yv;
  restVlv = restVlvs + (1 - yv)*(restVlvd - restVlvs);
  plv = (-restVlv + vlv)*Elv;

  % Expressions for the Systemic pressures and flows component
  p1s = ((v1s > 0)*((-rest_v1s + v1s)/c1s) + ~(v1s > 0)*(-rest_v1s/c1s));
  p2s = ((v2s > 0)*((-rest_v2s + v2s)/c2s) + ~(v2s > 0)*(-rest_v2s/c2s));
  qart_s = ((plv > p1s)*((-p1s + plv)/Rao_s) + ~(plv > p1s)*(0));
  q_mit = ((p2s > plv)*((-plv + p2s)/Rmit) + ~(p2s > plv)*(0));
  q1_s = (-p2s + p1s)/R1_s;

  % Expressions for the Ventricular volume component
  values(1) = -qart_s + q_mit;

  % Expressions for the Systemic volumes component
  values(2) = -q1_s + qart_s;
  values(3) = -q_mit + q1_s;
end