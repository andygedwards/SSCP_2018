# Gotran generated code for the  "systemic" model
from __future__ import division

def init_state_values(**values):
    """
    Initialize state values
    """
    # Imports
    import numpy as np
    from utils import Range

    # Init values
    # vlv=229, v1s=1439, v2s=1071
    init_values = np.array([229, 1439, 1071], dtype=np.float_)

    # State indices and limit checker
    state_ind = dict([("vlv",(0, Range())), ("v1s",(1, Range())), ("v2s",(2,\
        Range()))])

    for state_name, value in values.items():
        if state_name not in state_ind:
            raise ValueError("{0} is not a state.".format(state_name))
        ind, range = state_ind[state_name]
        if value not in range:
            raise ValueError("While setting '{0}' {1}".format(state_name,\
                range.format_not_in(value)))

        # Assign value
        init_values[ind] = value

    return init_values

def init_parameter_values(**values):
    """
    Initialize parameter values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Param values
    # restVlvd=28.0, restVlvs=19.0, Emaxlv=0.333, Eminlv=0.0133, R1_s=250,
    # Rao_s=0.01, Rmit=0.01, c1s=40, c2s=400, rest_v1s=0.0,
    # rest_v2s=0.0, bcl=600.0, t_atria=80, twitchperiod=300
    init_values = np.array([28.0, 19.0, 0.333, 0.0133, 250, 0.01, 0.01, 40,\
        400, 0.0, 0.0, 600.0, 80, 300], dtype=np.float_)

    # Parameter indices and limit checker
    param_ind = dict([("restVlvd", (0, Range())), ("restVlvs", (1, Range())),\
        ("Emaxlv", (2, Range())), ("Eminlv", (3, Range())), ("R1_s", (4,\
        Range())), ("Rao_s", (5, Range())), ("Rmit", (6, Range())), ("c1s",\
        (7, Range())), ("c2s", (8, Range())), ("rest_v1s", (9, Range())),\
        ("rest_v2s", (10, Range())), ("bcl", (11, Range())), ("t_atria", (12,\
        Range())), ("twitchperiod", (13, Range()))])

    for param_name, value in values.items():
        if param_name not in param_ind:
            raise ValueError("{0} is not a parameter.".format(param_name))
        ind, range = param_ind[param_name]
        if value not in range:
            raise ValueError("While setting '{0}' {1}".format(param_name,\
                range.format_not_in(value)))

        # Assign value
        init_values[ind] = value

    return init_values

def state_indices(*states):
    """
    State indices
    """
    state_inds = dict([("vlv", 0), ("v1s", 1), ("v2s", 2)])

    indices = []
    for state in states:
        if state not in state_inds:
            raise ValueError("Unknown state: '{0}'".format(state))
        indices.append(state_inds[state])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def parameter_indices(*params):
    """
    Parameter indices
    """
    param_inds = dict([("restVlvd", 0), ("restVlvs", 1), ("Emaxlv", 2),\
        ("Eminlv", 3), ("R1_s", 4), ("Rao_s", 5), ("Rmit", 6), ("c1s", 7),\
        ("c2s", 8), ("rest_v1s", 9), ("rest_v2s", 10), ("bcl", 11),\
        ("t_atria", 12), ("twitchperiod", 13)])

    indices = []
    for param in params:
        if param not in param_inds:
            raise ValueError("Unknown param: '{0}'".format(param))
        indices.append(param_inds[param])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def monitor_indices(*monitored):
    """
    Monitor indices
    """
    monitor_inds = dict([("t_lv", 0), ("yv", 1), ("Elv", 2), ("restVlv", 3),\
        ("plv", 4), ("p1s", 5), ("p2s", 6), ("qart_s", 7), ("q_mit", 8),\
        ("q1_s", 9), ("dvlv_dt", 10), ("dv1s_dt", 11), ("dv2s_dt", 12)])

    indices = []
    for monitor in monitored:
        if monitor not in monitor_inds:
            raise ValueError("Unknown monitored: '{0}'".format(monitor))
        indices.append(monitor_inds[monitor])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def rhs(states, t, parameters, values=None):
    """
    Compute the right hand side of the systemic ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 3)
    vlv, v1s, v2s = states

    # Assign parameters
    assert(len(parameters) == 14)
    restVlvd=parameters[0]; restVlvs=parameters[1]; Emaxlv=parameters[2];\
        Eminlv=parameters[3]; R1_s=parameters[4]; Rao_s=parameters[5];\
        Rmit=parameters[6]; c1s=parameters[7]; c2s=parameters[8];\
        rest_v1s=parameters[9]; rest_v2s=parameters[10]; bcl=parameters[11];\
        twitchperiod=parameters[13]

    # Init return args
    if values is None:
        values = np.zeros((3,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (3,)

    # Expressions for the Time varying elastance component
    t_lv = math.fmod(t, bcl)
    yv = (0.5 - 0.5*math.cos(2.0*math.pi*t_lv/twitchperiod) if t_lv <\
        twitchperiod else 0)
    Elv = Eminlv + (Emaxlv - Eminlv)*yv
    restVlv = restVlvs + (1 - yv)*(restVlvd - restVlvs)
    plv = (-restVlv + vlv)*Elv

    # Expressions for the Systemic pressures and flows component
    p1s = ((-rest_v1s + v1s)/c1s if v1s > 0 else -rest_v1s/c1s)
    p2s = ((-rest_v2s + v2s)/c2s if v2s > 0 else -rest_v2s/c2s)
    qart_s = ((-p1s + plv)/Rao_s if plv > p1s else 0)
    q_mit = ((-plv + p2s)/Rmit if p2s > plv else 0)
    q1_s = (-p2s + p1s)/R1_s

    # Expressions for the Ventricular volume component
    values[0] = -qart_s + q_mit

    # Expressions for the Systemic volumes component
    values[1] = -q1_s + qart_s
    values[2] = -q_mit + q1_s

    # Return results
    return values

def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the systemic ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 3)
    vlv, v1s, v2s = states

    # Assign parameters
    assert(len(parameters) == 14)
    restVlvd=parameters[0]; restVlvs=parameters[1]; Emaxlv=parameters[2];\
        Eminlv=parameters[3]; R1_s=parameters[4]; Rao_s=parameters[5];\
        Rmit=parameters[6]; c1s=parameters[7]; c2s=parameters[8];\
        rest_v1s=parameters[9]; rest_v2s=parameters[10]; bcl=parameters[11];\
        twitchperiod=parameters[13]

    # Init return args
    if monitored is None:
        monitored = np.zeros((13,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (13,)

    # Expressions for the Time varying elastance component
    monitored[0] = math.fmod(t, bcl)
    monitored[1] = (0.5 - 0.5*math.cos(2.0*math.pi*monitored[0]/twitchperiod)\
        if monitored[0] < twitchperiod else 0)
    monitored[2] = Eminlv + (Emaxlv - Eminlv)*monitored[1]
    monitored[3] = restVlvs + (1 - monitored[1])*(restVlvd - restVlvs)
    monitored[4] = (-monitored[3] + vlv)*monitored[2]

    # Expressions for the Systemic pressures and flows component
    monitored[5] = ((-rest_v1s + v1s)/c1s if v1s > 0 else -rest_v1s/c1s)
    monitored[6] = ((-rest_v2s + v2s)/c2s if v2s > 0 else -rest_v2s/c2s)
    monitored[7] = ((-monitored[5] + monitored[4])/Rao_s if monitored[4] >\
        monitored[5] else 0)
    monitored[8] = ((-monitored[4] + monitored[6])/Rmit if monitored[6] >\
        monitored[4] else 0)
    monitored[9] = (-monitored[6] + monitored[5])/R1_s

    # Expressions for the Ventricular volume component
    monitored[10] = -monitored[7] + monitored[8]

    # Expressions for the Systemic volumes component
    monitored[11] = -monitored[9] + monitored[7]
    monitored[12] = -monitored[8] + monitored[9]

    # Return results
    return monitored
