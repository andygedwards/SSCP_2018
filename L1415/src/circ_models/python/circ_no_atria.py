# Gotran generated code for the  "circ_no_atria" model
from __future__ import division

def init_state_values(**values):
    """
    Initialize state values
    """
    # Imports
    import numpy as np
    from utils import Range

    # Init values
    # vlv=332.42, vrv=226.5, v1s=1656, v2s=333, v1p=243, v2p=208
    init_values = np.array([332.42, 226.5, 1656, 333, 243, 208],\
        dtype=np.float_)

    # State indices and limit checker
    state_ind = dict([("vlv",(0, Range())), ("vrv",(1, Range())), ("v1s",(2,\
        Range())), ("v2s",(3, Range())), ("v1p",(4, Range())), ("v2p",(5,\
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
    from utils import Range

    # Param values
    # restVlvd=28.0, restVlvs=19.0, restVrvd=23.3, restVrvs=17.0,
    # Emaxlv=0.333, Emaxrv=0.0667, Eminlv=0.0133, Eminrv=0.0027,
    # restVlad=14.0, restVlas=13.0, restVrad=14.0, restVras=13.0,
    # Emaxla=0.07823464, Emaxra=0.03000767, Eminla=0.0711224,
    # Eminra=0.0272797, R1_s=246.9382, R2_s=51.4, Rao_s=0.5, Rmit=0.5,
    # c1s=32.5461712573, c2s=432.567206469, pext_s=0.0, rest_v1s=0.0,
    # rest_v2s=0.0, R1_p=5.0496, R2_p=5.0496, Rao_p=0.5, Rtric=0.5,
    # c1p=41.7248684095, c2p=50.0351139313, p_epi=0.0, pext_p=0.0,
    # rest_v1p=0.0, rest_v2p=0.0, bcl=600.0, t_atria=80,
    # twitchperiod=300
    init_values = np.array([28.0, 19.0, 23.3, 17.0, 0.333, 0.0667, 0.0133,\
        0.0027, 14.0, 13.0, 14.0, 13.0, 0.07823464, 0.03000767, 0.0711224,\
        0.0272797, 246.9382, 51.4, 0.5, 0.5, 32.5461712573, 432.567206469,\
        0.0, 0.0, 0.0, 5.0496, 5.0496, 0.5, 0.5, 41.7248684095,\
        50.0351139313, 0.0, 0.0, 0.0, 0.0, 600.0, 80, 300], dtype=np.float_)

    # Parameter indices and limit checker
    param_ind = dict([("restVlvd", (0, Range())), ("restVlvs", (1, Range())),\
        ("restVrvd", (2, Range())), ("restVrvs", (3, Range())), ("Emaxlv",\
        (4, Range())), ("Emaxrv", (5, Range())), ("Eminlv", (6, Range())),\
        ("Eminrv", (7, Range())), ("restVlad", (8, Range())), ("restVlas",\
        (9, Range())), ("restVrad", (10, Range())), ("restVras", (11,\
        Range())), ("Emaxla", (12, Range())), ("Emaxra", (13, Range())),\
        ("Eminla", (14, Range())), ("Eminra", (15, Range())), ("R1_s", (16,\
        Range())), ("R2_s", (17, Range())), ("Rao_s", (18, Range())),\
        ("Rmit", (19, Range())), ("c1s", (20, Range())), ("c2s", (21,\
        Range())), ("pext_s", (22, Range())), ("rest_v1s", (23, Range())),\
        ("rest_v2s", (24, Range())), ("R1_p", (25, Range())), ("R2_p", (26,\
        Range())), ("Rao_p", (27, Range())), ("Rtric", (28, Range())),\
        ("c1p", (29, Range())), ("c2p", (30, Range())), ("p_epi", (31,\
        Range())), ("pext_p", (32, Range())), ("rest_v1p", (33, Range())),\
        ("rest_v2p", (34, Range())), ("bcl", (35, Range())), ("t_atria", (36,\
        Range())), ("twitchperiod", (37, Range()))])

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
    state_inds = dict([("vlv", 0), ("vrv", 1), ("v1s", 2), ("v2s", 3),\
        ("v1p", 4), ("v2p", 5)])

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
    param_inds = dict([("restVlvd", 0), ("restVlvs", 1), ("restVrvd", 2),\
        ("restVrvs", 3), ("Emaxlv", 4), ("Emaxrv", 5), ("Eminlv", 6),\
        ("Eminrv", 7), ("restVlad", 8), ("restVlas", 9), ("restVrad", 10),\
        ("restVras", 11), ("Emaxla", 12), ("Emaxra", 13), ("Eminla", 14),\
        ("Eminra", 15), ("R1_s", 16), ("R2_s", 17), ("Rao_s", 18), ("Rmit",\
        19), ("c1s", 20), ("c2s", 21), ("pext_s", 22), ("rest_v1s", 23),\
        ("rest_v2s", 24), ("R1_p", 25), ("R2_p", 26), ("Rao_p", 27),\
        ("Rtric", 28), ("c1p", 29), ("c2p", 30), ("p_epi", 31), ("pext_p",\
        32), ("rest_v1p", 33), ("rest_v2p", 34), ("bcl", 35), ("t_atria",\
        36), ("twitchperiod", 37)])

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
        ("plv", 4), ("Erv", 5), ("restVrv", 6), ("prv", 7), ("p1s", 8),\
        ("p2s", 9), ("p1p", 10), ("p2p", 11), ("qart_s", 12), ("q_mit", 13),\
        ("q2_s", 14), ("q1_s", 15), ("qart_p", 16), ("q_tric", 17), ("q1_p",\
        18), ("dvlv_dt", 19), ("dvrv_dt", 20), ("dv1s_dt", 21), ("dv2s_dt",\
        22), ("dv1p_dt", 23), ("dv2p_dt", 24)])

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
    Compute the right hand side of the circ_no_atria ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 6)
    vlv, vrv, v1s, v2s, v1p, v2p = states

    # Assign parameters
    assert(len(parameters) == 38)
    restVlvd=parameters[0]; restVlvs=parameters[1]; restVrvd=parameters[2];\
        restVrvs=parameters[3]; Emaxlv=parameters[4]; Emaxrv=parameters[5];\
        Eminlv=parameters[6]; Eminrv=parameters[7]; R1_s=parameters[16];\
        Rao_s=parameters[18]; Rmit=parameters[19]; c1s=parameters[20];\
        c2s=parameters[21]; pext_s=parameters[22]; rest_v1s=parameters[23];\
        rest_v2s=parameters[24]; R1_p=parameters[25]; Rao_p=parameters[27];\
        Rtric=parameters[28]; c1p=parameters[29]; c2p=parameters[30];\
        p_epi=parameters[31]; pext_p=parameters[32]; rest_v1p=parameters[33];\
        rest_v2p=parameters[34]; bcl=parameters[35];\
        twitchperiod=parameters[37]

    # Init return args
    if values is None:
        values = np.zeros((6,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (6,)

    # Expressions for the Time varying elastance component
    t_lv = math.fmod(t, bcl)
    yv = (0.5 - 0.5*math.cos(2.0*math.pi*t_lv/twitchperiod) if t_lv <\
        twitchperiod else 0)
    Elv = Eminlv + (Emaxlv - Eminlv)*yv
    restVlv = restVlvs + (1 - yv)*(restVlvd - restVlvs)
    plv = p_epi + (-restVlv + vlv)*Elv
    Erv = Eminrv + (Emaxrv - Eminrv)*yv
    restVrv = restVrvs + (1 - yv)*(restVrvd - restVrvs)
    prv = p_epi + (-restVrv + vrv)*Erv

    # Expressions for the Pressures and flows component
    p1s = (-pext_s + (-rest_v1s + v1s)/c1s if v1s > 0 else -pext_s -\
        rest_v1s/c1s)
    p2s = (-pext_s + (-rest_v2s + v2s)/c2s if v2s > 0 else -pext_s -\
        rest_v2s/c2s)
    p1p = (-pext_p + (-rest_v1p + v1p)/c1p if v1p > 0 else -pext_p -\
        rest_v1p/c1p)
    p2p = (-pext_p + (-rest_v2p + v2p)/c2p if v2p > 0 else -pext_p -\
        rest_v2p/c2p)
    qart_s = ((-p1s + plv)/Rao_s if plv > p1s else 0)
    q_mit = ((-plv + p2p)/Rmit if p2p > plv else 0)
    q1_s = (-p2s + p1s)/R1_s
    qart_p = ((-p1p + prv)/Rao_p if prv > p1p else 0)
    q_tric = ((-prv + p2s)/Rtric if p2s > prv else 0)
    q1_p = (-p2p + p1p)/R1_p

    # Expressions for the Ventricular volumes component
    values[0] = -qart_s + q_mit
    values[1] = -qart_p + q_tric

    # Expressions for the Systemic volumes component
    values[2] = -q1_s + qart_s
    values[3] = -q_tric + q1_s

    # Expressions for the Pulmonary volumes component
    values[4] = -q1_p + qart_p
    values[5] = -q_mit + q1_p

    # Return results
    return values

def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the circ_no_atria ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 6)
    vlv, vrv, v1s, v2s, v1p, v2p = states

    # Assign parameters
    assert(len(parameters) == 38)
    restVlvd=parameters[0]; restVlvs=parameters[1]; restVrvd=parameters[2];\
        restVrvs=parameters[3]; Emaxlv=parameters[4]; Emaxrv=parameters[5];\
        Eminlv=parameters[6]; Eminrv=parameters[7]; R1_s=parameters[16];\
        R2_s=parameters[17]; Rao_s=parameters[18]; Rmit=parameters[19];\
        c1s=parameters[20]; c2s=parameters[21]; pext_s=parameters[22];\
        rest_v1s=parameters[23]; rest_v2s=parameters[24];\
        R1_p=parameters[25]; Rao_p=parameters[27]; Rtric=parameters[28];\
        c1p=parameters[29]; c2p=parameters[30]; p_epi=parameters[31];\
        pext_p=parameters[32]; rest_v1p=parameters[33];\
        rest_v2p=parameters[34]; bcl=parameters[35];\
        twitchperiod=parameters[37]

    # Init return args
    if monitored is None:
        monitored = np.zeros((25,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (25,)

    # Expressions for the Time varying elastance component
    monitored[0] = math.fmod(t, bcl)
    monitored[1] = (0.5 - 0.5*math.cos(2.0*math.pi*monitored[0]/twitchperiod)\
        if monitored[0] < twitchperiod else 0)
    monitored[2] = Eminlv + (Emaxlv - Eminlv)*monitored[1]
    monitored[3] = restVlvs + (1 - monitored[1])*(restVlvd - restVlvs)
    monitored[4] = p_epi + (-monitored[3] + vlv)*monitored[2]
    monitored[5] = Eminrv + (Emaxrv - Eminrv)*monitored[1]
    monitored[6] = restVrvs + (1 - monitored[1])*(restVrvd - restVrvs)
    monitored[7] = p_epi + (-monitored[6] + vrv)*monitored[5]

    # Expressions for the Pressures and flows component
    monitored[8] = (-pext_s + (-rest_v1s + v1s)/c1s if v1s > 0 else -pext_s -\
        rest_v1s/c1s)
    monitored[9] = (-pext_s + (-rest_v2s + v2s)/c2s if v2s > 0 else -pext_s -\
        rest_v2s/c2s)
    monitored[10] = (-pext_p + (-rest_v1p + v1p)/c1p if v1p > 0 else -pext_p\
        - rest_v1p/c1p)
    monitored[11] = (-pext_p + (-rest_v2p + v2p)/c2p if v2p > 0 else -pext_p\
        - rest_v2p/c2p)
    monitored[12] = ((-monitored[8] + monitored[4])/Rao_s if monitored[4] >\
        monitored[8] else 0)
    monitored[13] = ((-monitored[4] + monitored[11])/Rmit if monitored[11] >\
        monitored[4] else 0)
    monitored[14] = ((-monitored[7] + monitored[9])/R2_s if monitored[9] >\
        monitored[7] else 0)
    monitored[15] = (-monitored[9] + monitored[8])/R1_s
    monitored[16] = ((-monitored[10] + monitored[7])/Rao_p if monitored[7] >\
        monitored[10] else 0)
    monitored[17] = ((-monitored[7] + monitored[9])/Rtric if monitored[9] >\
        monitored[7] else 0)
    monitored[18] = (-monitored[11] + monitored[10])/R1_p

    # Expressions for the Ventricular volumes component
    monitored[19] = -monitored[12] + monitored[13]
    monitored[20] = -monitored[16] + monitored[17]

    # Expressions for the Systemic volumes component
    monitored[21] = -monitored[15] + monitored[12]
    monitored[22] = -monitored[17] + monitored[15]

    # Expressions for the Pulmonary volumes component
    monitored[23] = -monitored[18] + monitored[16]
    monitored[24] = -monitored[13] + monitored[18]

    # Return results
    return monitored
