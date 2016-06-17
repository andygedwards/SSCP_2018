# Gotran generated code for the  "circ_full" model
from __future__ import division

def init_state_values(**values):
    """
    Initialize state values
    """
    # Imports
    import numpy as np
    from utils import Range

    # Init values
    # vlv=93.42, vrv=70.136, vla=39.7717, vra=37.2625, v1s=764, v2s=1860,
    # v1p=112, v2p=101
    init_values = np.array([93.42, 70.136, 39.7717, 37.2625, 764, 1860, 112,\
        101], dtype=np.float_)

    # State indices and limit checker
    state_ind = dict([("vlv",(0, Range())), ("vrv",(1, Range())), ("vla",(2,\
        Range())), ("vra",(3, Range())), ("v1s",(4, Range())), ("v2s",(5,\
        Range())), ("v1p",(6, Range())), ("v2p",(7, Range()))])

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
    state_inds = dict([("vlv", 0), ("vrv", 1), ("vla", 2), ("vra", 3),\
        ("v1s", 4), ("v2s", 5), ("v1p", 6), ("v2p", 7)])

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
    monitor_inds = dict([("t_ventricle", 0), ("t_lv", 1), ("yv", 2), ("Elv",\
        3), ("restVlv", 4), ("plv", 5), ("Erv", 6), ("restVrv", 7), ("prv",\
        8), ("t_la", 9), ("ya", 10), ("Ela", 11), ("restVla", 12), ("pla",\
        13), ("Era", 14), ("restVra", 15), ("pra", 16), ("p1s", 17), ("p2s",\
        18), ("qart_s", 19), ("q_mit", 20), ("q2_s", 21), ("q1_s", 22),\
        ("p1p", 23), ("p2p", 24), ("qart_p", 25), ("q_tric", 26), ("q2_p",\
        27), ("q1_p", 28), ("dvlv_dt", 29), ("dvrv_dt", 30), ("dvla_dt", 31),\
        ("dvra_dt", 32), ("dv1s_dt", 33), ("dv2s_dt", 34), ("dv1p_dt", 35),\
        ("dv2p_dt", 36)])

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
    Compute the right hand side of the circ_full ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 8)
    vlv, vrv, vla, vra, v1s, v2s, v1p, v2p = states

    # Assign parameters
    assert(len(parameters) == 38)
    restVlvd, restVlvs, restVrvd, restVrvs, Emaxlv, Emaxrv, Eminlv, Eminrv,\
        restVlad, restVlas, restVrad, restVras, Emaxla, Emaxra, Eminla,\
        Eminra, R1_s, R2_s, Rao_s, Rmit, c1s, c2s, pext_s, rest_v1s,\
        rest_v2s, R1_p, R2_p, Rao_p, Rtric, c1p, c2p, p_epi, pext_p,\
        rest_v1p, rest_v2p, bcl, t_atria, twitchperiod = parameters

    # Init return args
    if values is None:
        values = np.zeros((8,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (8,)

    # Expressions for the Time varying elastance component
    t_ventricle = 200.0
    t_lv = (math.fmod(t - t_ventricle, bcl) if t > t_ventricle else 0)
    yv = (0.5 - 0.5*math.cos(2.0*math.pi*t_lv/twitchperiod) if t_lv <\
        twitchperiod else 0)
    Elv = Eminlv + (Emaxlv - Eminlv)*yv
    restVlv = restVlvs + (1 - yv)*(restVlvd - restVlvs)
    plv = p_epi + (-restVlv + vlv)*Elv
    Erv = Eminrv + (Emaxrv - Eminrv)*yv
    restVrv = restVrvs + (1 - yv)*(restVrvd - restVrvs)
    prv = p_epi + (-restVrv + vrv)*Erv
    t_la = (math.fmod(t - t_atria, bcl) if t > t_atria else 0)
    ya = (0.5 - 0.5*math.cos(2.0*math.pi*t_la/twitchperiod) if t_atria <\
        twitchperiod else 0)
    Ela = Eminla + (Emaxla - Eminla)*ya
    restVla = restVlas + (1 - ya)*(restVlad - restVlas)
    pla = p_epi + (-restVla + vla)*Ela
    Era = Eminra + (Emaxra - Eminra)*ya
    restVra = restVras + (1 - ya)*(restVrad - restVras)
    pra = p_epi + (-restVra + vra)*Era

    # Expressions for the Systemic pressures and flows component
    p1s = (-pext_s + (-rest_v1s + v1s)/c1s if v1s > 0 else -pext_s -\
        rest_v1s/c1s)
    p2s = (-pext_s + (-rest_v2s + v2s)/c2s if v2s > 0 else -pext_s -\
        rest_v2s/c2s)
    qart_s = ((-p1s + plv)/Rao_s if plv > p1s else 0)
    q_mit = ((-plv + pla)/Rmit if pla > plv else 0)
    q2_s = ((-pra + p2s)/R2_s if p2s > pra else 0)
    q1_s = (-p2s + p1s)/R1_s

    # Expressions for the Pulmonary pressures and flows component
    p1p = (-pext_p + (-rest_v1p + v1p)/c1p if v1p > 0 else -pext_p -\
        rest_v1p/c1p)
    p2p = (-pext_p + (-rest_v2p + v2p)/c2p if v2p > 0 else -pext_p -\
        rest_v2p/c2p)
    qart_p = ((-p1p + prv)/Rao_p if prv > p1p else 0)
    q_tric = ((-prv + pra)/Rtric if pra > prv else 0)
    q2_p = ((-pla + p2p)/R2_p if p2p > pla else 0)
    q1_p = (-p2p + p1p)/R1_p

    # Expressions for the Ventricular volumes component
    values[0] = -qart_s + q_mit
    values[1] = -qart_p + q_tric

    # Expressions for the Atrial volumes component
    values[2] = -q_mit + q2_p
    values[3] = -q_tric + q2_s

    # Expressions for the Systemic volumes component
    values[4] = -q1_s + qart_s
    values[5] = -q2_s + q1_s

    # Expressions for the Pulmonary volumes component
    values[6] = -q1_p + qart_p
    values[7] = -q2_p + q1_p

    # Return results
    return values

def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the circ_full ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 8)
    vlv, vrv, vla, vra, v1s, v2s, v1p, v2p = states

    # Assign parameters
    assert(len(parameters) == 38)
    restVlvd, restVlvs, restVrvd, restVrvs, Emaxlv, Emaxrv, Eminlv, Eminrv,\
        restVlad, restVlas, restVrad, restVras, Emaxla, Emaxra, Eminla,\
        Eminra, R1_s, R2_s, Rao_s, Rmit, c1s, c2s, pext_s, rest_v1s,\
        rest_v2s, R1_p, R2_p, Rao_p, Rtric, c1p, c2p, p_epi, pext_p,\
        rest_v1p, rest_v2p, bcl, t_atria, twitchperiod = parameters

    # Init return args
    if monitored is None:
        monitored = np.zeros((37,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (37,)

    # Expressions for the Time varying elastance component
    monitored[0] = 200.0
    monitored[1] = (math.fmod(t - monitored[0], bcl) if t > monitored[0] else\
        0)
    monitored[2] = (0.5 - 0.5*math.cos(2.0*math.pi*monitored[1]/twitchperiod)\
        if monitored[1] < twitchperiod else 0)
    monitored[3] = Eminlv + (Emaxlv - Eminlv)*monitored[2]
    monitored[4] = restVlvs + (1 - monitored[2])*(restVlvd - restVlvs)
    monitored[5] = p_epi + (-monitored[4] + vlv)*monitored[3]
    monitored[6] = Eminrv + (Emaxrv - Eminrv)*monitored[2]
    monitored[7] = restVrvs + (1 - monitored[2])*(restVrvd - restVrvs)
    monitored[8] = p_epi + (-monitored[7] + vrv)*monitored[6]
    monitored[9] = (math.fmod(t - t_atria, bcl) if t > t_atria else 0)
    monitored[10] = (0.5 -\
        0.5*math.cos(2.0*math.pi*monitored[9]/twitchperiod) if t_atria <\
        twitchperiod else 0)
    monitored[11] = Eminla + (Emaxla - Eminla)*monitored[10]
    monitored[12] = restVlas + (1 - monitored[10])*(restVlad - restVlas)
    monitored[13] = p_epi + (-monitored[12] + vla)*monitored[11]
    monitored[14] = Eminra + (Emaxra - Eminra)*monitored[10]
    monitored[15] = restVras + (1 - monitored[10])*(restVrad - restVras)
    monitored[16] = p_epi + (-monitored[15] + vra)*monitored[14]

    # Expressions for the Systemic pressures and flows component
    monitored[17] = (-pext_s + (-rest_v1s + v1s)/c1s if v1s > 0 else -pext_s\
        - rest_v1s/c1s)
    monitored[18] = (-pext_s + (-rest_v2s + v2s)/c2s if v2s > 0 else -pext_s\
        - rest_v2s/c2s)
    monitored[19] = ((-monitored[17] + monitored[5])/Rao_s if monitored[5] >\
        monitored[17] else 0)
    monitored[20] = ((-monitored[5] + monitored[13])/Rmit if monitored[13] >\
        monitored[5] else 0)
    monitored[21] = ((-monitored[16] + monitored[18])/R2_s if monitored[18] >\
        monitored[16] else 0)
    monitored[22] = (-monitored[18] + monitored[17])/R1_s

    # Expressions for the Pulmonary pressures and flows component
    monitored[23] = (-pext_p + (-rest_v1p + v1p)/c1p if v1p > 0 else -pext_p\
        - rest_v1p/c1p)
    monitored[24] = (-pext_p + (-rest_v2p + v2p)/c2p if v2p > 0 else -pext_p\
        - rest_v2p/c2p)
    monitored[25] = ((-monitored[23] + monitored[8])/Rao_p if monitored[8] >\
        monitored[23] else 0)
    monitored[26] = ((-monitored[8] + monitored[16])/Rtric if monitored[16] >\
        monitored[8] else 0)
    monitored[27] = ((-monitored[13] + monitored[24])/R2_p if monitored[24] >\
        monitored[13] else 0)
    monitored[28] = (-monitored[24] + monitored[23])/R1_p

    # Expressions for the Ventricular volumes component
    monitored[29] = -monitored[19] + monitored[20]
    monitored[30] = -monitored[25] + monitored[26]

    # Expressions for the Atrial volumes component
    monitored[31] = -monitored[20] + monitored[27]
    monitored[32] = -monitored[26] + monitored[21]

    # Expressions for the Systemic volumes component
    monitored[33] = -monitored[22] + monitored[19]
    monitored[34] = -monitored[21] + monitored[22]

    # Expressions for the Pulmonary volumes component
    monitored[35] = -monitored[28] + monitored[25]
    monitored[36] = -monitored[27] + monitored[28]

    # Return results
    return monitored
