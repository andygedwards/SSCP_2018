# Gotran generated code for the  "hybrid" model
from __future__ import division

def init_state_values(**values):
    """
    Initialize state values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Init values
    # TCa=0.0, TMon=0.0, AMATP=0.0, MATP=0.0385, MADPP=0.3846,
    # AwMADPP=0.5769, AsMADPP=0.0, AsMADP=0.0, AMADP=0.0, MADP=0.0,
    # M=0.0
    init_values = np.array([0.0, 0.0, 0.0, 0.0385, 0.3846, 0.5769, 0.0, 0.0,\
        0.0, 0.0, 0.0], dtype=np.float_)

    # State indices and limit checker
    state_ind = dict([("TCa",(0, Range())), ("TMon",(1, Range())),\
        ("AMATP",(2, Range())), ("MATP",(3, Range())), ("MADPP",(4,\
        Range())), ("AwMADPP",(5, Range())), ("AsMADPP",(6, Range())),\
        ("AsMADP",(7, Range())), ("AMADP",(8, Range())), ("MADP",(9,\
        Range())), ("M",(10, Range()))])

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
    # stretch=1.0, velocity=0.0, Ca_amplitude=1.45, Ca_diastolic=0.09,
    # start_time=5, tau1=20, tau2=110, ATP=4.0, F_physiol=1.0,
    # Fmax=0.5, N_v=10.0, TCaMax=50.0, TCaMin=0.0, TCa_stretch=1.0,
    # TMon_coop=2.0, TMon_pow=2.0, detachVel=10.0, k5_stretch=1.0,
    # k5_xb=1.5, k7_base=2.25, k7_force=1.0, k7_stretch=1.25, k_1=1.0,
    # k_10=1.0, k_11=1.0, k_12=0.05, k_13=1.0, k_14=1.0, k_2=1.0,
    # k_3=0.15, k_4=1.5, k_5=0.025, k_6=0.05, k_7=0.03, k_8=0.2,
    # k_9=1.0, k_m1=0.01, k_m3=0.015, k_m4=1.0, k_m5=0.008,
    # k_m6=0.02, k_m8=0.005, k_off=0.04, k_on=0.04, tm_off=0.035,
    # tm_on=0.012, v50=3.0
    init_values = np.array([1.0, 0.0, 1.45, 0.09, 5, 20, 110, 4.0, 1.0, 0.5,\
        10.0, 50.0, 0.0, 1.0, 2.0, 2.0, 10.0, 1.0, 1.5, 2.25, 1.0, 1.25, 1.0,\
        1.0, 1.0, 0.05, 1.0, 1.0, 1.0, 0.15, 1.5, 0.025, 0.05, 0.03, 0.2,\
        1.0, 0.01, 0.015, 1.0, 0.008, 0.02, 0.005, 0.04, 0.04, 0.035, 0.012,\
        3.0], dtype=np.float_)

    # Parameter indices and limit checker
    param_ind = dict([("stretch", (0, Range())), ("velocity", (1, Range())),\
        ("Ca_amplitude", (2, Range())), ("Ca_diastolic", (3, Range())),\
        ("start_time", (4, Range())), ("tau1", (5, Range())), ("tau2", (6,\
        Range())), ("ATP", (7, Range())), ("F_physiol", (8, Range())),\
        ("Fmax", (9, Range())), ("N_v", (10, Range())), ("TCaMax", (11,\
        Range())), ("TCaMin", (12, Range())), ("TCa_stretch", (13, Range())),\
        ("TMon_coop", (14, Range())), ("TMon_pow", (15, Range())),\
        ("detachVel", (16, Range())), ("k5_stretch", (17, Range())),\
        ("k5_xb", (18, Range())), ("k7_base", (19, Range())), ("k7_force",\
        (20, Range())), ("k7_stretch", (21, Range())), ("k_1", (22,\
        Range())), ("k_10", (23, Range())), ("k_11", (24, Range())), ("k_12",\
        (25, Range())), ("k_13", (26, Range())), ("k_14", (27, Range())),\
        ("k_2", (28, Range())), ("k_3", (29, Range())), ("k_4", (30,\
        Range())), ("k_5", (31, Range())), ("k_6", (32, Range())), ("k_7",\
        (33, Range())), ("k_8", (34, Range())), ("k_9", (35, Range())),\
        ("k_m1", (36, Range())), ("k_m3", (37, Range())), ("k_m4", (38,\
        Range())), ("k_m5", (39, Range())), ("k_m6", (40, Range())), ("k_m8",\
        (41, Range())), ("k_off", (42, Range())), ("k_on", (43, Range())),\
        ("tm_off", (44, Range())), ("tm_on", (45, Range())), ("v50", (46,\
        Range()))])

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
    state_inds = dict([("TCa", 0), ("TMon", 1), ("AMATP", 2), ("MATP", 3),\
        ("MADPP", 4), ("AwMADPP", 5), ("AsMADPP", 6), ("AsMADP", 7),\
        ("AMADP", 8), ("MADP", 9), ("M", 10)])

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
    param_inds = dict([("stretch", 0), ("velocity", 1), ("Ca_amplitude", 2),\
        ("Ca_diastolic", 3), ("start_time", 4), ("tau1", 5), ("tau2", 6),\
        ("ATP", 7), ("F_physiol", 8), ("Fmax", 9), ("N_v", 10), ("TCaMax",\
        11), ("TCaMin", 12), ("TCa_stretch", 13), ("TMon_coop", 14),\
        ("TMon_pow", 15), ("detachVel", 16), ("k5_stretch", 17), ("k5_xb",\
        18), ("k7_base", 19), ("k7_force", 20), ("k7_stretch", 21), ("k_1",\
        22), ("k_10", 23), ("k_11", 24), ("k_12", 25), ("k_13", 26), ("k_14",\
        27), ("k_2", 28), ("k_3", 29), ("k_4", 30), ("k_5", 31), ("k_6", 32),\
        ("k_7", 33), ("k_8", 34), ("k_9", 35), ("k_m1", 36), ("k_m3", 37),\
        ("k_m4", 38), ("k_m5", 39), ("k_m6", 40), ("k_m8", 41), ("k_off",\
        42), ("k_on", 43), ("tm_off", 44), ("tm_on", 45), ("v50", 46)])

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
    monitor_inds = dict([("beta", 0), ("Cai", 1), ("overlap", 2), ("k5", 3),\
        ("k7", 4), ("velFactor", 5), ("t1", 6), ("t2", 7), ("t3", 8), ("t4",\
        9), ("t5", 10), ("t6", 11), ("t7", 12), ("t8", 13), ("active", 14),\
        ("tCaTCa", 15), ("dTCa_dt", 16), ("dTMon_dt", 17), ("dAMATP_dt", 18),\
        ("dMATP_dt", 19), ("dMADPP_dt", 20), ("dAwMADPP_dt", 21),\
        ("dAsMADPP_dt", 22), ("dAsMADP_dt", 23), ("dAMADP_dt", 24),\
        ("dMADP_dt", 25), ("dM_dt", 26)])

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
    Compute the right hand side of the hybrid ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 11)
    TCa, TMon, AMATP, MATP, MADPP, AwMADPP, AsMADPP, AsMADP, AMADP, MADP, M =\
        states

    # Assign parameters
    assert(len(parameters) == 47)
    stretch=parameters[0]; velocity=parameters[1];\
        Ca_amplitude=parameters[2]; Ca_diastolic=parameters[3];\
        start_time=parameters[4]; tau1=parameters[5]; tau2=parameters[6];\
        ATP=parameters[7]; Fmax=parameters[9]; N_v=parameters[10];\
        TCa_stretch=parameters[13]; TMon_coop=parameters[14];\
        TMon_pow=parameters[15]; detachVel=parameters[16];\
        k5_stretch=parameters[17]; k5_xb=parameters[18];\
        k7_base=parameters[19]; k7_force=parameters[20];\
        k7_stretch=parameters[21]; k_1=parameters[22]; k_10=parameters[23];\
        k_11=parameters[24]; k_12=parameters[25]; k_13=parameters[26];\
        k_14=parameters[27]; k_2=parameters[28]; k_3=parameters[29];\
        k_4=parameters[30]; k_5=parameters[31]; k_6=parameters[32];\
        k_7=parameters[33]; k_8=parameters[34]; k_9=parameters[35];\
        k_m1=parameters[36]; k_m3=parameters[37]; k_m4=parameters[38];\
        k_m5=parameters[39]; k_m6=parameters[40]; k_m8=parameters[41];\
        k_off=parameters[42]; k_on=parameters[43]; tm_off=parameters[44];\
        tm_on=parameters[45]; v50=parameters[46]

    # Init return args
    if values is None:
        values = np.zeros((11,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (11,)

    # Expressions for the Equation for simulated calcium transient component
    beta = math.pow(tau1/tau2, -1/(-1 + tau1/tau2)) - math.pow(tau1/tau2,\
        -1/(1 - tau2/tau1))
    Cai = (Ca_diastolic + (Ca_amplitude -\
        Ca_diastolic)*(-math.exp((start_time - t)/tau2) +\
        math.exp((start_time - t)/tau1))/beta if t > start_time else\
        Ca_diastolic)

    # Expressions for the Troponin component
    tCaTCa = -k_off*TCa + k_on*math.pow(stretch, TCa_stretch)*(1.0 -\
        TCa)*(2.0 - AMATP - AwMADPP - M - MADP - MADPP - MATP)*Cai
    values[0] = tCaTCa

    # Expressions for the Tropomyosin component
    values[1] = -tm_off*TMon + tm_on*math.pow(1.0 + (TMon_coop +\
        stretch)*TMon, TMon_pow)*(1.0 - TMon)*TCa

    # Expressions for the Crossbridge component
    overlap = (-0.4666667 + 1.4666667*stretch if stretch < 1.0 else (1.0 if\
        stretch <= 1.1 else 2.61333337 - 1.4666667*stretch))
    k5 = k_5*TMon
    k7 = k_7*(k7_base - k7_stretch*stretch + math.fabs(velocity))/(1.0 +\
        k7_force*(1.0 - AMATP - AsMADPP - AwMADPP - M - MADP - MADPP -\
        MATP)*overlap/Fmax)
    velFactor = math.pow(math.fabs(velocity), N_v)/(math.pow(v50, N_v) +\
        math.pow(math.fabs(velocity), N_v))

    # Transition ratios for state variables
    t1 = -k_m1*AMATP + ATP*k_1*(1.0 - AMADP - AMATP - AsMADP - AsMADPP -\
        AwMADPP - M - MADP - MADPP - MATP)
    t2 = k_2*(1.0 + detachVel*velFactor)*AMATP
    t3 = k_3*MATP - k_m3*MADPP
    t4 = k_4*MADPP - k_m4*(1.0 + detachVel*velFactor)*AwMADPP
    t5 = -k_m5*AsMADPP + ((1.0 + k5_xb*(1.0 - AMATP - AwMADPP - M - MADP -\
        MADPP - MATP))*(1.0 + k5_xb*(1.0 - AMATP - AwMADPP - M - MADP - MADPP\
        - MATP)))*(0.4 + k5_stretch*stretch)*AwMADPP*k5
    t6 = k_6*AsMADPP - k_m6*AsMADP
    t7 = AsMADP*k7
    t8 = k_8*AMADP - k_m8*(1.0 - AMADP - AMATP - AsMADP - AsMADPP - AwMADPP -\
        M - MADP - MADPP - MATP)

    # State variables of actin-myosin complex
    values[2] = -t2 + t1
    values[3] = -t3 + ATP*k_14*M + t2
    values[4] = -t4 + k_13*AsMADPP*velFactor + t3
    values[5] = -t5 + t4
    values[6] = -t6 - k_13*AsMADPP*velFactor + t5
    values[7] = -t7 - k_11*AsMADP*velFactor + t6
    values[8] = -t8 - k_10*AMADP*velFactor + t7
    values[9] = -k_12*MADP + k_10*AMADP*velFactor + k_11*AsMADP*velFactor
    values[10] = k_12*MADP + k_9*(1.0 - AMADP - AMATP - AsMADP - AsMADPP -\
        AwMADPP - M - MADP - MADPP - MATP)*velFactor - ATP*k_14*M

    # Return results
    return values

def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the hybrid ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 11)
    TCa, TMon, AMATP, MATP, MADPP, AwMADPP, AsMADPP, AsMADP, AMADP, MADP, M =\
        states

    # Assign parameters
    assert(len(parameters) == 47)
    stretch=parameters[0]; velocity=parameters[1];\
        Ca_amplitude=parameters[2]; Ca_diastolic=parameters[3];\
        start_time=parameters[4]; tau1=parameters[5]; tau2=parameters[6];\
        ATP=parameters[7]; F_physiol=parameters[8]; Fmax=parameters[9];\
        N_v=parameters[10]; TCa_stretch=parameters[13];\
        TMon_coop=parameters[14]; TMon_pow=parameters[15];\
        detachVel=parameters[16]; k5_stretch=parameters[17];\
        k5_xb=parameters[18]; k7_base=parameters[19];\
        k7_force=parameters[20]; k7_stretch=parameters[21];\
        k_1=parameters[22]; k_10=parameters[23]; k_11=parameters[24];\
        k_12=parameters[25]; k_13=parameters[26]; k_14=parameters[27];\
        k_2=parameters[28]; k_3=parameters[29]; k_4=parameters[30];\
        k_5=parameters[31]; k_6=parameters[32]; k_7=parameters[33];\
        k_8=parameters[34]; k_9=parameters[35]; k_m1=parameters[36];\
        k_m3=parameters[37]; k_m4=parameters[38]; k_m5=parameters[39];\
        k_m6=parameters[40]; k_m8=parameters[41]; k_off=parameters[42];\
        k_on=parameters[43]; tm_off=parameters[44]; tm_on=parameters[45];\
        v50=parameters[46]

    # Init return args
    if monitored is None:
        monitored = np.zeros((27,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (27,)

    # Expressions for the Equation for simulated calcium transient component
    monitored[0] = math.pow(tau1/tau2, -1/(-1 + tau1/tau2)) -\
        math.pow(tau1/tau2, -1/(1 - tau2/tau1))
    monitored[1] = (Ca_diastolic + (Ca_amplitude -\
        Ca_diastolic)*(-math.exp((start_time - t)/tau2) +\
        math.exp((start_time - t)/tau1))/monitored[0] if t > start_time else\
        Ca_diastolic)

    # Expressions for the Troponin component
    monitored[15] = -k_off*TCa + k_on*math.pow(stretch, TCa_stretch)*(1.0 -\
        TCa)*(2.0 - AMATP - AwMADPP - M - MADP - MADPP - MATP)*monitored[1]
    monitored[16] = monitored[15]

    # Expressions for the Tropomyosin component
    monitored[17] = -tm_off*TMon + tm_on*math.pow(1.0 + (TMon_coop +\
        stretch)*TMon, TMon_pow)*(1.0 - TMon)*TCa

    # Expressions for the Crossbridge component
    monitored[2] = (-0.4666667 + 1.4666667*stretch if stretch < 1.0 else (1.0 if\
        stretch <= 1.1 else 2.61333337 - 1.4666667*stretch))
    monitored[3] = k_5*TMon
    monitored[4] = k_7*(k7_base - k7_stretch*stretch +\
        math.fabs(velocity))/(1.0 + k7_force*(1.0 - AMATP - AsMADPP - AwMADPP\
        - M - MADP - MADPP - MATP)*monitored[2]/Fmax)
    monitored[5] = math.pow(math.fabs(velocity), N_v)/(math.pow(v50, N_v) +\
        math.pow(math.fabs(velocity), N_v))

    # Transition ratios for state variables
    monitored[6] = -k_m1*AMATP + ATP*k_1*(1.0 - AMADP - AMATP - AsMADP -\
        AsMADPP - AwMADPP - M - MADP - MADPP - MATP)
    monitored[7] = k_2*(1.0 + detachVel*monitored[5])*AMATP
    monitored[8] = k_3*MATP - k_m3*MADPP
    monitored[9] = k_4*MADPP - k_m4*(1.0 + detachVel*monitored[5])*AwMADPP
    monitored[10] = -k_m5*AsMADPP + ((1.0 + k5_xb*(1.0 - AMATP - AwMADPP - M\
        - MADP - MADPP - MATP))*(1.0 + k5_xb*(1.0 - AMATP - AwMADPP - M -\
        MADP - MADPP - MATP)))*(0.4 +\
        k5_stretch*stretch)*AwMADPP*monitored[3]
    monitored[11] = k_6*AsMADPP - k_m6*AsMADP
    monitored[12] = AsMADP*monitored[4]
    monitored[13] = k_8*AMADP - k_m8*(1.0 - AMADP - AMATP - AsMADP - AsMADPP\
        - AwMADPP - M - MADP - MADPP - MATP)

    # State variables of actin-myosin complex
    monitored[18] = -monitored[7] + monitored[6]
    monitored[19] = -monitored[8] + ATP*k_14*M + monitored[7]
    monitored[20] = -monitored[9] + k_13*AsMADPP*monitored[5] + monitored[8]
    monitored[21] = -monitored[10] + monitored[9]
    monitored[22] = -monitored[11] - k_13*AsMADPP*monitored[5] + monitored[10]
    monitored[23] = -monitored[12] - k_11*AsMADP*monitored[5] + monitored[11]
    monitored[24] = -monitored[13] - k_10*AMADP*monitored[5] + monitored[12]
    monitored[25] = -k_12*MADP + k_10*AMADP*monitored[5] +\
        k_11*AsMADP*monitored[5]
    monitored[26] = k_12*MADP + k_9*(1.0 - AMADP - AMATP - AsMADP - AsMADPP -\
        AwMADPP - M - MADP - MADPP - MATP)*monitored[5] - ATP*k_14*M

    # Active force
    monitored[14] = F_physiol*(1.0 - AMATP - AsMADPP - AwMADPP - M - MADP -\
        MADPP - MATP)*monitored[2]/Fmax

    # Return results
    return monitored
