# Gotran generated code for the  "rice_lmbda" model
from __future__ import division

def init_state_values(**values):
    """
    Initialize state values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Init values
    # lmbda=1.02702702703, TRPNCaH=0.130660965615, TRPNCaL=0.0147730085064,
    # N=0.99999783454, N_NoXB=0.999999959256, P_NoXB=4.07437173989e-08,
    # XBpostr=1.81017564384e-06, XBprer=3.049496488e-07,
    # lmbda_a=1.02702702703
    init_values = np.array([1.02702702703, 0.130660965615, 0.0147730085064,\
        0.99999783454, 0.999999959256, 4.07437173989e-08, 1.81017564384e-06,\
        3.049496488e-07, 1.02702702703], dtype=np.float_)

    # State indices and limit checker
    state_ind = dict([("lmbda",(0, Range())), ("TRPNCaH",(1, Range())),\
        ("TRPNCaL",(2, Range())), ("N",(3, Range())), ("N_NoXB",(4,\
        Range())), ("P_NoXB",(5, Range())), ("XBpostr",(6, Range())),\
        ("XBprer",(7, Range())), ("lmbda_a",(8, Range()))])

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
    # Qfapp=6.25, Qgapp=2.5, Qgxb=6.25, Qhb=6.25, Qhf=6.25, fapp=0.5,
    # gapp=0.07, gslmod=6, gxb=0.07, hb=0.4, hbmdc=0, hf=2,
    # hfmdc=5, sigman=1, sigmap=8, xbmodsp=1, Cp=0.002, KSE=1,
    # SEon=1, SLmax=2.4, SLmin=1.4, SLrest=1.85, SLset=1.9,
    # b_ff=50.0, fixed_afterload=0, kxb=222.0, visc=1.62162162162,
    # Ca_amplitude=1.45, Ca_diastolic=0.09, start_time=5, tau1=20,
    # tau2=110, TmpC=24, eps_0=0.003783783, len_hbare=0.1,
    # len_thick=1.65, len_thin=1.2, Qkn_p=1.6, Qkoff=1.3, Qkon=1.5,
    # Qkp_n=1.6, kn_p=0.5, koffH=0.025, koffL=0.25, koffmod=1,
    # kon=0.05, kp_n=0.05, nperm=15, perm50=0.5, xPsi=4.0
    init_values = np.array([6.25, 2.5, 6.25, 6.25, 6.25, 0.5, 0.07, 6, 0.07,\
        0.4, 0, 2, 5, 1, 8, 1, 0.002, 1, 1, 2.4, 1.4, 1.85, 1.9, 50.0, 0,\
        222.0, 1.62162162162, 1.45, 0.09, 5, 20, 110, 24, 0.003783783, 0.1,\
        1.65, 1.2, 1.6, 1.3, 1.5, 1.6, 0.5, 0.025, 0.25, 1, 0.05, 0.05, 15,\
        0.5, 4.0], dtype=np.float_)

    # Parameter indices and limit checker
    param_ind = dict([("Qfapp", (0, Range())), ("Qgapp", (1, Range())),\
        ("Qgxb", (2, Range())), ("Qhb", (3, Range())), ("Qhf", (4, Range())),\
        ("fapp", (5, Range())), ("gapp", (6, Range())), ("gslmod", (7,\
        Range())), ("gxb", (8, Range())), ("hb", (9, Range())), ("hbmdc",\
        (10, Range())), ("hf", (11, Range())), ("hfmdc", (12, Range())),\
        ("sigman", (13, Range())), ("sigmap", (14, Range())), ("xbmodsp",\
        (15, Range())), ("Cp", (16, Range())), ("KSE", (17, Range())),\
        ("SEon", (18, Range())), ("SLmax", (19, Range())), ("SLmin", (20,\
        Range())), ("SLrest", (21, Range())), ("SLset", (22, Range())),\
        ("b_ff", (23, Range())), ("fixed_afterload", (24, Range())), ("kxb",\
        (25, Range())), ("visc", (26, Range())), ("Ca_amplitude", (27,\
        Range())), ("Ca_diastolic", (28, Range())), ("start_time", (29,\
        Range())), ("tau1", (30, Range())), ("tau2", (31, Range())), ("TmpC",\
        (32, Range())), ("eps_0", (33, Range())), ("len_hbare", (34,\
        Range())), ("len_thick", (35, Range())), ("len_thin", (36, Range())),\
        ("Qkn_p", (37, Range())), ("Qkoff", (38, Range())), ("Qkon", (39,\
        Range())), ("Qkp_n", (40, Range())), ("kn_p", (41, Range())),\
        ("koffH", (42, Range())), ("koffL", (43, Range())), ("koffmod", (44,\
        Range())), ("kon", (45, Range())), ("kp_n", (46, Range())), ("nperm",\
        (47, Range())), ("perm50", (48, Range())), ("xPsi", (49, Range()))])

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
    state_inds = dict([("lmbda", 0), ("TRPNCaH", 1), ("TRPNCaL", 2), ("N",\
        3), ("N_NoXB", 4), ("P_NoXB", 5), ("XBpostr", 6), ("XBprer", 7),\
        ("lmbda_a", 8)])

    indices = []
    for state in states:
        if state not in state_inds:
            raise ValueError("Unknown state: '{0}'".format(state))
        indices.append(state_inds[state])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def parameter_indices(**params):
    """
    Parameter indices
    """
    param_inds = dict(("Qfapp", 0), ("Qgapp", 1), ("Qgxb", 2), ("Qhb", 3),\
        ("Qhf", 4), ("fapp", 5), ("gapp", 6), ("gslmod", 7), ("gxb", 8),\
        ("hb", 9), ("hbmdc", 10), ("hf", 11), ("hfmdc", 12), ("sigman", 13),\
        ("sigmap", 14), ("xbmodsp", 15), ("Cp", 16), ("KSE", 17), ("SEon",\
        18), ("SLmax", 19), ("SLmin", 20), ("SLrest", 21), ("SLset", 22),\
        ("b_ff", 23), ("fixed_afterload", 24), ("kxb", 25), ("visc", 26),\
        ("Ca_amplitude", 27), ("Ca_diastolic", 28), ("start_time", 29),\
        ("tau1", 30), ("tau2", 31), ("TmpC", 32), ("eps_0", 33),\
        ("len_hbare", 34), ("len_thick", 35), ("len_thin", 36), ("Qkn_p",\
        37), ("Qkoff", 38), ("Qkon", 39), ("Qkp_n", 40), ("kn_p", 41),\
        ("koffH", 42), ("koffL", 43), ("koffmod", 44), ("kon", 45), ("kp_n",\
        46), ("nperm", 47), ("perm50", 48), ("xPsi", 49))

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
    monitor_inds = dict([("eps_xb", 0), ("fappT", 1), ("gapslmd", 2),\
        ("gappT", 3), ("hfmd", 4), ("hbmd", 5), ("hfT", 6), ("hbT", 7),\
        ("gxbmd", 8), ("gxbT", 9), ("SSXBpostr", 10), ("Fnordv", 11),\
        ("force", 12), ("active", 13), ("ppforce", 14), ("preload", 15),\
        ("afterload", 16), ("total_force", 17), ("dlmbda", 18), ("beta", 19),\
        ("Cai", 20), ("konT", 21), ("koffLT", 22), ("koffHT", 23),\
        ("dTRPNCaL", 24), ("dTRPNCaH", 25), ("Tropreg", 26), ("permtot", 27),\
        ("inprmt", 28), ("kn_pT", 29), ("kp_nT", 30), ("dXBpostr", 31), ("P",\
        32), ("dXBprer", 33), ("dutyfrac", 34), ("SL", 35), ("sovr_ze", 36),\
        ("sovr_cle", 37), ("len_sovr", 38), ("SOVFThick", 39), ("SOVFThin",\
        40), ("dlmbda_dt", 41), ("dTRPNCaH_dt", 42), ("dTRPNCaL_dt", 43),\
        ("dN_dt", 44), ("dN_NoXB_dt", 45), ("dP_NoXB_dt", 46),\
        ("dXBpostr_dt", 47), ("dXBprer_dt", 48), ("dlmbda_a_dt", 49)])

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
    Compute the right hand side of the rice_lmbda ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 9)
    lmbda, TRPNCaH, TRPNCaL, N, N_NoXB, P_NoXB, XBpostr, XBprer, lmbda_a =\
        states

    # Assign parameters
    assert(len(parameters) == 50)
    Qfapp, Qgapp, Qgxb, Qhb, Qhf, fapp, gapp, gslmod, gxb, hb, hbmdc, hf,\
        hfmdc, sigman, sigmap, xbmodsp, Cp, KSE, SEon, SLmax, SLmin, SLrest,\
        SLset, b_ff, fixed_afterload, kxb, visc, Ca_amplitude, Ca_diastolic,\
        start_time, tau1, tau2, TmpC, eps_0, len_hbare, len_thick, len_thin,\
        Qkn_p, Qkoff, Qkon, Qkp_n, kn_p, koffH, koffL, koffmod, kon, kp_n,\
        nperm, perm50, xPsi = parameters

    # Init return args
    if values is None:
        values = np.zeros((9,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (9,)

    # Expressions for the Sarcomere geometry component
    SL = SLrest*lmbda
    sovr_ze = (len_thick/2 if len_thick/2 < SL/2 else SL/2)
    sovr_cle = (-SL/2 + len_thin if -SL/2 + len_thin > len_hbare/2 else\
        len_hbare/2)
    len_sovr = sovr_ze - sovr_cle
    SOVFThick = 2*len_sovr/(-len_hbare + len_thick)
    SOVFThin = len_sovr/len_thin

    # Expressions for the Thin filament regulation and crossbridge cycling
    # rates component
    eps_xb = 0.5*lmbda - 0.5*lmbda_a
    fappT = fapp*xbmodsp*math.pow(Qfapp, -37/10 + TmpC/10)
    gapslmd = 1 + gslmod*(1 - SOVFThick)
    gappT = gapp*xbmodsp*math.pow(Qgapp, -37/10 + TmpC/10)*gapslmd
    hfmd = math.exp(-hfmdc*(eps_xb*eps_xb)*math.copysign(1.0,\
        eps_xb)/(eps_0*eps_0))
    hbmd = math.exp(hbmdc*(eps_xb*eps_xb)*math.copysign(1.0,\
        eps_xb)/(eps_0*eps_0))
    hfT = hf*xbmodsp*math.pow(Qhf, -37/10 + TmpC/10)*hfmd
    hbT = hb*xbmodsp*math.pow(Qhb, -37/10 + TmpC/10)*hbmd
    gxbmd = (math.exp(sigmap*(eps_xb*eps_xb)/(eps_0*eps_0)) if eps_xb < 0 else\
        math.exp(sigman*(eps_xb*eps_xb)/(eps_0*eps_0)))
    gxbT = gxb*xbmodsp*math.pow(Qgxb, -37/10 + TmpC/10)*gxbmd

    # Expressions for the Equation for simulated calcium transient component
    beta = -math.pow(tau1/tau2, -1/(1 - tau2/tau1)) + math.pow(tau1/tau2,\
        -1/(-1 + tau1/tau2))
    Cai = ((-Ca_diastolic + Ca_amplitude)*(math.exp((start_time - t)/tau1) -\
        math.exp((start_time - t)/tau2))/beta + Ca_diastolic if t >\
        start_time else Ca_diastolic)

    # Expressions for the Ca binding to troponin to thin filament regulation
    # component
    konT = kon*math.pow(Qkon, -37/10 + TmpC/10)
    koffLT = koffL*koffmod*math.pow(Qkoff, -37/10 + TmpC/10)
    koffHT = koffH*koffmod*math.pow(Qkoff, -37/10 + TmpC/10)
    dTRPNCaL = (1 - TRPNCaL)*Cai*konT - TRPNCaL*koffLT
    dTRPNCaH = (1 - TRPNCaH)*Cai*konT - TRPNCaH*koffHT
    Tropreg = (1 - SOVFThin)*TRPNCaL + SOVFThin*TRPNCaH
    permtot = math.sqrt(math.fabs(1.0/(1 + math.pow(perm50/Tropreg, nperm))))
    inprmt = (1.0/permtot if 1.0/permtot < 100 else 100)
    values[2] = dTRPNCaL
    values[1] = dTRPNCaH
    kn_pT = kn_p*math.pow(Qkn_p, -37/10 + TmpC/10)*permtot
    kp_nT = kp_n*math.pow(Qkp_n, -37/10 + TmpC/10)*inprmt

    # Expressions for the Regulation and crossbridge cycling state equations
    # component
    values[4] = -N_NoXB*kn_pT + P_NoXB*kp_nT
    values[5] = -P_NoXB*kp_nT + N_NoXB*kn_pT
    dXBpostr = -XBpostr*hbT - XBpostr*gxbT + XBprer*hfT
    P = 1 - XBprer - XBpostr - N
    values[3] = P*kp_nT - N*kn_pT
    dXBprer = XBpostr*hbT - XBprer*hfT - XBprer*gappT + P*fappT
    values[6] = dXBpostr
    values[7] = dXBprer

    # Expressions for the Normalised active and passive force component
    SSXBpostr = fapp*hf/(gapp*hb + gapp*gxb + fapp*hb + gxb*hf + fapp*gxb +\
        fapp*hf)
    Fnordv = eps_0*kxb*SSXBpostr
    force = kxb*((eps_xb + eps_0)*XBpostr + XBprer*eps_xb)*SOVFThick
    active = force/Fnordv
    ppforce = Cp*b_ff*(-0.5 + 0.5*(lmbda*lmbda))*math.exp(b_ff*((-0.5 +\
        0.5*(lmbda*lmbda))*(-0.5 + 0.5*(lmbda*lmbda))))
    preload = Cp*b_ff*(-0.5 +\
        0.5*(SLset*SLset)/(SLrest*SLrest))*math.exp(b_ff*((-0.5 +\
        0.5*(SLset*SLset)/(SLrest*SLrest))*(-0.5 +\
        0.5*(SLset*SLset)/(SLrest*SLrest))))
    afterload = (KSE*(SLset - SL) if SEon == 1 else fixed_afterload)
    total_force = afterload - ppforce - active + preload
    dlmbda = (total_force/visc if SL > SLmin and SL <= SLmax else 0)
    values[0] = dlmbda

    # Expressions for the Mean strain of strongly bound states component
    dutyfrac = (fappT*gxbT + fappT*hbT + fappT*hfT)/(gxbT*hfT + fappT*gxbT +\
        fappT*hbT + gappT*gxbT + fappT*hfT + gappT*hbT)
    values[8] = -2*xPsi*(-0.5*lmbda + 0.5*lmbda_a)*fappT/dutyfrac

    # Return results
    return values

def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the rice_lmbda ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 9)
    lmbda, TRPNCaH, TRPNCaL, N, N_NoXB, P_NoXB, XBpostr, XBprer, lmbda_a =\
        states

    # Assign parameters
    assert(len(parameters) == 50)
    Qfapp, Qgapp, Qgxb, Qhb, Qhf, fapp, gapp, gslmod, gxb, hb, hbmdc, hf,\
        hfmdc, sigman, sigmap, xbmodsp, Cp, KSE, SEon, SLmax, SLmin, SLrest,\
        SLset, b_ff, fixed_afterload, kxb, visc, Ca_amplitude, Ca_diastolic,\
        start_time, tau1, tau2, TmpC, eps_0, len_hbare, len_thick, len_thin,\
        Qkn_p, Qkoff, Qkon, Qkp_n, kn_p, koffH, koffL, koffmod, kon, kp_n,\
        nperm, perm50, xPsi = parameters

    # Init return args
    if monitored is None:
        monitored = np.zeros((50,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (50,)

    # Expressions for the Sarcomere geometry component
    monitored[35] = SLrest*lmbda
    monitored[36] = (len_thick/2 if len_thick/2 < monitored[35]/2 else\
        monitored[35]/2)
    monitored[37] = (len_thin - monitored[35]/2 if len_thin - monitored[35]/2 >\
        len_hbare/2 else len_hbare/2)
    monitored[38] = monitored[36] - monitored[37]
    monitored[39] = 2*monitored[38]/(-len_hbare + len_thick)
    monitored[40] = monitored[38]/len_thin

    # Expressions for the Thin filament regulation and crossbridge cycling
    # rates component
    monitored[0] = 0.5*lmbda - 0.5*lmbda_a
    monitored[1] = fapp*xbmodsp*math.pow(Qfapp, -37/10 + TmpC/10)
    monitored[2] = 1 + gslmod*(1 - monitored[39])
    monitored[3] = gapp*xbmodsp*math.pow(Qgapp, -37/10 + TmpC/10)*monitored[2]
    monitored[4] =\
        math.exp(-hfmdc*(monitored[0]*monitored[0])*math.copysign(1.0,\
        monitored[0])/(eps_0*eps_0))
    monitored[5] =\
        math.exp(hbmdc*(monitored[0]*monitored[0])*math.copysign(1.0,\
        monitored[0])/(eps_0*eps_0))
    monitored[6] = hf*xbmodsp*math.pow(Qhf, -37/10 + TmpC/10)*monitored[4]
    monitored[7] = hb*xbmodsp*math.pow(Qhb, -37/10 + TmpC/10)*monitored[5]
    monitored[8] =\
        (math.exp(sigmap*(monitored[0]*monitored[0])/(eps_0*eps_0)) if\
        monitored[0] < 0 else\
        math.exp(sigman*(monitored[0]*monitored[0])/(eps_0*eps_0)))
    monitored[9] = gxb*xbmodsp*math.pow(Qgxb, -37/10 + TmpC/10)*monitored[8]

    # Expressions for the Equation for simulated calcium transient component
    monitored[19] = -math.pow(tau1/tau2, -1/(1 - tau2/tau1)) +\
        math.pow(tau1/tau2, -1/(-1 + tau1/tau2))
    monitored[20] = (Ca_diastolic + (-Ca_diastolic +\
        Ca_amplitude)*(math.exp((start_time - t)/tau1) - math.exp((start_time\
        - t)/tau2))/monitored[19] if t > start_time else Ca_diastolic)

    # Expressions for the Ca binding to troponin to thin filament regulation
    # component
    monitored[21] = kon*math.pow(Qkon, -37/10 + TmpC/10)
    monitored[22] = koffL*koffmod*math.pow(Qkoff, -37/10 + TmpC/10)
    monitored[23] = koffH*koffmod*math.pow(Qkoff, -37/10 + TmpC/10)
    monitored[24] = (1 - TRPNCaL)*monitored[20]*monitored[21] -\
        TRPNCaL*monitored[22]
    monitored[25] = (1 - TRPNCaH)*monitored[20]*monitored[21] -\
        TRPNCaH*monitored[23]
    monitored[26] = (1 - monitored[40])*TRPNCaL + TRPNCaH*monitored[40]
    monitored[27] = math.sqrt(math.fabs(1.0/(1 +\
        math.pow(perm50/monitored[26], nperm))))
    monitored[28] = (1.0/monitored[27] if 1.0/monitored[27] < 100 else 100)
    monitored[43] = monitored[24]
    monitored[42] = monitored[25]
    monitored[29] = kn_p*math.pow(Qkn_p, -37/10 + TmpC/10)*monitored[27]
    monitored[30] = kp_n*math.pow(Qkp_n, -37/10 + TmpC/10)*monitored[28]

    # Expressions for the Regulation and crossbridge cycling state equations
    # component
    monitored[45] = P_NoXB*monitored[30] - N_NoXB*monitored[29]
    monitored[46] = N_NoXB*monitored[29] - P_NoXB*monitored[30]
    monitored[31] = -XBpostr*monitored[9] - XBpostr*monitored[7] +\
        XBprer*monitored[6]
    monitored[32] = 1 - XBprer - XBpostr - N
    monitored[44] = monitored[30]*monitored[32] - N*monitored[29]
    monitored[33] = -XBprer*monitored[6] + monitored[1]*monitored[32] +\
        XBpostr*monitored[7] - XBprer*monitored[3]
    monitored[47] = monitored[31]
    monitored[48] = monitored[33]

    # Expressions for the Normalised active and passive force component
    monitored[10] = fapp*hf/(gapp*hb + gapp*gxb + fapp*hb + gxb*hf + fapp*gxb\
        + fapp*hf)
    monitored[11] = eps_0*kxb*monitored[10]
    monitored[12] = kxb*((eps_0 + monitored[0])*XBpostr +\
        XBprer*monitored[0])*monitored[39]
    monitored[13] = monitored[12]/monitored[11]
    monitored[14] = Cp*b_ff*(-0.5 + 0.5*(lmbda*lmbda))*math.exp(b_ff*((-0.5 +\
        0.5*(lmbda*lmbda))*(-0.5 + 0.5*(lmbda*lmbda))))
    monitored[15] = Cp*b_ff*(-0.5 +\
        0.5*(SLset*SLset)/(SLrest*SLrest))*math.exp(b_ff*((-0.5 +\
        0.5*(SLset*SLset)/(SLrest*SLrest))*(-0.5 +\
        0.5*(SLset*SLset)/(SLrest*SLrest))))
    monitored[16] = (KSE*(SLset - monitored[35]) if SEon == 1 else\
        fixed_afterload)
    monitored[17] = monitored[16] - monitored[13] - monitored[14] +\
        monitored[15]
    monitored[18] = (monitored[17]/visc if monitored[35] > SLmin and\
        monitored[35] <= SLmax else 0)
    monitored[41] = monitored[18]

    # Expressions for the Mean strain of strongly bound states component
    monitored[34] = (monitored[1]*monitored[7] + monitored[1]*monitored[6] +\
        monitored[1]*monitored[9])/(monitored[3]*monitored[9] +\
        monitored[1]*monitored[7] + monitored[3]*monitored[7] +\
        monitored[1]*monitored[6] + monitored[6]*monitored[9] +\
        monitored[1]*monitored[9])
    monitored[49] = -2*xPsi*(-0.5*lmbda +\
        0.5*lmbda_a)*monitored[1]/monitored[34]

    # Return results
    return monitored
