import math

def force_transient(t):
    tau1=20
    tau2=110
    start_time = 5
    F_amplitude = 1.0
    F_diastolic = 0.0
    F = F_amplitude*(-math.exp((start_time - t)/tau1) +\
        math.exp((start_time - t)/tau2)) if t > start_time else F_diastolic
    return F
