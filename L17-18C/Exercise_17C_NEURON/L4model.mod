TITLE L4model.mod   sodium + potassium
 
COMMENT
  The exercise 5 from lecture 4
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
}
 
NEURON {
        SUFFIX L4model
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        RANGE gnabar, gkbar, h, gna
}
 
PARAMETER {
        gnabar = 0.8 (S/cm2)	<0,1e9>
        gkbar = 0.1 (S/cm2)	<0,1e9>
        Vam = -72 (mV)
        Vbm = -10 (mV)
        Vah = -80 (mV)
        Vbh = -20 (mV)
        dam = 15 (mV)
        dah = -5 (mV)
        dbm = -15 (mV)
        dbh = 5 (mV)
}
 
STATE {
        h
}
 
ASSIGNED {
        v (mV)
        gna (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        minf
        ena
        ek
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*minf*h
        ina = gna*(v - ena)
        ik = gkbar*(v - ek)      
}
 
 
INITIAL {
        LOCAL  alpha, beta
        alpha = exp((v-Vam)/dam)
        beta = exp((v-Vbm)/dbm)
        minf = alpha/(alpha+beta)
        h = 0.5
}

DERIVATIVE states {  
        LOCAL  alpha, beta
        alpha = exp((v-Vam)/dam)
        beta = exp((v-Vbm)/dbm)
        minf = alpha/(alpha+beta)

        alpha = exp((v-Vah)/dah)
        beta = exp((v-Vbh)/dbh)
        h' = alpha*(1-h) - h*beta
}
 
UNITSON
