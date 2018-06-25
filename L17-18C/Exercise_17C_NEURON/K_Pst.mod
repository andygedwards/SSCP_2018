:Comment : The persistent component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21


NEURON	{
	SUFFIX K_Pst
	USEION k READ ek WRITE ik
	RANGE gK_Pstbar, gK_Pst, m, h, ik, offm, slom, offmt, slomt, taummin, taumdiff1, taumdiff2, offh, sloh, offht1, offht2, sloht, tauhmean, tauhdiff1, tauhdiff2
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gK_Pstbar = 0.00001 (S/cm2)
	offm = -11 (mV)
	slom = 12 (mV)
	offmt = -10 (mV)
	slomt = 38.46153846 (mV)
	taummin = 1.25 (ms)
        taumdiff1 = 175.03 (ms)
        taumdiff2 = 13 (ms)
	offh = -64 (mV)
	sloh = 11 (mV)
	offht1 = -65 (mV)
	offht2 = -85 (mV)
	sloht = 48 (mV)
	tauhmean = 360 (ms)
	tauhdiff1 = 1010 (ms)
        tauhdiff2 = 24 (ms/mV)

}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gK_Pst	(S/cm2)
	mInf
	mTau
	hInf
	hTau
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gK_Pst = gK_Pstbar*m*m*h
	ik = gK_Pst*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
  LOCAL qt, thresh
  qt = 2.3^((34-21)/10)
  thresh = offmt-slomt/2*log(taumdiff1/taumdiff2)
	UNITSOFF
		mInf =  (1/(1 + exp((offm-v)/slom)))
                if(v<thresh){
		    mTau =  (taummin+taumdiff1*exp(-(offmt-v)/slomt))/qt
                } else {
                    mTau = ((taummin+taumdiff2*exp((offmt-v)/slomt)))/qt
                }
		hInf =  1/(1 + exp(-(offh-v)/sloh))
		hTau =  (tauhmean+(tauhdiff1-tauhdiff2*(offht1-v))*exp(-((offht2-v)/sloht)^2))/qt
	UNITSON
}
