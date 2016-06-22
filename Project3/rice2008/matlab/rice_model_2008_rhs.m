function [values] = rice_model_2008_rhs(t, states, parameters)
  % Compute the right hand side of the rice_model_2008 ODE

  % Assign states
  if length(states)~=11
    error('Expected the states array to be of size 11.');
  end
  intf=states(1); SL=states(2); TRPNCaL=states(3); TRPNCaH=states(4);...
    N_NoXB=states(5); P_NoXB=states(6); N=states(7); XBpostr=states(8);...
    XBprer=states(9); xXBprer=states(10); xXBpostr=states(11);

  % Assign parameters
  if length(parameters)~=56
    error('Expected the parameters array to be of size 56.');
  end
  Qfapp=parameters(1); Qgapp=parameters(2); Qgxb=parameters(3);...
    Qhb=parameters(4); Qhf=parameters(5); fapp=parameters(6);...
    gapp=parameters(7); gslmod=parameters(8); gxb=parameters(9);...
    hb=parameters(10); hbmdc=parameters(11); hf=parameters(12);...
    hfmdc=parameters(13); sigman=parameters(14); sigmap=parameters(15);...
    xbmodsp=parameters(16); KSE=parameters(17); PCon_c=parameters(18);...
    PCon_t=parameters(19); PExp_c=parameters(20); PExp_t=parameters(21);...
    SEon=parameters(22); SL_c=parameters(23); SLmax=parameters(24);...
    SLmin=parameters(25); SLrest=parameters(26); SLset=parameters(27);...
    fixed_afterload=parameters(28); kxb_normalised=parameters(29);...
    massf=parameters(30); visc=parameters(31); Ca_amplitude=parameters(32);...
    Ca_diastolic=parameters(33); start_time=parameters(34);...
    tau1=parameters(35); tau2=parameters(36); TmpC=parameters(37);...
    len_hbare=parameters(38); len_thick=parameters(39);...
    len_thin=parameters(40); x_0=parameters(41); Qkn_p=parameters(42);...
    Qkoff=parameters(43); Qkon=parameters(44); Qkp_n=parameters(45);...
    kn_p=parameters(46); koffH=parameters(47); koffL=parameters(48);...
    koffmod=parameters(49); kon=parameters(50); kp_n=parameters(51);...
    nperm=parameters(52); perm50=parameters(53); xPsi=parameters(54);

  % Init return args
  values = zeros(11, 1);

  % Expressions for the Sarcomere geometry component
  sovr_ze = ((len_thick/2 < SL/2)*(len_thick/2) + ~(len_thick/2 <...
    SL/2)*(SL/2));
  sovr_cle = ((len_thin - SL/2 > len_hbare/2)*(len_thin - SL/2) + ~(len_thin...
    - SL/2 > len_hbare/2)*(len_hbare/2));
  len_sovr = -sovr_cle + sovr_ze;
  SOVFThick = 2*len_sovr/(len_thick - len_hbare);
  SOVFThin = len_sovr/len_thin;

  % Expressions for the Thin filament regulation and crossbridge cycling
  % rates component
  fappT = fapp*xbmodsp*Qfapp^(-37/10 + TmpC/10);
  gapslmd = 1 + gslmod*(1 - SOVFThick);
  gappT = gapp*xbmodsp*Qgapp^(-37/10 + TmpC/10)*gapslmd;
  hfmd = exp(-hfmdc*xXBprer^2*sign(xXBprer)/x_0^2);
  hbmd = exp(hbmdc*(-x_0 + xXBpostr)^2*sign(-x_0 + xXBpostr)/x_0^2);
  hfT = hf*xbmodsp*Qhf^(-37/10 + TmpC/10)*hfmd;
  hbT = hb*xbmodsp*Qhb^(-37/10 + TmpC/10)*hbmd;
  gxbmd = ((xXBpostr < x_0)*(exp(sigmap*(x_0 - xXBpostr)^2/x_0^2)) +...
    ~(xXBpostr < x_0)*(exp(sigman*(-x_0 + xXBpostr)^2/x_0^2)));
  gxbT = gxb*xbmodsp*Qgxb^(-37/10 + TmpC/10)*gxbmd;

  % Expressions for the Normalised active and passive force component
  SSXBpostr = fapp*hf/(fapp*gxb + fapp*hb + fapp*hf + gapp*gxb + gapp*hb +...
    gxb*hf);
  Fnordv = kxb_normalised*x_0*SSXBpostr;
  force = kxb_normalised*(XBpostr*xXBpostr + XBprer*xXBprer)*SOVFThick;
  active = force/Fnordv;
  ppforce_t = PCon_t*(-1 + exp(PExp_t*abs(SLrest - SL)))*sign(-SLrest + SL);
  ppforce_c = ((SL > SL_c)*(PCon_c*(-1 + exp(PExp_c*abs(SL_c - SL)))) + ~(SL...
    > SL_c)*(0));
  ppforce = ppforce_c + ppforce_t;
  preload = PCon_t*(-1 + exp(PExp_t*abs(SLrest - SLset)))*sign(SLset - SLrest);
  afterload = ((SEon == 1)*(KSE*(SLset - SL)) + ~(SEon ==...
    1)*(fixed_afterload));
  dSL = ((SL > SLmin & SL <= SLmax)*((visc*(SLset - SL) + intf)/massf) + ~(SL...
    > SLmin & SL <= SLmax)*(0));
  values(1) = -active - ppforce + afterload + preload;
  values(2) = dSL;

  % Expressions for the Equation for simulated calcium transient component
  beta = (tau1/tau2)^(-1/(-1 + tau1/tau2)) - (tau1/tau2)^(-1/(1 - tau2/tau1));
  Cai = ((t > start_time)*(Ca_diastolic + (Ca_amplitude -...
    Ca_diastolic)*(-exp((start_time - t)/tau2) + exp((start_time -...
    t)/tau1))/beta) + ~(t > start_time)*(Ca_diastolic));

  % Expressions for the Ca binding to troponin to thin filament regulation
  % component
  konT = kon*Qkon^(-37/10 + TmpC/10);
  koffLT = koffL*koffmod*Qkoff^(-37/10 + TmpC/10);
  koffHT = koffH*koffmod*Qkoff^(-37/10 + TmpC/10);
  dTRPNCaL = -TRPNCaL*koffLT + (1 - TRPNCaL)*Cai*konT;
  dTRPNCaH = -TRPNCaH*koffHT + (1 - TRPNCaH)*Cai*konT;
  Tropreg = (1 - SOVFThin)*TRPNCaL + SOVFThin*TRPNCaH;
  permtot = sqrt(abs(1.0/(1 + (perm50/Tropreg)^nperm)));
  inprmt = ((1.0/permtot < 100)*(1.0/permtot) + ~(1.0/permtot < 100)*(100));
  values(3) = dTRPNCaL;
  values(4) = dTRPNCaH;
  kn_pT = kn_p*Qkn_p^(-37/10 + TmpC/10)*permtot;
  kp_nT = kp_n*Qkp_n^(-37/10 + TmpC/10)*inprmt;

  % Expressions for the Regulation and crossbridge cycling state equations
  % component
  values(5) = P_NoXB*kp_nT - N_NoXB*kn_pT;
  values(6) = N_NoXB*kn_pT - P_NoXB*kp_nT;
  dXBpostr = XBprer*hfT - XBpostr*gxbT - XBpostr*hbT;
  P = 1 - N - XBpostr - XBprer;
  values(7) = P*kp_nT - N*kn_pT;
  dXBprer = P*fappT + XBpostr*hbT - XBprer*gappT - XBprer*hfT;
  values(8) = dXBpostr;
  values(9) = dXBprer;

  % Expressions for the Mean strain of strongly bound states component
  dutyprer = (fappT*gxbT + fappT*hbT)/(fappT*gxbT + fappT*hbT + fappT*hfT +...
    gappT*gxbT + gappT*hbT + gxbT*hfT);
  dutypostr = fappT*hfT/(fappT*gxbT + fappT*hbT + fappT*hfT + gappT*gxbT +...
    gappT*hbT + gxbT*hfT);
  dxXBprer = dSL/2 + xPsi*((-x_0 - xXBprer + xXBpostr)*hbT -...
    fappT*xXBprer)/dutyprer;
  dxXBpostr = dSL/2 + xPsi*(x_0 - xXBpostr + xXBprer)*hfT/dutypostr;
  values(10) = dxXBprer;
  values(11) = dxXBpostr;
end