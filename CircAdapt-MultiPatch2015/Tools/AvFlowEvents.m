function [tAvOpen, tAvClose, tAvDecay, vAvMaxE, vAvMaxA]=AvFlowEvents(P, ValveAv, tArtClose)

% Function definition
tCycle = P.General.tCycle;
ModC= @(ti) mod(ti,tCycle);
Dt= P.General.Dt;

% Atrioventricular valve flow
q    = GetFt('Valve','q',ValveAv);
qDot = GetFt('Valve','qDot',ValveAv);
t    = (tCycle/length(q))*[0:length(q)-1]';

% Atrial Activation
C= GetFt('Patch','C',P.Cavity.Name(GetFt('Valve','iNodeProx',ValveAv)));
CDot= GetFt('Patch','CDot',P.Cavity.Name(GetFt('Valve','iNodeProx',ValveAv)));
WA= C/max(C); % weight function of atrial contractility
%     tCR   = tZero(C-0.1*max(C),t,tLArtCloseM); %rise of atrial C
tCR= tZero(WA-0.1,t,tArtClose); %rise of atrial C
tCMax = tZero(CDot,t,tCR); % peak atrial activation time

% Atrioventricular valve flow
tSmear=0.04*tCycle; nSmear=tSmear/Dt; %slight smear of flow
qAv    = LowPass(q        ,nSmear);
qDot   = LowPass(qDot     ,nSmear);
qE1    = LowPass(q.*(1-WA),nSmear); %initial estimate E-Peak
qE1Dot = conv2([qE1(end);qE1;qE1(1)],[+0.5;0;-0.5],'valid')/Dt;

% Analysis of E-peak, timing of flow events
tAvClose= tZero(qAv-0.01*max(qAv),t,tCMax);
t1= tZero(qE1-0.2*max(qE1),t,tArtClose); % after mitral valve opening
tAvPeakE= tZero(qDot,t,t1); % E-peak flow
qAvMaxE = Ft(qE1,t,tAvPeakE); % early peak flow
t1= tZero(qE1-0.2*qAvMaxE,t,tArtClose); % after valve opening
t2= tZero(qE1-0.5*qAvMaxE,t,t1);
t3= tZero(qE1-0.8*qAvMaxE,t,t2);
tAvOpen = ModC(t1-0.33*ModC(t3-t1)); % valve opening
vAvMaxE = max(q(round(length(q)/2):end))/GetFt('Valve','AOpen',ValveAv); % E-peak velocity

% Diastasis moment, estimate tail of E-Peak
Diastasis= ( ModC(t-tAvPeakE)> 0 & ...
    ModC(t-tAvPeakE)<ModC(tCR-tAvPeakE) );
qDotMin= min(qE1Dot.*Diastasis); % maximum decaying slope in diastasis
RgDiast= find(qDotMin==qE1Dot.*Diastasis); RgDiast=RgDiast(1);
tDiast= t(RgDiast); qDiast=qAv(RgDiast); % time and flow at diastasis
tD=ModC(t-tDiast); tDMax= -1.4*qDiast/qDotMin; %1.4 correction
Rg= find( tD>0 & tD<tDMax );
qE=qAv; qE(Rg)= qDiast*(1-tD(Rg)/tDMax);
W0= ( ModC(t-tDiast-tDMax) < ModC(t-tAvClose) );
qE=qE.*(1-W0); %zero tail of qE

qE=min(qAv,qE); qA=qAv-qE; % estimate of E and A component
qADot= conv2([qA(end);qA;qA(1)],[+0.5;0;-0.5],'valid')/Dt;

tAvDecay  = ModC(tDiast+tDMax); % end of E-peak
tAvPeakA  = tZero(qADot,t,tCR+tSmear+0.01); % time of A-peak
qAvMaxA   = Ft(qAv,t,tAvPeakA); % Flow at A-Peak
vAvMaxA   = max(q(1:round(length(q)/2)))/GetFt('Valve','AOpen',ValveAv); %A-peak velocity

return

%==========================================================================
% Auxiliary functions 
%==========================================================================

function tZeroC= tZero(q,t,t0) %first sign change after t0
%q(t)= a signal array of a complete cycle
%t   = time array
%t0  = starting time moment (double)
%output=
n= length(t); T= (t(end)-t(1))*n/(n-1);
Signq= double(q>0);
ZeroCross= conv2([Signq(end);Signq],[+1;-1],'valid');
RgZero1= find(ZeroCross~= 0);
if isempty(RgZero1)
    tZeroC= 0;
else
    RgZero2= RgZero1+1;
    qp= [q(end);q]; tp= [-2*t(1)+t(2);t];
    q1= qp(RgZero1);q2= qp(RgZero2);
    di= q1./(q1-q2);
    tZeroX= tp(RgZero1)+di.*(tp(RgZero2)-tp(RgZero1));
    d= tZeroX-t0-floor((tZeroX-t0)/T)*T;
    Rg= find(d== min(d));
    tZeroC= tZeroX(Rg(1));
end
return

function ft0= Ft(f,t,t0) %interpolated function value at t= t0
f= [f(end);f;f(1)];
t= [2*t(1)-t(2);t;2*t(end)-t(end-1)];
ft0= interp1(t,f,t0);
return

function fOut=LowPass(fIn,nLf)
% nLf=wave length cuttoff
% makes function cyclic and smears with exp(-x2)
F1=fft(fIn);
n=length(fIn); w=[0:n-1]';
w0=n/nLf;
Filter=2*exp(-(w/w0).^2); Filter(1)=1;
fOut=real(ifft(F1.*Filter));
return

