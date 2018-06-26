function [tArtOpen, tArtClose, tArtMax, vArtMax, vArtMed]=ArtFlowEvents(P, ValveArt)

% Function definition
tCycle = P.General.tCycle;
ModC= @(ti) mod(ti,tCycle);
Dt= P.General.Dt;

% Arterial valve flow
q    = GetFt('Valve','q',ValveArt);
qDot = GetFt('Valve','qDot',ValveArt);
t    = (tCycle/length(q))*[0:length(q)-1]';

% Ventricular activation
C    = mean(P.Patch.C,2);
tAct = tCycle*mean(C)/max(C);

% Maximum flow velocity
qMax= max(q);
vArtMax= qMax/GetFt('Valve','AOpen',ValveArt); % peak flow velocity
Rg= find(q>0.85*qMax); tAux= t(Rg(1));
tArtMax= tZero(qDot,t,tAux); % time of maximum flow
tRg= ModRound(t(Rg)-tArtMax,tCycle);
M= [0*tRg+1,tRg,tRg.^2];
ab= M\q(Rg); % best quadratic fit
qMax= ab(1)-ab(2)^2/(4*ab(3));
tArtMax= tArtMax-0.5*ab(2)/ab(3); % improvement of peak time

% Opening valve extrapolation slope 10%-40% of max flow
Rg= find(0.20*qMax<q & q<0.7*qMax & qDot>0);
tRg= ModRound(t(Rg)-tArtMax,tCycle);
M= [0*tRg+1,tRg];
ab= M\q(Rg);
tArtOpen= ModC(tArtMax-ab(1)/ab(2));

% Closure valve
Rg= find(0.05*max(q)<q & q<0.6*max(q) & qDot<0);
tRg= ModRound(t(Rg)-tArtMax,tCycle);
M= [0*tRg+1,tRg];
ab= M\q(Rg);
tArtClose= ModC(tArtMax-ab(1)/ab(2));
tEj= ModC(tArtClose-tArtOpen);

% Median flow (velocity) during ejection
Rg2= find(q>0.1*qMax); %+++++++++++++++++
Aux= -sort(-q(Rg2));
n= tEj/(2*Dt); % half duration of ejection
n= max(1,min(length(Rg2)-2,n));
ni= floor(n); f= n-ni;
qArtMed= Aux(ni)+(Aux(ni+1)-Aux(ni))*f; % median flow
vArtMed= qArtMed/GetFt('Valve','AOpen',ValveArt); % median flow velocity

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

function tm= ModRound(t,T)
%modulus around zero value continuous (-0.5 -> 0 -> +0.5)
tm= t-round(t/T)*T;
return

