function [nVAvalveOp,nVAvalveCl,nAVvalveOp,nAVvalveCl] = ValveEvents(P,nREF,LR)

t           = (1:length(P.t))';
tCycle      = t(end);
ModC        = @(ti) mod(ti,length(t));
reref       = [nREF:length(t) 1:nREF-1];

if LR == 'L'
    ind     = 1;
elseif LR == 'R'
    ind     = 2;
end

VAvalve     = {'LvSyArt' 'RvPuArt'};
AVvalve     = {'LaLv' 'RaRv'};
Atr         = {'La1' 'Ra1'};

%% Ventriculoarterial valve flow
q_aux       = GetFt('Valve','q',VAvalve(ind));
q           = q_aux(reref,:);
qDot_aux    = GetFt('Valve','qDot',VAvalve(ind));
qDot        = qDot_aux(reref,:);

% nAVO and nAVC extrapolation slopes 10%-40% of max flow
qMax        = max(q);
nqMax       = find(q==qMax);
Rg          = find(0.10*qMax<q & q<0.4*qMax & qDot>0);
tRg         = ModRound(t(Rg)-nqMax,t(end));
M           = [0*tRg+1,tRg];
ab          = M\q(Rg);
nVAvalveOp  = round(ModC(nqMax-ab(1)/ab(2)));
Rg          = find(0.10*qMax<q & q<0.4*qMax & qDot<0);
tRg         = ModRound(t(Rg)-nqMax,t(end));
M           = [0*tRg+1,tRg];
ab          = M\q(Rg);
nVAvalveCl  = round(ModC(nqMax-ab(1)/ab(2)));

% figure; 
% subplot(3,1,1), plot(t,q,'b'), hold on, plot(t([nVAvalveOp,nVAvalveCl]),q([nVAvalveOp,nVAvalveCl]),'ob')
% subplot(3,1,2), plot(t,qDot,'b'), hold on, plot(t([nVAvalveOp,nVAvalveCl]),qDot([nVAvalveOp,nVAvalveCl]),'ob')

%% Atrioventricular valve flow
q_aux       = GetFt('Valve','q',AVvalve(ind));
q           = q_aux(reref,:);
qDot_aux    = GetFt('Valve','qDot',AVvalve(ind));
qDot        = qDot_aux(reref,:);

% Atrial Activation
C_aux       = GetFt('Patch','C',Atr(ind));
C           = C_aux(reref,:);
CDot_aux    = GetFt('Patch','CDot',Atr(ind));
CDot        = CDot_aux(reref,:);
WA          = C/max(C); % weight function of atrial contractility
tCR         = tZero(WA-0.1,t,nVAvalveCl); % rise of atrial C
tCMax       = tZero(CDot,t,tCR); % peak atrial activation time

% Atrioventricular valve flow
tSmear      = 0.02*tCycle; nSmear=tSmear; % no smear of flow
qAv         = LowPass(q        ,nSmear);
qDot        = LowPass(qDot     ,nSmear);
qE1         = LowPass(q.*(1-WA),nSmear); %initial estimate E-Peak

% Analysis of E-peak, timing of flow events
nAVvalveCl  = round(tZero(qAv-0.01*max(qAv),t,tCMax));
t1          = tZero(qE1-0.2*max(qE1),t,nVAvalveCl); % after mitral valve opening
tAvPeakE    = tZero(qDot,t,t1); % E-peak flow
qAvMaxE     = Ft(qE1,t,tAvPeakE); % early peak flow
t1          = tZero(qE1-0.2*qAvMaxE,t,nVAvalveCl); % after valve opening
t2          = tZero(qE1-0.5*qAvMaxE,t,t1);
t3          = tZero(qE1-0.8*qAvMaxE,t,t2);
nAVvalveOp  = round(ModC(t1-0.33*ModC(t3-t1))); % valve opening
% vAvMaxE     = max(q(round(length(q)/2):end))/GetFt('Valve','AOpen',ValveAv); % E-peak velocity

% subplot(3,1,1), plot(t,q,'r'), hold on, plot(t([nAVvalveOp,nAVvalveCl]),q([nAVvalveOp,nAVvalveCl]),'or')
% subplot(3,1,2), plot(t,qDot,'r'), hold on, plot(t([nAVvalveOp,nAVvalveCl]),qDot([nAVvalveOp,nAVvalveCl]),'or')

%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tm= ModRound(t,T)
%modulus around zero value continuous (-0.5 -> 0 -> +0.5)
tm= t-round(t/T)*T;
return

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

