function ValveqDot
% function ValveqDot
% Node pressure differences -> valve flow acceleration qDot
% Currently, mitral valve (LaLv) and tricuspid valve (RaRv) close
% completely with backflow during diastole, still open to improvement.
% Theo Arts, Maastricht University, Oct 30, 2011

global P

iNodeProx=P.Valve.iNodeProx;
iNodeDist=P.Valve.iNodeDist;
nt=size(P.t,1); Col1=ones(nt,1);
AOpen = P.Valve.AOpen(Col1,:); % open valve cross-section
ALeak = P.Valve.ALeak(Col1,:); % closed valve cross-section

% Ws=P.Valve.AvWalls; Vs=P.Valve.AvValves;
Ws=[7 9]; Vs=[5 2]; % Ws=[]; Vs=[]; % empty arrays turns off papillary muscle function
T   = P.Wall.T(:,Ws); %wall tension related to pap. muscle
DADT= P.Wall.DADT(:,Ws); %wall stiffness
Am0 = P.Wall.Am0(:,Ws); %wall zero-stress area
Diast= tanh(100*max(0.001,T.*DADT./Am0-0.1).^2 ); % diastole->1.0
ALeak(:,Vs)= max( cat( 3, 0.1*AOpen(:,Vs).*Diast, ALeak(:,Vs) ), [ ], 3 ); % concatenate in third dimension so method still works on whole time series 

% Valve qDot
q     = P.Valve.q;
Len   = P.Valve.Len  ; % effective length of flow channel
rhob  = 1050         ; % density of blood
Dp    = P.Node.p(:,iNodeProx)-P.Node.p(:,iNodeDist); % pressure drop
iCavProx= P.Node.iCavity(P.Valve.iNodeProx); % Valve->prox Node->iCavity 
iCavDist= P.Node.iCavity(P.Valve.iNodeDist); % Valve->dist Node->iCavity
AProx = P.Cavity.A(:,iCavProx); % proximal area
ADist = P.Cavity.A(:,iCavDist); % distal area
AExt  = min(AProx,ADist); % minimum external orifice area, fu(t)
AOpenEff = min(AOpen,AExt); % open flow orifice area, fu(t)
ALeakEff = min(ALeak,AExt); % leak flow orifice area, fu(t)

% Slow valve closure to stabilize solution Diff.Eq.
Reverse  = double(AOpenEff<ALeakEff); % if valve acts backward
Normal  = (1-Reverse);
AFw= Normal.*AOpenEff+Reverse.*ALeakEff; 
ABw= Normal.*ALeakEff+Reverse.*AOpenEff;
APr= Normal.*AProx+Reverse.*ADist;
ADi= Normal.*ADist+Reverse.*AProx;
qV = (1-2*Reverse).* q;
DpV= (1-2*Reverse).*Dp;
Dp0= 20; % pressure [Pa]  as safety for zero division

x= 40.0*rhob*qV.*abs(qV)./(AFw.^2); % Factor ~ 1/closure time
y= DpV;
r= sqrt(x.^2+y.^2+Dp0^2);
yPos= double(y>0); % positive forward pressure
xPos= double(x>0); % positive forward flow
Closed= (1-xPos).*(1-yPos); % Closed for backward flow & pressure
Closing= xPos .*(1-yPos); % Closing phase by backward pressure
Open   = yPos; % Fully open if pressure forward
AClosing= sqrt(x./r).*(AFw-ABw)+ABw; % Soft closure described by
% Closure by a continuous function
A= Closed.*ABw + Open.*AFw + Closing.*AClosing;

% Bernouilli pressure drop
v1= qV./APr; v2= qV./A; v3= qV./ADi; % velocities
DpB= (0.5*rhob) * (1-2*Reverse) .* (...
    double(qV>0).*max(v2.^2-v1.^2,0) - double(qV<0).*max(v2.^2-v3.^2,0) );
L= 1.5*rhob*( Len(Col1,:)./A+ 0.5*(1./sqrt(APr)+1./sqrt(ADi))); % inertia
P.Valve.L=L; % inertia
P.Valve.qDot= (Dp-DpB)./L; % flow derivative

end

