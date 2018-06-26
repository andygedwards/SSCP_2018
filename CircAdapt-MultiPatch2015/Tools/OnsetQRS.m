function [OnsQRS]=OnsetQRS(P);

% Function definition
tCycle = P.General.tCycle;
ModC= @(ti) mod(ti,tCycle);

% Contractilities
CL    = sum([GetFt('Patch','C','Lv') GetFt('Patch','C','Sv') GetFt('Patch','C','Rv')],2);
CDotL = sum([GetFt('Patch','CDot','Lv') GetFt('Patch','CDot','Sv') GetFt('Patch','CDot','Rv')],2);
% Maximum contractilities
CMax= max(CL);

t    = (tCycle/length(CL))*[0:length(CL)-1]';

% Onset activation extrapolation slope 10%-40% of max contractility
tCLMax= t(find(CL == max(CL)));
RgL= find(0.10*CMax<CL & CL<0.40*CMax & CDotL>0);
tRgL= ModRound(t(RgL)-tCLMax,tCycle);
ML= [0*tRgL+1,tRgL];
abL= ML\CL(RgL);
OnsActL= ModC(tCLMax-abL(1)/abL(2));

OnsQRS=OnsActL;

return


%==========================================================================
% Auxiliary functions 
%==========================================================================

function tm= ModRound(t,T);
%modulus around zero value continuous (-0.5 -> 0 -> +0.5)
tm= t-round(t/T)*T;
return

