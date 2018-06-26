function AdaptRestP
% function AdaptRestP
% Simulates adaptation of vessels to hemodynamics at rest
% Adaptation of vessel cross-section
% Theo Arts, Maastricht University, Oct 30, 2011


global P
save PTemp P; %saves last intermediate solution

% 1st part like Adapt0P controls systemic blood pressure and flow
% by adjustment of circulatory blood volume and peripheral resistance
%
% Assessment of stationarity of flows
FlowVec=mean(P.Valve.q);
% test on presence of FlowVec with right size
if isfield(P.Adapt, 'FlowVec');
    FlowVecPrev= P.Adapt.FlowVec;
    if length(FlowVecPrev)~=length(FlowVec);
        FlowVecPrev= 0*FlowVec;
    end
else
    FlowVecPrev= 0*FlowVec;
end
% flow instationarity error
Errq= std(FlowVec-FlowVecPrev)/P.General.q0; % instationarity
iqSy= strmatch('SyVenRa',P.Valve.Name); % index systemic venous flow
ErrqSy= log(FlowVec(iqSy)/P.General.q0); %relative flow error
ErrFlow=sqrt(Errq.^2+ErrqSy.^2);
% display of normalized flow
P.Adapt.FlowVec= FlowVec;
str=[];
for i=1:6
    str= [str,' ',P.Valve.Name{i}];
end
disp(['Flow/q0: ',str]);
disp(FlowVec(1:6)/P.General.q0);
disp(['Relative FlowError 1000x= ',num2str(round(1000*ErrFlow))]);
% ---- End flow error calculation

if (Errq<0.1); % Adaptation only with sufficient stationarity

    % === recording stationarity variables before and after heart beat ====
    P.Adapt.In=[P.Adapt.In; P2StatVec('In')]; %record 1st sample of beat
    
    % ===== Pressure and flow control
    pMean= mean(GetFt('Cavity','p','SyArt'));
    qMean= mean(GetFt('Valve','q','SyVenRa'));
    Fac=0.5; %Feedback factor blood pressure control
    FacpControl= (pMean/P.General.p0)^Fac; % used for volume depletion
    P.General.FacpControl=FacpControl;
    p0AV= GetFt('ArtVen','p0AV','Sy'); % Peripheral resistance control
    PutFt('ArtVen','p0AV','Sy',(p0AV/FacpControl)*(qMean/P.General.q0)^Fac);

    % Estimate AV-delay
    P.General.TauAv=0.1765*P.General.tCycle;
  
    %=== Adapt ArtVen diameters
    ArtVenAdapt('All','Diameter');

    if 1;% Adapt diameter of valves to the connected blood vessel
        ValveNames={'SyVenRa','RvPuArt','PuVenLa','LvSyArt'};
        CavityNames={'SyVen','PuArt','PuVen','SyArt'};
        S=mean(GetFt('Cavity','A',CavityNames));
        PutFt('Valve','AOpen',ValveNames,S);
        % mitral and tricuspid valve are larger
        PutFt('Valve','AOpen',{'RaRv','LaLv'},1.5*S([2,4]));
    end
    P.Adapt.Out=[P.Adapt.Out; P2StatVec('Out') ]; %record last sample
    % for fast steady state procedure
end

% ===== end recording stationarity variables

% === if Error==small, ending is made faster
if size(P.Adapt.Out,1)>1; % Escape if steady state is reached
    
    ErrVec= log( P.Adapt.Out(end,:)./P.Adapt.In(end,:) );
    disp(['Stationarity error= ',...
        num2str(round(1000*sqrt(mean(ErrVec .^2 ))))] );
    if sqrt(mean(ErrVec .^ 2 ))<0.001
        P.General.tEnd= P.General.tEnd-0.5*(P.General.tEnd-P.t(end));
    end
    %=== ERROR criterium on flow stationarity
    
    % === Faster Steady State
    if P.Adapt.Fast;
        Vec= SteadyStateP;
        StatVec2P(Vec); % sets new initial conditions to StrInOut-variables
    end
    
end;

% After adaptation the initial condition is set for next beat.
P2SVar; % load physiologic data in record of state variables P.SVar
%SVarDot(0,P.SVar(end,:)',[]); % Most compact P-structure, contains starting
%conditions only.

end

%===== Stationarity functions for AdaptRest =====
function Vec= P2StatVec(InOut) % vector of variables used for reaching stationarity
global P;
% VecV=GetFt('Cavity','V','All');
VecV=GetFt('Cavity','V',{'SyArt','SyVen','PuArt','PuVen','La','Ra'});
VecA=GetFt('ArtVen','A0',{'Sy','Pu'}); VecA=VecA(:)';
Facp=P.General.FacpControl;
p0AVSy=GetFt('ArtVen','p0AV','Sy');
if strmatch(InOut,'In')
    Vec=[VecV(  1,:),VecA,Facp,p0AVSy];
else
    Vec=[VecV(end,:),VecA,Facp,p0AVSy];
end
end

function StatVec2P(Vec) % Inverse of P2StatPar
% fills values at the end of the arrays of state variables
global P;
PutFt('Cavity','V',{'SyArt','SyVen','PuArt','PuVen','La','Ra'},...
    Vec(1:6)); % Volumes used getting stationarity
A=Vec(6+[1,3;2,4]);
PutFt('ArtVen','A0',{'Sy','Pu'},A);
P.General.FacpControl=Vec(end-1);
PutFt('ArtVen','p0AV','Sy',Vec(end));
end

