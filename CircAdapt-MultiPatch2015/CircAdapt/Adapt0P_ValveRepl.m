function Adapt0P_ValveRepl
%function Adapt0P;

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
disp(num2str([FlowVec(1:6)/P.General.q0],'%11.4f'));
disp(['Relative FlowError 1000x= ',num2str(round(1000*ErrFlow))]);
% ---- End flow error calculation

if (Errq<0.4); % Adaptation only with sufficient stationarity

    % === recording stationarity variables before and after heart beat ====
    P.Adapt.In=[P.Adapt.In; P2StatVec('In')]; %record 1st sample of beat

    % ===== Pressure and flow control
    pMean= mean(GetFt('Cavity','p','SyArt'));
    qMean= mean(GetFt('Valve','q','SyVenRa'));
    Fac=0.5; %Feedback factor blood pressure control
    FacpControl= (pMean/P.General.p0)^Fac; % used for volume depletion     
    p0AV= GetFt('ArtVen','p0AV','Sy'); % Peripheral resistance control
    
    if P.General.PressFlowContr;%++++ yes/no pressure and flow control +++++++++++++++
        PutFt('ArtVen','p0AV','Sy',(p0AV/FacpControl)*(qMean/P.General.q0)^Fac);
        P.General.FacpControl=FacpControl;
    end

    % Estimate AV-delay
    P.General.TauAv=0.1765*P.General.tCycle;
    
    P.Adapt.Out=[P.Adapt.Out; P2StatVec('Out') ]; %record last sample
    % for fast steady state procedure
end

if P.Adapt.Fast;
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
        
    end
end

if isfield(P,'ValveRepl')
    P.ValveRepl.Count = P.ValveRepl.Count +1;
    ValveName = P.ValveRepl.ValveName;
    dALeak = P.ValveRepl.dALeak;
    dBeat = P.ValveRepl.dBeat;
    Count = P.ValveRepl.Count;
    if rem(Count-1,dBeat)==0
        PutFt('Valve','ALeak',ValveName,max(GetFt('Valve','ALeak',ValveName)-dALeak,2.6471e-10));
    end
end;

% get the initial condition for next beat, P is most compact information to
% start the simulation
P2SVar; % load physiologic data in record of state variables P.SVar
%SVarDot(0,P.SVar(end,:)',[]); % Compact P-structure

end

%===== Stationarity functions for Adapt0P =====
function Vec= P2StatVec(InOut) % vector of variables used for reaching stationarity
global P;
% VecV=GetFt('Cavity','V','All');
VecV=GetFt('Cavity','V',{'SyArt','SyVen','PuArt','PuVen','La','Ra'});
Facp=P.General.FacpControl;
p0AVSy=GetFt('ArtVen','p0AV','Sy');
if strmatch(InOut,'In')
    Vec=[VecV(1,:),Facp,p0AVSy];
else
    Vec=[VecV(end,:),Facp,p0AVSy];
end
end

function StatVec2P(Vec) % Inverse of P2StatPar
% fills values at the end
global P;
PutFt('Cavity','V',{'SyArt','SyVen','PuArt','PuVen','La','Ra'},Vec(1:6));
P.General.FacpControl=Vec(end-1);
PutFt('ArtVen','p0AV','Sy',Vec(end));
end

