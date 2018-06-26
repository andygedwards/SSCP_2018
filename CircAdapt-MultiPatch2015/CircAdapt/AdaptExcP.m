function AdaptExcP

global P
save PTemp P; %saves last intermediate solution

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
    StatVecIn=P2StatVec('In');
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
    
    if ErrFlow<0.1
        % === recording stationarity variables before and after heart beat ====
        P.Adapt.In=[P.Adapt.In; StatVecIn]; %record 1st sample of beat
       
        %=== Adapt ArtVen wall thickness and Patches
        ArtVenAdapt('All','WallVolume');
        %PatchAdapt('All',{'WallArea','EcmStress'});
        PatchAdapt('All','All');
        BagAdapt;
        
        P.Adapt.Out=[P.Adapt.Out; P2StatVec('Out') ]; %record last sample
    end

end % ---- end adaptation

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

% get the initial condition for next beat, P is most compact information to
% start the simulation
P2SVar; % load physiologic data in record of state variables P.SVar
%SVarDot(0,P.SVar(end,:)',[]); % Compact P-structure

end

%===== Stationarity functions for AdaptRest =====
function Vec= P2StatVec(InOut) % vector of variables used for reaching stationarity
global P;

i= strmatch(InOut,{'In','Out'}); A=[1,size(P.t,1)]; it=A(i);
VecA=GetFt('ArtVen','AWall','All'); VecA=VecA(:)';
VecSfPas=GetFt('Patch','SfPas','All');
VecAmRef=GetFt('Patch','AmRef','All');
VecVWall=GetFt('Patch','VWall','All');
% Adapt0P variables:
Aux=GetFt('Cavity','V',{'SyArt','SyVen','PuArt','PuVen','La','Ra'});
VecV=Aux(it,:);
Facp=P.General.FacpControl;
p0AVSy=GetFt('ArtVen','p0AV','Sy');
Vec=[VecA,VecSfPas,VecAmRef,VecVWall,VecV,Facp,p0AVSy];
end

function StatVec2P(Vec) % Inverse of P2StatPar
% fills values at the end
global P;
nP=P.Patch.n; nCV=6; nAV=P.ArtVen.n;
i=cumsum([0,2*nAV,nP,nP,nP,nCV,1,1]);
j= [i(1:end-1)+1; i(2:end)];
VRg=@(i) Vec(j(1,i):j(2,i)); %successive index ranges
PutFt( 'ArtVen','AWall','All',reshape(VRg(1),[2,nAV]));
PutFt( 'Patch','SfPas','All',VRg(2) );
PutFt( 'Patch','AmRef','All',VRg(3) );
PutFt( 'Patch','VWall','All',VRg(4) );
% Adapt0P variables:
PutFt( 'Cavity','V',{'SyArt','SyVen','PuArt','PuVen','La','Ra'},VRg(5) );
P.General.FacpControl=VRg(6);
PutFt('ArtVen','p0AV','Sy',VRg(7));

end

function BagAdapt % pericardial Bag adaptation to pAdapt reference
global P
pMax = max(P.Bag.p);
Fac_pMax= P.Bag.pAdapt./pMax;
P.Bag.VRef= P.Bag.VRef/(Fac_pMax.^(0.3/P.Bag.k));
disp(['Pericard adaptation: ',num2str(log(Fac_pMax))]);
end

