function PatchAdapt(StrPatch,AdaptType)
% function PatchAdapt(StrPatch,AdaptType);
% StrPatch= array of Patch names, e.g. {'Lv1','Sv1'}
% AdaptType= {'WallVolume','WallArea','EcmStress'} indicates type of adaptation

global P

% Determine Patch and AdaptType indices
iPatch= Str2Index(StrPatch ,P.Patch.Name);
iAdapt= Str2Index(AdaptType,{'WallVolume','WallArea','EcmStress'});

Adapt=[0,0,0]         ; Adapt(iAdapt)=1;
AdaptWallVolume=Adapt(1);
AdaptWallArea  =Adapt(2);
AdaptEcmStress =Adapt(3);

Lsi     = P.Patch.Lsi; %sarcomere length
SfEcm   = max(0,P.Patch.SfEcm); % ECM fiber stress
SfTit   = max(0,P.Patch.SfPasT-P.Patch.SfEcm); % Myocyte passive stress
SfAct   = max(0,P.Patch.Sf-P.Patch.SfPasT); % Myocyte active stress
SfMyo   = max(0,P.Patch.Sf-P.Patch.SfEcm); % Total myocyte stress
Ef      = P.Patch.Ef; % tissue fiber strain
EfDot   = P.Patch.LsiDot./Lsi; % ~fiber strain rate
CDot    = P.Patch.CDot/max(P.Patch.C); % ~rate of contraction activation
LSe     = P.Patch.Ls-Lsi; % length series elastic element
 
%Target values
SfEcmMaxT= P.Patch.Adapt.SfPasMax(iPatch); % Ecm-stress
SfActMaxT= P.Patch.Adapt.FacSfAct(iPatch) .* P.Patch.SfAct(iPatch); %active stress
SfPasActT= P.Patch.Adapt.SfPasAct(iPatch); % ~integrin stress
% SfTitActT= (0.0165*Atrium+0.0332*Ventricle)*Sarc.SfAct;
SfTitActT= 0.05*P.Patch.Adapt.FacSfAct(iPatch).*P.Patch.SfAct(iPatch);
LsPasActT= 2.10;% setting integrin stress weighted sarcomere length
 
%sensed variables, ratio to target value
SfEcmMax= max(SfEcm)./SfEcmMaxT;
SfActMax= max(SfAct)./SfActMaxT;
SfPasAct= max(SfEcm .* SfAct./(max(100,SfEcm+SfAct)))./SfPasActT;
SfTitAct= max(SfTit .* SfAct./(SfTit+SfAct))./SfTitActT;
LsPasAct= sum(SfEcm.*SfAct.*Lsi)./sum(SfEcm.*SfAct)./LsPasActT;
% factor inversion for convenient mathematical manipulation
FacSfEcm   = 1./SfEcmMax; %1 
FacSfAct   = 1./SfActMax; %2
FacSfPasAct= 1./SfPasAct; %3
FacSfTitAct= 1./SfTitAct; %4
FacLsPasAct= 1./LsPasAct; %6

Eff2Geom=[3,0,0; 0,1,5; 0,0,-5]; % general geometry <-> local remodeling

% Calculate adaptation factors
switch 1 % various trials,
    % otherwise no change
    case 1 %1236 good correlation, best convergence
        Sens=log([FacSfEcm;FacSfAct;FacSfPasAct;FacLsPasAct]); % sensed variables
        Sens2Eff=[...
            -0.1917   -0.1279   -0.0906
            -0.2617    0.0999    0.1401
            +0.2021    0.2066   -0.0573
            -0.1855   -1.8441    0.5308];
        Aux= Eff2Geom'*Sens2Eff'*Sens;
    case 2 %1346 good correlation, better convergence
        Sens=log([FacSfEcm;FacSfPasAct;FacSfTitAct;FacLsPasAct]); % sensed variables
        Sens2Eff=[...
            -0.0170   -0.2198   -0.2201
            -0.0225    0.3120    0.0906
            -0.7362    0.4042    0.5689
            +5.0679   -4.8282   -3.6713];
        Aux= Eff2Geom'*Sens2Eff'*Sens;
    otherwise %no action

end

% Clipping of Fac around 1 with range +/-Clip 
a= 0.10; % gain of adaptation feedback. 1.0 represents best fit matrix%+++++++++++++++
ClipFac= @(x,Clip) exp(Clip*tanh(log(x)/Clip));
Clip=0.10;
FacVWall  = ClipFac(exp(a*Aux(1,:)),Clip);
FacAmRef  = ClipFac(exp(a*Aux(2,:)),Clip);
FacSfPas  = ClipFac(exp(a*Aux(3,:)),Clip);

%=== Fixation Septal/Lv wall area (adaptation is critical)+++++++++
% AmRefLS=[sum(GetFt('Patch','AmRef','Lv')) sum(GetFt('Patch','AmRef','Sv'))];
% LvSv=2.0;
% Fac= LvSv*AmRefLS(2)/AmRefLS(1);
% FacAmRef([3,4])= FacAmRef([3,4]).*[Fac,1/Fac].^0.5;
% a=0.1; b=0.6; mat= ( (1-a)*eye(2)+a*[1;1]*[b,1-b] );
% PutFt('Patch','AmRef',{'Lv1','Sv1'},AmRefLS * mat);
%---


%=== Carrying out adaptation

Error=[];

if AdaptWallVolume
    Vw0= P.Patch.VWall;
    Vw1= Vw0 .* FacVWall;
    P.Patch.VWall=Vw1;
    Error= [Error;abs(log(FacVWall))];
end

if AdaptWallArea
    %Strain softening of passive matrix by creep
%     Sheet.AmRef= Sheet.AmRef*FacAmRef;
    P.Patch.AmRef= P.Patch.AmRef .* FacAmRef;
    Error= [Error;abs(log(FacAmRef))];
end

if AdaptEcmStress 
    % Shifting active sarcomere range over passive matrix
    % Passive stress increase with same sarc length is equivalent with
    % decrease of sarcomere length with same passive stress
    P.Patch.SfPas= P.Patch.SfPas .* FacSfPas;
    Error= [Error; abs(log(FacSfPas))];
    % Softening by changing dLsPas seems very OK
    %     Sarc.dLsPas= Sarc.dLsPas * FacSfPas.^-0.3;
end

ErrorT= sqrt(mean(Error.^2,2));
str=[];
for i=1:P.Patch.n
    str=[str,' ',P.Patch.Name{i}];
end
disp(str);
disp(['VWall ',num2str(round(1000*[log(FacVWall),ErrorT(1)]),...
    '%+5.3d')]);
disp(['AmRef ',num2str(round(1000*[log(FacAmRef),ErrorT(2)]),...
    '%+5.3d')]);
disp(['SfPas ',num2str(round(1000*[log(FacSfPas),ErrorT(3)]),...
    '%+5.3d')]);

end


function Index=Str2Index(Str,Names)
n=size(Names,2);
Index=[];
if ischar(Str) %if single string
    if strcmp(Str,'All')
        Index=1:n;
    else
        Index=strmatch(Str,Names);
    end
else % if array of strings
    for i=1:length(Str)
        Index= [Index,strmatch(Str{i},Names)];
    end
end
end

