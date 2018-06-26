function PNew
% based on ParRef
% should be made stand-alone
% Theo Arts, Maastricht University, Oct 30, 2011

global P;P=[];


% load ParRef;
P.General.rhob=1050;
P.General.q0=45e-6;
P.General.p0=12000;
P.General.tCycle=0.600;
P.General.FacpControl=1;
P.General.dTauAv=0;
pLv=19000;
pRv=7500;
pLa=2800;
pRa=1400;

%=============================
qExc=3*P.General.q0;
VStroke= P.General.q0*P.General.tCycle;
VStrokeExc=1.5*VStroke;


%===== INPUT: naming ArtVen, TriSeg and Chambers
P.ArtVen.Name = {'Sy','Pu'};
P.Chamber.Name= {'La','Ra'};
P.TriSeg.Name = {'v'};
%===============================================

jCavity=0; jWall=0; jPatch=0; jNode=0;
P.Cavity.Name= {}; P.Cavity.n= 0;
P.Wall.Name  = {}; P.Wall.n  = 0; P.Wall.nPatch=[];
P.Patch.Name = {}; P.Patch.n = 0;
P.Node.Name  = {}; P.Node.n  = 0;

% ===== ArtVen -> Cavity
P.ArtVen.n = length(P.ArtVen.Name);
str={'Art','Ven'};
NameCav={}; iCav=[];
NameWall={}; iWall=[];
for i=1:P.ArtVen.n
    str1=P.ArtVen.Name{i};
    iCav=[iCav,jCavity+1];
    iWall=[iWall,jWall+1];
    for j=1:2
        jCavity=jCavity+1;
        NameCav=[NameCav,[str1,str{j}]];
        jWall=jWall+1;
        NameWall=[NameWall,[str1,str{j}]];
    end
end
P.ArtVen.iCavity=iCav;
P.Cavity.Name= [P.Cavity.Name,NameCav];
P.Cavity.n   = P.Cavity.n+length(NameCav);
P.ArtVen.iWall=iWall;
P.Wall.Name= [P.Wall.Name,NameWall];
P.Wall.n   = jWall;
P.Wall.nPatch=[P.Wall.nPatch,zeros(size(NameWall))];

% ===== Chamber -> Cavity, Wall
P.Chamber.n = length(P.Chamber.Name);
NameCav ={}; iCav =[];
NameWall={}; iWall=[];
for i=1:P.Chamber.n
    str1=P.Chamber.Name{i};
    jCavity=jCavity+1;
    iCav=[iCav,jCavity];
    jWall=jWall+1;
    iWall=[iWall,jWall];
    NameCav = [NameCav,str1];
    NameWall= [NameWall,str1];
end
P.Chamber.iCavity= iCav;
P.Cavity.Name    = [P.Cavity.Name,NameCav];
P.Cavity.n       = P.Cavity.n+length(NameCav);
P.Chamber.iWall  = iWall;
P.Wall.Name      = [P.Wall.Name,NameWall];
P.Wall.n         = jWall;
P.Wall.nPatch    =[P.Wall.nPatch,ones(size(NameWall))];

% ===== TriSeg-> Cavity, Wall
P.TriSeg.n = length(P.TriSeg.Name);
strC={'L','R'};
strW={'L','S','R'};
NameCav ={}; iCav= [];
NameWall={}; iWall=[];
for i=1:P.TriSeg.n
    str1=P.TriSeg.Name{i};
    iCav =[iCav,jCavity+1];
    iWall=[iWall,jWall+1];
    for j=1:2
        jCavity=jCavity+1;
        NameCav=[NameCav,[strC{j},str1]];
    end
    for j=1:3
        jWall=jWall+1;
        NameWall=[NameWall,[strW{j},str1]];
    end
end
P.TriSeg.iCavity= iCav;
P.Cavity.Name   = [P.Cavity.Name,NameCav];
P.Cavity.n      = P.Cavity.n+length(NameCav);
P.TriSeg.iWall  = iWall;
P.Wall.Name     = [P.Wall.Name,NameWall];
P.Wall.n        = jWall;
P.Wall.nPatch   =[P.Wall.nPatch,ones(size(NameWall))];

% === Bag = passively elastic bag like pericardium
P.Bag.Name={'Peri'};
P.Bag.n= size(P.Bag.Name,1);
P.Bag.iWall= jWall+(1:P.Bag.n);
jWall= jWall+P.Bag.n;
P.Wall.Name=[P.Wall.Name,P.Bag.Name];
P.Wall.n   =jWall;
P.Wall.nPatch=[P.Wall.nPatch,zeros(1,P.Bag.n)];

% Cavity-> Nodes (Default)
NameNode ={}; iNode =[];
for i=1:P.Cavity.n
    str1=P.Cavity.Name{i};
    jNode=jNode+1;
    iNode=[iNode,jNode];
    NameNode= [NameNode,str1];
end
P.Cavity.iNode  = iNode;
P.Node.iCavity  = iNode; % transfer of cavity properties to related node
P.Node.Name     = [P.Node.Name,NameNode];
P.Node.n        = P.Node.n+length(NameNode);
        
% ===== INPUT: VALVE with Node connections
ProxDist=[...
{'SyVen','Ra'   };...
{'Ra'   ,'Rv'   };...
{'Rv'   ,'PuArt'};...
{'PuVen','La'   };...
{'La'   ,'Lv'   };...
{'Lv'   ,'SyArt'};...
{'La'   ,'Ra'   };...
{'Lv'   ,'Rv'   };...
{'SyArt','PuArt'}];
nValve= size(ProxDist,1);
NameValve ={};
iProx=[]; iDist=[];
for i=1:nValve
    StrProx=ProxDist{i,1};
    StrDist=ProxDist{i,2};
    iP=strmatch(StrProx,P.Node.Name);
    iD=strmatch(StrDist,P.Node.Name);
    iProx=[iProx,iP]; iDist=[iDist,iD];
    Str=[StrProx,StrDist];
    NameValve=[NameValve,Str]; 
end
P.Valve.Name      = NameValve;
P.Valve.n         = nValve;
P.Valve.iNodeProx = iProx;
P.Valve.iNodeDist = iDist;

% ===== INPUT number of Patches per Wall

P.Patch.n=sum(P.Wall.nPatch);
iPatch=[1,1+cumsum(P.Wall.nPatch(1:end-1))];
P.Wall.iPatch=iPatch;%link to patch
NamePatch={};
for i=1:P.Wall.n
    n=P.Wall.nPatch(i);
    for j=1:n
        jPatch=jPatch+1;
        NamePatch=[NamePatch,[P.Wall.Name{i},num2str(j)]];
    end
end
P.Patch.Name= NamePatch;
P.Patch.dT=zeros(1,P.Patch.n);

%============= end backbone definition of structure =========


PutFt('Patch','Lsi'             ,'All',1.9);
PutFt('Patch','C'               ,'All',0.001);
PutFt('Patch','ActivationDelay' ,'All',[0.0;0.0]);
PutFt('Patch','LsRef'           ,'All',2.0);
PutFt('Patch','Ls0Pas'          ,'All',1.8);
PutFt('Patch','dLsPas'          ,'All',0.6);
PutFt('Patch','SfPas'           ,'All',[50,50,22,22,22]*1000);%+++++++++++
PutFt('Patch','Lsi0Act'         ,'All',1.51);
PutFt('Patch','LenSeriesElement','All',0.04);
PutFt('Patch','SfAct'           ,'All',[84,84,120,120,120]*1e3);
PutFt('Patch','vMax'            ,'All',[14,14,7,7,7]);
PutFt('Patch','TimeAct'         ,'All',[.15,.15,.42,.42,.42]);
PutFt('Patch','TR'              ,'All',[.4,.4,.25,.25,.25]);
PutFt('Patch','TD'              ,'All',[.4,.4,.25,.25,.25]);
PutFt('Patch','CRest'           ,'All',0.0);

% VWall=VStrokeExc*[3.3,3.9,5.8,2.4,7.4].*[pLa,pRa,pLv,pLv,pRv]./P.Patch.SfAct;
% AmRef=(VStrokeExc*[4.1,4.1,9.6,2.3,12.3]).^(2/3);
VWall=VStrokeExc*[3.3,3.9,5.8,2.4,7.4].*[pLa,pRa,pLv,pLv,pRv]./P.Patch.SfAct;
AmRef=(VStrokeExc*[4.1,4.1,9.0,3.0,12.0]).^(2/3);
PutFt('Patch','VWall','All',VWall);
PutFt('Patch','AmRef','All',AmRef);


%======= Adaptation setpoints

% ArtVen.Adapt
Adapt=[];
Mat=zeros(2,P.ArtVen.n);
Adapt.WallStress=Mat+500e3;
Adapt.vFlowMean =Mat+0.17;
Adapt.vImpact   =Mat+3;
P.ArtVen.Adapt  =Adapt;

Adapt=[];
Vec=zeros(1,P.Patch.n);
Adapt.LsBe     =Vec+2.3000;
Adapt.LsEe     =Vec+1.7500;
Adapt.SfPasMax =Vec;
iPatch=Str2iD({'La','Ra','Lv','Sv','Rv'},P.Patch.Name);%++++
Adapt.SfPasMax(iPatch)=[60000,60000,6000,6000,6000];
Adapt.SfPasAct(iPatch)=[6800,6800,4800,4800,4800];
Adapt.FacSfAct(iPatch)=[0.35,0.35,0.61,0.61,0.61];
P.Patch.Adapt=Adapt;
%==============

%===================

iC=P.ArtVen.iCavity; iCav=[iC;iC+1];
PutFt('ArtVen','k','All',[8;10]);
PutFt('ArtVen','Len','All',[0.4,0.2;0.4,0.2]);
PutFt('ArtVen','A0','All',P.General.q0./P.ArtVen.Adapt.vFlowMean);
PutFt('ArtVen','p0','All',[[12;0.12],[1.8;0.5]]*1e3);
PutFt('ArtVen','AWall','All',[0.23,0.18; 0.09,0.11].*P.ArtVen.A0);

PutFt('ArtVen','p0AV',{'Sy','Pu'},[P.General.p0,1500]);
PutFt('ArtVen','q0AV','All',P.General.q0);
PutFt('ArtVen','kAV',{'Sy','Pu'},[1.0,2.0]);
%++++++++++++++++++++++
% ArtVen specials
% Cavity: starting volumes
HCav={'La','Ra','Lv','Rv'};
nCav=P.Cavity.n;
P.Cavity.V=zeros(1,nCav);
PutFt('Cavity','V',HCav,VStrokeExc*[0.5,0.41,1.00,0.72]);
iCavArt=P.ArtVen.iCavity;
iCavVen=iCavArt+1;
V=P.ArtVen.A0 .* P.ArtVen.Len;
P.Cavity.V(iCavArt)=V(1,:);
P.Cavity.V(iCavVen)=V(2,:);

%TriSeg
PutFt('TriSeg','V','v',0.52*VStroke);
PutFt('TriSeg','Y','v',0.80*VStroke^(1/3));
P.TriSeg.Tau = 0.005;

% Peri
VWall=sum(GetFt('Patch','VWall',{'La1','Ra1','Lv1','Sv1','Rv1'}));
VCav=sum(GetFt('Cavity','V',{'La','Ra','Lv','Rv'}));
P.Bag.VRef  = (67/53)*(VCav+VWall);
P.Bag.k     = 10.0;

A0=GetFt('ArtVen','A0','Sy'); A=A0(1);
PutFt('Valve','q','All',0.0);
PutFt('Valve','AOpen','All',A);
PutFt('Valve','ALeak','All',A*1e-6);
PutFt('Valve','Len','All',sqrt(A));
PutFt('Valve','AOpen',{'RaRv','LaLv'},...
    1.5* GetFt('Valve','AOpen',{'RaRv','LaLv'}) );
PutFt('Valve','ALeak',{'SyVenRa','PuVenLa'},...
    GetFt('Valve','AOpen',{'SyVenRa','PuVenLa'}) );
PutFt('Valve','AOpen',{'LaRa','LvRv','SyArtPuArt'},...
    GetFt('Valve','ALeak',{'LaRa','LvRv','SyArtPuArt'}) );

% Wall: AmDead
ALv=sum(GetFt('Valve','AOpen',{'LvSyArt','LaLv'}));
ARv=sum(GetFt('Valve','AOpen',{'RvPuArt','RaRv'}));
ALa=sum(GetFt('Valve','AOpen',{'PuVenLa','LaLv'}));
ARa=sum(GetFt('Valve','AOpen',{'SyVenRa','RaRv'}));
PutFt('Wall','AmDead',{'Lv','Rv','Ra','La'},[ALv,ARv,ALa,ARa]); 

P.t=0;

%======= Adaptation setpoints
P.Bag.pAdapt= 100; % pericardial pressure

Adapt=[];
Vec=zeros(1,P.Cavity.n);
Adapt.WallStress=Vec;
Adapt.vFlowMean =Vec;
Adapt.vImpact   =Vec;
iCav=Str2iD({'SyArt','SyVen','PuArt','PuVen'},P.Cavity.Name);%++++
Adapt.WallStress(iCav)= 500e3;
Adapt.vFlowMean(iCav) = 0.17;
Adapt.vImpact(iCav)   = 3;
P.Cavity.Adapt=Adapt;

Adapt=[];
Vec=zeros(1,P.Patch.n);
Adapt.LsBe     =Vec+2.3000;
Adapt.LsEe     =Vec+1.7500;
Adapt.SfPasMax =Vec;
iPatch=Str2iD({'La','Ra','Lv','Sv','Rv'},P.Patch.Name);%++++
Adapt.SfPasMax(iPatch)=[60000,60000,6000,6000,6000];
Adapt.SfPasAct(iPatch)=[6800,6800,4800,4800,4800];
Adapt.FacSfAct(iPatch)=[0.35,0.35,0.61,0.61,0.61];
P.Patch.Adapt=Adapt;

%
%================= General information ====================
P.General.ScaleVqY=[1e-5,1e-4,1e-1];
P.General.Dt=0.002;
P.General.tEnd=36.0;
P.General.tCycleRest=0.8500;
P.General.tCycle= 0.8500;
P.General.TimeFac=1;
P.General.TauAv=0.16;
P.General.q0=85e-6;
P.General.p0=12200;
P.General.PressFlowContr=1;
end

function iD= Str2iD(Str,Name)
n=length(Str);
iD=zeros(1,n);
for i=1:n;
    iD(i)=strmatch(Str{i},Name);
end
end

