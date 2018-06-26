% Conversion of Par-structure to network of Chambers, tubes and valves with
% nodes.
clear;
load ParRef
global P;

% NameNode={};
% for i=1:P.Cavity.n
% NameNode=[NameNode,[NameCavity{i},'1']];
% end


%====== BACKBONE OF STRUCTURE ===============
jCavity=0; jWall=0; jPatch=0; jNode=0;
NameCavity={};
NameNode={};
NameWall={};
nPWall=[]; %number of myocardial patches per wall

%ARTVEN
Name={'Syst','Pulm'};
NameCav={'SystArt','SystVen','PulmArt','PulmVen'};
n=length(Name);
ArtVen.n=n;
ArtVen.Name=Name;
%Cavities
iCavity=jCavity+(1:2:2*n); jCavity=iCavity(end)+1;
ArtVen.iCavity=iCavity; %link to cavity
NameCavity=[NameCavity,NameCav];
%Walls
iWall=jWall+(1:2:2*n); jWall=iWall(end)+1;
ArtVen.iWall=iWall; %link to wall
NameWall=[NameWall,NameCav];
nPWall=[nPWall,zeros(1,2*n)];

%TRISEG
Name={'Ventricles'};
n=length(Name);
TriSeg.n=n;
TriSeg.Name=Name;
NameW={'Lv','Sv','Rv'};
%Cavities
iCavity=jCavity+(1:2:2*n); jCavity=iCavity(end)+1;
TriSeg.iCavity=iCavity;%link to cavity
NameCavity=[NameCavity,NameW([1,3])];
%Walls
iWall=jWall+(1:3:3*n); jWall=iWall(end)+2;
TriSeg.iWall=iWall;%link to wall
NameWall=[NameWall,NameW];
nPWall=[nPWall,ones(1,3*n)];

%CHAMBER
Name={'La','Ra'};
n=length(Name);
Chamber.n=n;
Chamber.Name=Name;
%Cavities
iCavity=jCavity+(1:n); jCavity=iCavity(end);
Chamber.iCavity=iCavity;%link to cavity
NameCavity=[NameCavity,Name];
%Walls
iWall=jWall+(1:n); jWall=iWall(end);
Chamber.iWall=iWall;%link to wall
NameWall  =[NameWall,Name];
nPWall=[nPWall,ones(1,n)];

%NODE
NameNode=NameCavity;

%PERICARDIUM
Name={'Peri'};
n=length(Name);
Peri.n=n;
Peri.Name=Name;
%Walls
iWall=jWall+(1:n); jWall=iWall(end);
Peri.iWall=iWall;%link to wall
NameWall  =[NameWall,Name];
nPWall=[nPWall,0];

%VALVE
Name={'PulmVen','LAv','SystArt','SystVen','RAv','PulmArt','ASD','VSD','DUCT'};
n=length(Name);
Valve.n=n;
Valve.Name=Name;

%===================
%Memory allocations in Cavity, Wall and Patch
%Cavities
Cavity.n=jCavity;
Cavity.Name=NameCavity;
Cavity.iNode=1:Cavity.n;
%Walls
Wall.n=jWall;
Wall.Name=NameWall;
Wall.nPatch=nPWall;%number of patches
%Patches
iPatch=[1,1+cumsum(Wall.nPatch(1:end-1))];
Wall.iPatch=iPatch;%link to patch
Patch.n=sum(nPWall);
NamePatch={};
for i=1:Wall.n
    n=Wall.nPatch(i);
    for j=1:n
        NamePatch=[NamePatch,[NameWall{i},num2str(j)]];
    end
end
Patch.Name=NamePatch;
%Nodes
Node.n=Cavity.n;
Node.Name= NameNode;

P=[];
P.t      = Par.t;
P.tDot   = ones(size(P.t));
P.General= Par.General;
P.ArtVen = ArtVen;
P.TriSeg = TriSeg;
P.Chamber= Chamber;
P.Peri   = Peri;
P.Cavity = Cavity;
P.Wall   = Wall;
P.Patch  = Patch;
P.Valve  = Valve;
P.Node   = Node;

%============= end backbone definition of structure =========

% ArtVen
nt=length(Par.TubeLArt.V);

% initialization matrices and vectors
nArtVen=P.ArtVen.n;
MatArtVen=zeros(nt,nArtVen);
VecArtVen=zeros(1,nArtVen);

nTriSeg=P.TriSeg.n;
MatTriSeg=zeros(nt,nTriSeg);
VecTriSeg=zeros(1,nTriSeg);

nChamber=P.Chamber.n;
MatChamber=zeros(nt,nChamber);
VecChamber=zeros(1,nChamber);

nCavity=P.Cavity.n;
MatCavity=zeros(nt,nCavity);
VecCavity=zeros(1,nCavity);

nWall=P.Wall.n;
MatWall=zeros(nt,nWall);
VecWall=zeros(1,nWall);

nPatch=P.Patch.n;
MatPatch=zeros(nt,nPatch);
VecPatch=zeros(1,nPatch);

MatPeri=zeros(nt,1);

nValve=P.Valve.n;
MatValve=zeros(nt,nValve);
VecValve=zeros(1,nValve);

nNode=P.Node.n;
MatNode=zeros(nt,nNode);
VecNode=zeros(1,nNode);

P.Peri.p       = MatPeri;

P.Cavity.V     = MatCavity;
P.Cavity.VDot  = MatCavity;
P.Cavity.p     = MatCavity;
P.Cavity.Z     = MatCavity;
P.Cavity.A     = MatCavity;
P.Wall.AmDead  = VecWall;
P.Wall.VWall   = VecWall;
P.Wall.Am      = MatWall;
P.Wall.Cm      = MatWall;
P.Wall.Am0     = MatWall;
P.Wall.DADT    = MatWall;
P.Wall.T       = MatWall;
P.Wall.pTrans  = MatWall;

P.Patch.Lsi             = MatPatch;
P.Patch.LsiDot          = MatPatch;
P.Patch.C               = MatPatch;
P.Patch.CDot            = MatPatch;
P.Patch.VWall           = VecPatch;
P.Patch.AmRef           = VecPatch;
%Sarcomere
P.Patch.ActivationDelay =[VecPatch;VecPatch];
P.Patch.Ef              = MatPatch;
P.Patch.LsRef           = VecPatch;
P.Patch.Ls0Pas          = VecPatch;
P.Patch.dLsPas          = VecPatch;
P.Patch.SfPas           = VecPatch;
P.Patch.Lsi0Act         = VecPatch;
P.Patch.LenSeriesElement= VecPatch;
P.Patch.SfAct           = VecPatch;
P.Patch.vMax            = VecPatch;
P.Patch.TimeAct         = VecPatch;
P.Patch.TR              = VecPatch;
P.Patch.TD              = VecPatch;
P.Patch.CRest           = VecPatch;
P.Patch.Ls              = MatPatch;
P.Patch.SfPasT          = MatPatch;
P.Patch.Sf              = MatPatch;
P.Patch.DSfDEf          = MatPatch;
P.Patch.Am              = MatPatch;
P.Patch.T               = MatPatch;
P.Patch.DADT            = MatPatch;

P.Valve.q    = MatValve;
P.Valve.qDot = MatValve;
P.Valve.AOpen= VecValve;
P.Valve.ALeak= VecValve;
P.Valve.Len  = VecValve;

P.Node.q = MatNode;
P.Node.p = MatNode;
P.Node.Y = MatNode;

% ==== Data transfer Par->P ====
% ArtVen and related cavities
for i=1:nArtVen
    VecMat=zeros(2,nArtVen);
    iCav0 =P.ArtVen.iCavity(i);
    iWall0=P.ArtVen.iWall(i);
    P.ArtVen.rhob=1050;
    NamePar={'LArt','RVen','RArt','LVen'};
    for j=0:1
        jCavity=iCav0+j;
        Name=NamePar{jCavity};
        Tube=Par.(['Tube',Name]);
        P.Cavity.V(:,jCavity)= Tube.V;
 
        jWall=iWall0+j;
        Name=P.Wall.Name{jWall};        

        P.ArtVen.k(j+1,i)    = Tube.k;
        P.ArtVen.Len(j+1,i)  = Tube.Len;
        P.ArtVen.AWall(j+1,i)= Tube.AWall;
        P.ArtVen.p0(j+1,i )  = Tube.p0;
        P.ArtVen.A0(j+1,i )  = Tube.A0;

        P.ArtVen.p0AV= [Par.LRp.R*P.General.q0,P.General.pDropPulm];
        P.ArtVen.q0AV= P.General.q0*[1,1];
        P.ArtVen.kAV= [1.0, 2.0];

    end
end

%TriSeg and related cavities, walls and patches
P.TriSeg.V   = Par.Sv.V;
P.TriSeg.VDot= ones(size(Par.Sv.V));
P.TriSeg.Y   = Par.Sv.Y;
P.TriSeg.YDot= ones(size(Par.Sv.Y));
P.TriSeg.Tau=Par.Sv.Tau;
for i=1:nTriSeg %each TriSeg
    % cavity: volumes
    iCav0 =P.TriSeg.iCavity(i);
    jWall0=P.TriSeg.iWall(i);
    for j=0:1
        jCavity=iCav0+j;
        Name=P.Cavity.Name{jCavity};
        Aux=Par.(Name);%read
        P.Cavity.V(:,jCavity)=Aux.V;
    end
    % patch: sarcomere Lsi and C
    for j=0:2 %each wall of TriSeg
        jWall=jWall0+j;
        nPatch= P.Wall.nPatch(jWall);
        RgPatch=P.Wall.iPatch(jWall)+(0:nPatch-1);
        Vec= ones(1,nPatch)/nPatch;
        Name=P.Wall.Name{jWall};
        Sarc=Par.(Name).Sarc;%read
        P.Wall.AmDead(jWall)  = Par.(Name).AmDead;
        P.Wall.VWall(jWall)   = Par.(Name).VWall;
        P.Patch.Lsi(:,RgPatch)= Sarc.Lsi;
        P.Patch.C(:,RgPatch)  = Sarc.C;
        P.Patch.VWall(RgPatch)= Par.(Name).VWall*Vec;
        P.Patch.AmRef(RgPatch)= Par.(Name).AmRef*Vec;

        P.Patch.ActivationDelay(:,RgPatch) =Sarc.ActivationDelay([-1,0]+end);
        P.Patch.LsRef(RgPatch)           = Sarc.LsRef;
        P.Patch.Ls0Pas(RgPatch)          = Sarc.Ls0Pas;
        P.Patch.dLsPas(RgPatch)          = Sarc.dLsPas;
        P.Patch.SfPas(RgPatch)           = Sarc.SfPas;
        P.Patch.Lsi0Act(RgPatch)         = Sarc.Lsi0Act;
        P.Patch.LenSeriesElement(RgPatch)= Sarc.LenSeriesElement;
        P.Patch.SfAct(RgPatch)           = Sarc.SfAct;
        P.Patch.vMax(RgPatch)            = Sarc.vMax;
        P.Patch.TimeAct(RgPatch)         = Sarc.TimeAct;
        P.Patch.TR(RgPatch)              = Sarc.TR;
        P.Patch.TD(RgPatch)              = Sarc.TD;
        P.Patch.CRest(RgPatch)           = Sarc.CRest;

    end
end

%Chamber and related cavities, walls and patches
for i=1:nChamber
    % cavity: volumes
    jCavity=P.Chamber.iCavity(i);
    Name=P.Cavity.Name{jCavity};
    Aux=Par.(Name);
    P.Cavity.V(:,jCavity)=Aux.V;

    % patch: sarcomere Lsi and C
    jWall  = P.Chamber.iWall(i);
    nPatch = P.Wall.nPatch(jWall);
    RgPatch= P.Wall.iPatch(jWall)+(0:nPatch-1);
    Vec    = ones(1,nPatch)/nPatch;
    Name   = P.Wall.Name{jWall};
    Sarc   = Par.(Name).Sarc;
    P.Wall.AmDead(jWall)  = Par.(Name).AmDead;
    P.Wall.VWall(jWall)   = Par.(Name).VWall;
    P.Patch.Lsi(:,RgPatch)= Sarc.Lsi;
    P.Patch.C(:,RgPatch)  = Sarc.C;
    P.Patch.VWall(RgPatch)= Par.(Name).VWall*Vec;
    P.Patch.AmRef(RgPatch)= Par.(Name).AmRef*Vec;

    P.Patch.ActivationDelay(:,RgPatch) =Sarc.ActivationDelay([-1,0]+end);
    P.Patch.LsRef(RgPatch)           = Sarc.LsRef;
    P.Patch.Ls0Pas(RgPatch)          = Sarc.Ls0Pas;
    P.Patch.dLsPas(RgPatch)          = Sarc.dLsPas;
    P.Patch.SfPas(RgPatch)           = Sarc.SfPas;
    P.Patch.Lsi0Act(RgPatch)         = Sarc.Lsi0Act;
    P.Patch.LenSeriesElement(RgPatch)= Sarc.LenSeriesElement;
    P.Patch.SfAct(RgPatch)           = Sarc.SfAct;
    P.Patch.vMax(RgPatch)            = Sarc.vMax;
    P.Patch.TimeAct(RgPatch)         = Sarc.TimeAct;
    P.Patch.TR(RgPatch)              = Sarc.TR;
    P.Patch.TD(RgPatch)              = Sarc.TD;
    P.Patch.CRest(RgPatch)           = Sarc.CRest;

end

% Cavity
P.Cavity.iNode= 1:P.Cavity.n;
P.Node.iCavity= 1:P.Cavity.n;% for now ++++++

% Peri
P.Peri.VRef  = Par.Peri.VRef  ;
P.Peri.k     = Par.Peri.k     ;
P.Peri.pAdapt= Par.Peri.pAdapt;

% Valve
ParName={'LVen','LAv','LArt','RVen','RAv','RArt','ASD','VSD','DUCT'};
for i=1:nValve
    Aux=Par.(['Valve',ParName{i}]);
    P.Valve.q(:,i)  =Aux.q;
    P.Valve.AOpen(i)=Aux.AOpen;
    P.Valve.ALeak(i)=Aux.ALeak(1,:);
    P.Valve.Len(i)  =Aux.Len;
end
P.Valve.iNodeProx=[4,7,5,2,8,6,7,5,1];
P.Valve.iNodeDist=[7,5,1,8,6,3,8,6,3];

%================= SVar communication ====================
P.General.ScaleVqY=[1e-5,1e-4,1e-1];%????++++++

