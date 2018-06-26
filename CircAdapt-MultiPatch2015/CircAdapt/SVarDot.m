function OutDot=SVarDot(tDummy,SVar,flag)
% function OutDot=SVarDot(tDummy,SVar,flag)
% State Variables -> Time derivatives of State Variables
% ready to use in MatLab function odeXXX() to solve set of Diff Eq
% Theo Arts, Maastricht University, Oct 30, 2011

%====test by execution of
%Par2P; P2SVar; SVDt=SVarDot(0,P.SVar',0);
%====

global P

P.SVar= SVar'; % store state variables SVar
SVar2P; % state variables SVar -> physiologic representation P.xx
P.tDot=ones(size(P.t)); % time derivative of time = 1
MemAlloc; % necessary memory allocations
ArtVenV2p; % arteries to veins compliant network
PatchWallA2T; % patch and wall: Am= Am0+T*DADT
ChamberV2p; % Chamber, pTrans and Wall.Am,T
TriSegV2p; % TriSeg, pTrans and Wall.Am,T
Wall2Patch2Sarc; % filling Sarc with LsiDot,CDot
BagV2p; % pTrans of Pericardium
pCavity; % Determines p in cavities by adding P.Wall.pTrans(:,xx)
pNodeVDot; % Node pressures and cavity VDot's
ValveqDot; % flow derivatives qDot
P2SVarDot; % transfer of derivatives to P-structure

OutDot= real(P.SVarDot)'; % odeXX requires SVarDot to be a row vector

end

