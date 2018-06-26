function PatchWallA2T
% function PatchWallA2T
% State variables of the sarcomere -> Patch area Am is written as a
% linear function of wall tension T: Am(T)=Am0+T*DADT
% Theo Arts, Maastricht University, Oct 30, 2011

global P;

% Patch Am= Am0+DADT*T, provides Am0 and DADT
AmRef=SparseDiag(P.Patch.AmRef); % midwall area for ref sarcomere length 2mu
LsRef=SparseDiag(P.Patch.LsRef); % ref sarcomere length 2mu
VWall=SparseDiag(P.Patch.VWall); % wall volume

Lsi   = P.Patch.Lsi; % unloaded sarc length= state variable sarc length due to passive stress only. No series elastic length.
Lambda= Lsi/LsRef; %extension factor
Am    = Lambda.^2*AmRef; %actual midwall area %area resulting from Ls = lsi only. Zero length series elastic element.
Ef    = log(Lambda); %natural fiber strain with series elastic element length = 0.

P.Patch.Ef= Ef; % set strain in patch module. % based on sarcomere length Lsi.
SarcEf2Sf; % sarcomere strain->stress % calculate the stiffness of all patches based on Ls = Lsi.
Sf   = P.Patch.Sf; % sarcomere stress % sarcomere passive stress at length Ls = Lsi.
T    = (Sf *(0.5 *VWall)) ./ Am   ; % tension % wall tension at length Ls = Lsi.
DADT = (Am.^2 ./ max(P.Patch.DSfDEf-2*Sf,1e3))/(0.25*VWall); %compliance %at length Ls = Lsi.
P.Patch.T   = T; % not used!
P.Patch.DADT= DADT;
P.Patch.Am0 = Am - T.*DADT; % zero load midwall area % true unloaded patch area with effect of passive stress removed.

% Wall is composed of patches: Also for wall: Am(T)=Am0+DADT*T
for iWall=1:P.Wall.n
    iPatch= (P.Wall.iPatch(iWall)-1)+(1:P.Wall.nPatch(iWall));
    P.Wall.VWall(iWall)  = sum(P.Patch.VWall(iPatch));
    P.Wall.Am0(:,iWall)  = P.Wall.AmDead(iWall)+sum(P.Patch.Am0(:,iPatch),2);
    P.Wall.DADT(:,iWall) = sum(P.Patch.DADT(:,iPatch),2);
end
end

