function SVar2P
% function SVar2P
% Transfer of scaled state variables in P.SVar to P-structure (SI-units)
% Theo Arts, Maastricht University, Oct 30, 2011

global P

ScaleVqY=P.General.ScaleVqY; % scaling factors for volume, flow, distance
ScV=ScaleVqY(1); Scq=ScaleVqY(2); ScY=ScaleVqY(3); Sc0=1;
nCav  =P.Cavity.n; % number of cavities
nValve=P.Valve.n; % number of valves
nPatch=P.Patch.n; % number of patches with representative sarcomere
% finding indices for value transfer
a=cumsum([0,1,nCav,nValve,nPatch,nPatch,1,1]);
iB=a(1:end-1)+1; iE=a(2:end); % successive begin and end indices

% value transfer
i=1;
P.t        =     P.SVar(:,iB(i):iE(i)); i=i+1;
P.Cavity.V = ScV*P.SVar(:,iB(i):iE(i)); i=i+1;
P.Valve.q  = Scq*P.SVar(:,iB(i):iE(i)); i=i+1;
P.Patch.C  =     P.SVar(:,iB(i):iE(i)); i=i+1;
P.Patch.Lsi=     P.SVar(:,iB(i):iE(i)); i=i+1;
P.TriSeg.V = ScV*P.SVar(:,iB(i):iE(i)); i=i+1;
P.TriSeg.Y = ScY*P.SVar(:,iB(i):iE(i)); i=i+1;

end    

