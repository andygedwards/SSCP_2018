function ChamberV2p
% function ChamberV2p
% A chamber is a cavity encapsulated in a single myocardial wall.
% Calculates: Cavity volume V ->
% myofiber stress Sf, wall tension T and cavity pressure p
% using the locally linearized T(Am) relation
% Theo Arts, Maastricht University, Oct 30, 2011


global P;
Chamber=P.Chamber; % Chamber structure

rhob    = 1050;
iCavity = Chamber.iCavity; % index of related cavity
iWall   = Chamber.iWall  ; % index of related wall
V    = max(0,P.Cavity.V(:,iCavity)); %cavity volumes
VWall= P.Wall.VWall(iWall); % wall volumes
nt   = size(V,1); % nt=number of time samples,
Vm   = max(0,V) + repmat(VWall/2,[nt,1]);% numerical safety with V<0 
Cm   = (4/3*pi./Vm).^(1/3); % mid wall curvature
Am   = (4*pi)./Cm.^2; % midwall area
Am0  = P.Wall.Am0(:,iWall); % zero tension midwall area
DADT = P.Wall.DADT(:,iWall); % wall compliance

Am=max(Am,Am0);% buckling with T<0

T    = (Am-Am0)./DADT; % wall tension
pTrans= 2*Cm.*T; % transmural pressure

% wall properties
P.Wall.T(:,iWall)     = T; % wall tension
P.Wall.Cm(:,iWall)    = Cm; % curvature=1/radius
P.Wall.Am(:,iWall)    = Am; % wall area
P.Wall.pTrans(:,iWall)= pTrans; % transmural pressure

 % Cavity impedance properties, needed to make node connection
Len= 2*Vm.^(1/3); % cavity length
A  = ( V + 0.1*repmat(VWall,[nt,1]) ) ./Len; % cross-sectional area
P.Cavity.A(:,iCavity) = A; % cross-sectional area for valve 
P.Cavity.Z(:,iCavity) = 0.2*sqrt(rhob*Len./abs(DADT))./A;
%         wave resistance to node

end

