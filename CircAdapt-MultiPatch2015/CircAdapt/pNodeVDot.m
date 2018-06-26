function pNodeVDot
% function pNodeVDot
% Pressures in cavities p -> pressure in Nodes: p
%                         -> flows to cavities: VDot
% Theo Arts, Maastricht University, Oct 30, 2011

global P
P.Node.q=0*P.Node.q; P.Node.Y=P.Node.q;% zero initialization
iNodeProx=P.Valve.iNodeProx; % nodes proximal to valves
iNodeDist=P.Valve.iNodeDist; % nodes distal to valves
nValve=P.Valve.n; % number of valves

% flow of valves to/from nodes
for i=1:nValve
    iProx=iNodeProx(i);
    iDist=iNodeDist(i);
    q=P.Valve.q(:,i); % valve flow
    P.Node.q(:,iProx)=P.Node.q(:,iProx)-q; % valve flow out of nodes
    P.Node.q(:,iDist)=P.Node.q(:,iDist)+q; % valve flow into nodes
end

% flow of ArtVen object to/from nodes
p0AV= SparseDiag(P.ArtVen.p0AV); % reference AV pressure drop
q0AV= SparseDiag(P.ArtVen.q0AV); % reference AV flow
kAV = P.ArtVen.kAV; % exponent non-linearity AV 'resistance'
iCav= P.ArtVen.iCavity; % related cavities
iNodeP= P.Cavity.iNode(iCav); % nodes related to arterial cavity
iNodeD= P.Cavity.iNode(iCav+1); % nodes related to venous cavity
Dp= P.Cavity.p(:,iCav)-P.Cavity.p(:,iCav+1); % AV cavity pressure drop
Col1= ones(size(Dp,1),1);
q=(Dp/p0AV).^kAV(Col1,:)*q0AV; % AV flow
P.ArtVen.q=q;

%Blood pressure control
Facq=P.General.FacpControl; % FacpControl= current p / pRef
P.Node.q(:,iNodeP)=P.Node.q(:,iNodeP)-q*Facq; % flow out of node
P.Node.q(:,iNodeD)=P.Node.q(:,iNodeD)+q/Facq; % flow into node

% Blood pressure control

% node pressure derived for single node-cavity connection
iCavNode= P.Cavity.iNode; % iNode pointed by cavity
Y= 1./P.Cavity.Z; % should be a summation, not needed now
Yp= P.Cavity.p .* Y; % should be a summation, not needed now
P.Node.Y(:,iCavNode)=Y; % internal Y conductance of flow source
P.Node.q(:,iCavNode)= P.Node.q(:,iCavNode)+Yp; % short circuit flow of node
P.Node.p= P.Node.q./P.Node.Y; % node pressure by solving NodeFlow q=0

% flow into all Cavities, connected to nodes
P.Cavity.VDot= (P.Node.p(:,iCavNode)-P.Cavity.p)./P.Cavity.Z;
end

function X=SparseDiag(x)
X=sparse(diag(x(:)));
end

