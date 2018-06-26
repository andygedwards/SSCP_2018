function ArtVenV2p
%function ArtVenV2p
% Elastic Art/Ven hemodynamics
% Flow waves enter Art and Ven. Art and Ven are connected by peripheral
% resistance, representing the microcirculation.
% volume V-> transmural pressure pTrans(V) and wave impedance Z(V)
% Theo Arts, Maastricht University, Oct 30, 2011

global P;

ArtVen = P.ArtVen; % ArtVen structure
rhob   = P.General.rhob;
% indices referring to realted cavities and walls
iC= [P.ArtVen.iCavity; P.ArtVen.iCavity+1]; iCavity=iC(:);
iW= [P.ArtVen.iWall  ; P.ArtVen.iWall+1]  ; iWall  =iW(:);

V      = P.Cavity.V(:,iCavity); % cavity volumes (state variable)
% SparseDiag function is used for efficient column-wise multiplication
Len  = SparseDiag(ArtVen.Len  ); % repesentative length of blood vessels
AWall= SparseDiag(ArtVen.AWall); % wall cross-section of blood vessel
kd3m1= SparseDiag(ArtVen.k/3-1); % stiffness exponential
p0   = SparseDiag(ArtVen.p0   ); % working pressure
A0   = SparseDiag(ArtVen.A0   ); % working cross-section
Half = SparseDiag([.5,.5,.5,.5]); % 0.5 diagonal matrix

A      = max(1e-20,V/Len); % vessel cross-section, buckling included
ANorm  = A/AWall; % cross-section normalized to wall volume
p      = exp( log( max(1e-20, ANorm+0.5)/(A0/AWall+Half ) ) * kd3m1) * p0; %p(A)
Z      = sqrt( (p./(ANorm.*(ANorm+0.5 )) * (rhob*kd3m1/AWall.^2) ) ) ; %Z(A)

P.Wall.pTrans(:,iWall)=p; % transmural pressure
P.Cavity.Z(:,iCavity) =Z; % wave impedance
P.Cavity.A(:,iCavity) =A; % cross-section

end

% function X=SparseDiag(x)
% X=sparse(diag(x(:)));
% end

