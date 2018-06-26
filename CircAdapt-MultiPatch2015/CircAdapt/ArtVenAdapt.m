function ArtVenAdapt(StrAV,AdaptType)
% function ArtVenAdapt(StrAV,AdaptType);
% Adaptation of Diameter and Wall thickness of Art and Ven to
% pressure and flow.
% StrAV= array of ArtVen names, e.g. {'Sy','Pu'}
% AdaptType= {'Diameter', 'WallVolume'} indicates type of adaptation
% Theo Arts, Maastricht University, Oct 30, 2011

global P

% Determine ArtVen indices
iAV   = Str2Index(StrAV ,P.ArtVen.Name); % read string -> index
nAV   = size(iAV,2);

% read type(s) of adaptation
iAdapt= Str2Index(AdaptType,{'Diameter','WallVolume'});
Adapt=[0,0]; Adapt(iAdapt)=1;
AdaptDiameter=Adapt(1); AdaptWallVolume=Adapt(2); % active if Adapt==1

szAV=[2,nAV]; szRow=[1,2*nAV]; % shape of arry in P.ArtVen and P.Cavity
iCavArt = P.ArtVen.iCavity(:,iAV); iCavVen=iCavArt+1;
iCav    = reshape([iCavArt;iCavVen],szRow); % corresponding Cavity indices in row
q    = P.ArtVen.q(:,reshape([iAV;iAV],szRow)); % Microcirculatory flow
A    = P.Cavity.A(:,iCav); % cavity cross-section
p    = P.Cavity.p(:,iCav); % cavity pressure
Z    = P.Cavity.Z(:,iCav); % cavity wave impedance
vM   = reshape(mean(q./A),szAV); % mean flow velocity
AMean= reshape(mean(A),szAV);
pMean= reshape( ZerosFind(p,log(A/SparseDiag(AMean(:)))),szAV ) ; % pMean=p(AMean)

AWall  = P.ArtVen.AWall(:,iAV); % wall cross-section
p0     = P.ArtVen.p0(:,iAV); % working pressure
A0     = P.ArtVen.A0(:,iAV); % working cress-section
vImpact   = P.ArtVen.Adapt.vImpact(:,iAV); % assumed body impact velocity
vFlowMean = P.ArtVen.Adapt.vFlowMean(:,iAV); % target value flow velocity

pMax= max(p) + max((Z.*A)*SparseDiag(vImpact)); % maximum pressure
AMax= max(A);
pMax      =reshape(pMax,szAV);
AMax      =reshape(AMax,szAV);
WallStress= 3*pMax.*(0.5+AMax./AWall); % maximum wall stress

Fac_vFlow    = vM./vFlowMean; % mean velocity/target value
FacWallStress= WallStress./P.ArtVen.Adapt.WallStress(:,iAV); % max wall stress/target value
Facp0        = pMean./p0; % mean pressure/target value
FacA0        = AMean./A0; % mean cross-section/target value
%---

%=== Carrying out adaptation
FacV=ones(size(A0));
if AdaptDiameter; % adapts diameter to flow
    a= 0.25; % feedback factor, low value-> slow but stabile
    FacV = reshape(Fac_vFlow.^a,szRow); % volume stabilization during adaptation
    A0   = A0 .* (FacA0.*Fac_vFlow).^a; % set working cross-section
    p0   = p0 .* Facp0.^a ; % set working pressure
end

if AdaptWallVolume % adapts wall thickness to pressure
    a=0.5; % feedback factor, low value-> slow but stabile
    AWall= AWall .* FacWallStress.^a;
end
%---
P.Cavity.V(end,iCav)= P.Cavity.V(end,iCav).*FacV(:)'; 
% SVar stabilization
P.ArtVen.A0(:,iAV)   =A0; % working cross-section
P.ArtVen.p0(:,iAV)   =p0; % working pressure
P.ArtVen.AWall(:,iAV)=AWall; % wall cross-section

disp('relative deviation of:')
disp(['vFlow     : ', num2str(round(1000*log(Fac_vFlow(:)'     )),'%+6.3d')]);
disp(['WallStress: ', num2str(round(1000*log(FacWallStress(:)')),'%+6.3d')]);

end

function Index=Str2Index(Str,Names)
% conversion of Name(s) to indices in structure
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

function pZero= ZerosFind(Y,X)
% X and Y(X) are samples points, ordered in columns
% find Y(X) for X=0 for each column by polinomal interpolation
sz=size(X);
pZero=zeros(1,sz(2));
Col1= ones(sz(1),1);
for i=1:sz(2)
    x= X(:,i);
    y= Y(:,i);
    Coef=pinv([Col1,x,x.^2,x.^3])*y; % 3rd order polinomal
    pZero(i)=Coef(1);
end
end

