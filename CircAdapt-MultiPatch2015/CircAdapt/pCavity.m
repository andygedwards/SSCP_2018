function pCavity
% function pCavity
% pericardium encloses a number of walls and cavities. Pressures are
% calculated by adding transmural pressures
% Theo Arts, Maastricht University, Oct 30, 2011

global P

% Pericardium
P.Bag.p=P.Wall.pTrans(:,P.Bag.iWall);
%TriSeg +Pericard pressure
iCavity=P.TriSeg.iCavity;
iWall  =P.TriSeg.iWall;
P.Cavity.p(:,iCavity)  =P.Bag.p-P.Wall.pTrans(:,iWall); %pLV
P.Cavity.p(:,iCavity+1)=P.Bag.p+P.Wall.pTrans(:,iWall+2); %pRV
%Chamber (La, Ra)+Pericard pressure
iCavity=P.Chamber.iCavity; Row0=zeros(1,P.Chamber.n);
iWall  =P.Chamber.iWall;
P.Cavity.p(:,iCavity)=P.Bag.p(:,1+Row0)+P.Wall.pTrans(:,iWall);
%ArtVen, not enclosed by other walls
iC=P.ArtVen.iCavity; iCavity=[iC;iC+1];
iW=P.ArtVen.iWall  ; iWall=[iW;iW+1];
P.Cavity.p(:,iCavity)=P.Wall.pTrans(:,iWall);

end

