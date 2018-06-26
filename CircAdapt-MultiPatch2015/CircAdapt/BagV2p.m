function BagV2p
% function BagV2p
% A passive elastic bag like the pricardiam.
% Pericard transmural pressure P.Wall.pTrans = fu(enclosed volume)
% No cavity, only a wall
% k= non-linearity exponential
% Theo Arts, Maastricht University, Oct 30, 2011

global P;

% works for a single Bag
% pericardium enclosed cardiac wall and cavities
iCav =[P.TriSeg.iCavity+(0:1),P.Chamber.iCavity]; % enclosed cavities
iWall=[P.TriSeg.iWall+(0:2),P.Chamber.iWall]; % enclosed walls
V= sum(P.Cavity.V(:,iCav),2)+sum(P.Wall.VWall(iWall)); % total enclosed volume
VNorm   = V/P.Bag.VRef; % normalized volume
P.Wall.pTrans(:,P.Bag.iWall)= P.Bag.pAdapt * VNorm .^ P.Bag.k; 

end

